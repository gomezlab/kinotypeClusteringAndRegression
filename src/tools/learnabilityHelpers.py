import pandas as pd
import numpy as np
import functools as ft
import itertools as it
import matplotlib.pyplot as plt
import sklearn
from sklearn.preprocessing import LabelBinarizer as LB
from sklearn.model_selection import train_test_split as tts
import sklearn.metrics as skm
from sklearn.base import clone
import multiprocessing as mp
from copy import copy
from gensim import corpora
from collections import defaultdict
from gensim.models import TfidfModel
from collections import Counter
import scipy.sparse as sp


def get_go_annotations_series(path_to_file):
    ''' Reads a csv file of {Kinase: "[Annotations]"} and returns a formatted pd.Series
    input: path_to_file -- "./example_folder/example.csv"
    returns: pd.Series with Kinases as index, list of tokens/phrases as outputs
    '''
    dat = pd.read_csv(path_to_file, header=None) # read in the data
    dat.columns = ['Kinase', 'GO Labels']
    dat.set_index('Kinase', inplace=True)
    dat = dat['GO Labels'] # convert to pd.Series

    # csv saves a list as a string with a list inside
    # this is a poor but working way to turn that string back to a list
    def csv_string_to_list(i):
        return list(map(lambda x: x.lstrip(), i.replace('[', '').replace(']', '').replace('\'', '').split(',')))


    dat = dat.apply(csv_string_to_list)
    return dat

def convert_go_annotations_to_one_hot(ser, return_binarizer=True):
    ''' Takes a series object and returns the one-hot encoded version as well as the LabelBinarizer object
    inputs: ser -- series of {Kinase: [Annotations]}
            return_binarizer -- whether or not the binarizer used is returned
    outputs: lb -- a fitted LabelBinarizer
            ser2 -- a pd.Series of {Kinase: [BinarizedVectors]}
    '''

    # create a label LabelBinarizer object to binarize
    lb = LB()

    # fit the binarizer to unique token/phrases in the sentences (lists) of ser
    lb.fit(list(set(filter(lambda x: x is not '', [token for sentence in ser for token in sentence]))))

    # sum converts a list of [00000100000, 00010000000, 00000000001]  vectors into 00100100001 vectors
    ser2 = ser.apply(lb.transform).apply(sum)

    # return results as specified by return_binarizer
    if(return_binarizer):
        return lb, ser2
    else:
        return ser2



def add_cluster_labels(path_to_file, ser):
    ''' Reads a csv file (that was probably written from R) into a series of values joined with the series data input
    inputs: path_to_file -- "./example_folder/example.csv"
            ser -- a series of GO data, as output from convert_go_annotations_to_one_hot
    outputs: pd.Series with Kinases as index, list of tokens/phrases as outputs
    '''

    dat = pd.read_csv(path_to_file, header=0, sep='\t')
    dat.set_index('names', inplace=True)
    n, m = dat.shape # save shape for drop amount
    dat = dat.join(other=ser) # SQL left join
    dat.dropna(axis=0, inplace=True)
    
    print('Dropped ', n-dat.shape[0], 'kinases due to zero length post-processing')

    # no need to verify size of output, left join is always the same as input.
    # may be good to add a NaN checker later
    return dat

def get_tts(dat, X_col_name='GO Labels', Y_col_name='cluster', test_size=0.2, random_state=None):
    ''' Gets the train_test_split on dat
    inputs: dat -- DataFrame of kinases with cluster number and GO Labels attached, as output by add_cluster_labels
            X_col_name -- name of the column(s) to use for the predictor variables
            Y_col_name -- name of the column to use for the response variable
            test_size -- portion of data to place in the test set
            random_state -- random seed used for train_test_split
    '''
    # split into training & testing
    X_train, X_test, y_train, y_test = tts(dat[X_col_name], dat[Y_col_name], test_size=test_size, random_state=random_state)

    # convert to sklearn-friendly format (lists of arrays)
    X_train = X_train.values.tolist()
    X_test = X_test.values.tolist()
    y_train = y_train.values.tolist()
    y_test = y_test.values.tolist()

    return (X_train, X_test, y_train, y_test)

def score_model(y_true, y_predict, metric='accuracy', kwargs={}):
    ''' Scores a model according to the desired loss function
    inputs: y_predict -- predicted labels
            y_true -- the true labels
            metric -- this is the fun part. There are many good ways to score, but each has its bias/upside!
                        one of: 'accuracy' - use accuracy_score (default, can also be weighted)
                                'matthews' - use matthews_corrcoef (good for unbalanced)
                                'cohen' - use cohen_kappa_score (can be weighted with empirical probabilities)
                                'confusion' - just give the whole confusion_matrix -- might take a long while
                                'report' - ditto, but whole classification_report -- might take a little while
                                'f1' - f1_score (can also be weighted/averaged)
                                'hinge' - hinge_loss (can also be weighted/averaged)
                                'prf' - precision_recall_fscore_support (per class)
            kwargs -- a kwarg dict needed/for customizability for the scorers
    returns: out -- either a scalar score or an iterable/array with detailed scoring information
    '''

    # this stores the scorers
    score_dict = {'accuracy':skm.accuracy_score, 'matthews':skm.matthews_corrcoef,
                'cohen':skm.cohen_kappa_score, 'confusion':skm.confusion_matrix,
                'report':skm.classification_report, 'f1':skm.f1_score,
                'hinge':skm.hinge_loss, 'prf':skm.precision_recall_fscore_support}

    scorer = score_dict[metric]

    score = scorer(y_true, y_predict, **kwargs)

    return score

def validate_learnability(n_run, dat_df, clf, X_col_name='GO Labels', Y_col_name='cluster', test_size = 0.3, parallel=None, metrics=['accuracy'], scorer_kwargs={}):
    ''' takes the output of an add_cluster_labels and calls clf.fit(), get_tts, and score_model n_run times
    inputs: dat_df -- the input dat with 'cluster' response and X_col_name predictor
            clf -- a classifier with .fit() and .predict() methods
            X_col_name -- a string (or list of strings) corresponding to the predictors in dat_df
            Y_col_name -- string for the column of the response variable in dat_df (default = 'cluster')
            test_size -- proportion of dat_df to use in the holdout sets
            parallel -- how many cores to use, default None is sequential operation on a single core
            metric -- the metric to call for score_model
    output: out -- an array of length n_run of scores from score_model on each of the holdout sets
    '''
    num_metrics = len(metrics)

    if parallel is None:
        clf_local = clone(clf) # clone original clf for parallel read/write safety
        output = []
        
        # make sure we normalize the confusion matrix for any missing classes
        if 'confusion' in metrics or 'prf' in metrics:
            # find out & sort the classes
            classes = sorted(list(set(dat_df[Y_col_name])))

            # copy scorer kwargs and add the classes
            confusion_scorer_kwargs = copy(scorer_kwargs)
            confusion_scorer_kwargs.update({'labels':classes})
        
        for i in range(n_run):
            # split, fit, score, append
            X_train, X_test, y_train, y_test = get_tts(dat_df, 
                                                       X_col_name=X_col_name, 
                                                       Y_col_name=Y_col_name, 
                                                       test_size=test_size, 
                                                       random_state = int(id(mp.current_process())+i)%2**32)
            
            # it's possible that only one class will show up in the y_train
            # this will throw an error, so we re-sample in case this happens
            rerun=2
            while len(np.unique(y_train)) < 1:
                X_train, X_test, y_train, y_test = get_tts(dat_df, 
                                                       X_col_name=X_col_name, 
                                                       Y_col_name=Y_col_name, 
                                                       test_size=test_size, 
                                                       random_state = int(id(mp.current_process())*rerun+i)%2**32)
                rerun+=1
            
            # fit the model
            clf.fit(X_train, y_train)

            # loop through the metrics
            scores = []
            for i in range(num_metrics):
                # check if current metric is confusion
                if(metrics[i]) == 'confusion' or (metrics[i]) == 'prf':
                    scores.append(score_model(y_test, clf.predict(X_test), metric=metrics[i], kwargs=confusion_scorer_kwargs))
                else:
                    scores.append(score_model(y_test, clf.predict(X_test), metric=metrics[i], kwargs=scorer_kwargs))

            # append a tuple if more than one metric, else just a float
            if(num_metrics > 1):
                output.append(tuple(scores))
            else:
                output.append(scores[0])
    else:
        n_run_array = [n_run//parallel]*parallel # number of runs per core

        # add one if needed to some cores
        overflow = n_run % parallel
        if overflow > 0:
            n_run_array[:overflow] = n_run_array[:overflow] + [1]*overflow

        # create a handy kwarg dict --- note that parallel is None, so we run sequentially in parallel
        kwarg_dict = {'dat_df':dat_df, 'clf':clf, 'X_col_name':X_col_name, 'Y_col_name':Y_col_name, 'test_size':test_size, 'parallel':None, 'metrics':metrics, 'scorer_kwargs':scorer_kwargs}
        val_learn = ft.partial(validate_learnability, **kwarg_dict)

        # create a pool, map the run values
        p = mp.Pool(parallel)
        output = p.map(val_learn, n_run_array)
        p.close()

        output = list(it.chain.from_iterable(output))

    return output
    
def get_svm_coeffs_for_cluster(svm, cluster_num):
    num_classes = len(svm.classes_)
    combos = list(it.combinations(range(num_classes),2))
    
    # append positive and negative occurrences 
    idx_locs_positive = [x for x, y in enumerate(np.array(list(it.combinations(range(num_classes), 2)))) if cluster_num == y[0]]
    idx_locs_negative = [x for x, y in enumerate(np.array(list(it.combinations(range(num_classes), 2)))) if cluster_num == y[1]]
    
    return np.vstack([svm.coef_[idx_locs_positive,:], -svm.coef_[idx_locs_negative,:]])

def get_cnb_coeffs_for_cluster(cnb, cluster_num):
    num_classes = len(t.class_count_)
    probas = np.exp(- cnb.feature_log_prob_)
    
    idx_locs_positive = [i for i in range(num_classes) if i != cluster_num]
    idx_locs_negative = [cluster_num]
    
    return np.vstack([probas[idx_locs_positive,:], -probas[idx_locs_negative,:]])

def generate_kinase_labels(path_to_synonyms='../data/goData/go_synonym_data.txt',
                           path_to_kinase_network='../data/kin/kin_anscombe_weighted.csv', 
                           path_to_alias_spreadsheet='../data/ref/KINASESmasterlist_w_Aliases.xlsx', 
                           path_to_stopwords='../data/goData/stopwords.csv', 
                           path_to_process_list='../data/goData/go_biological_processes.txt',
                           use_go_subset=False,
                           path_to_go_subset='../data/goData/go_subset.csv',
                           out_path = '../data/goData/kinase_go_processes.csv'):
    '''
    Generates a labeled csv file of {Kinase: [Annotations]} using the inputs provided.
    '''
    df = pd.read_csv(path_to_synonyms, header=None, sep='\t', low_memory=False,)
    df.columns = ['ID', 'Gene/Product', 'Name', 'GO Class Labels', 'Synonyms']

    def splitter(x):
        try:
            temp = x.split('|')
        except:
            temp = []
        return temp

    df['Synonyms'] = df['Synonyms'].apply(splitter)

    kinase_network_df = pd.read_csv(path_to_kinase_network, header=None, sep='\t')
    known_kinases = list(set(kinase_network_df[0]) | set(kinase_network_df[1]))

    alias = pd.read_excel(path_to_alias_spreadsheet, header = 0)
    kin_map = alias.set_index('Uniprot Protein')['MS Gene'].to_dict()
    all_aliases = alias.set_index('Uniprot Protein')['Aliases (Conservative)'].dropna().apply(lambda x: x.split(',')).to_dict()

    go_dat = {}

    for k in known_kinases:
        temp = df[df['Gene/Product']==k]
        if(temp.shape[0] == 0):
            temp = df[df['Synonyms'].apply(lambda x: k in x)]
            if(temp.shape[0] == 0):
                temp = df[df['Gene/Product']==kin_map[k]]
                if(temp.shape[0] == 0):
                    temp = df[df['Synonyms'].apply(lambda x: kin_map[k] in x)]
                    if(temp.shape[0]==0):
                        r = all_aliases.get(k)
                        if(r is not None):
                            for a in r:
                                temp = df[df['Gene/Product']==a]
                                if(temp.shape[0] ==0):
                                    temp = df[df['Synonyms'].apply(lambda x: a in x)]
                                    if(temp.shape[0]==0):
                                        pass
                                    else:
                                        go_dat[k] = temp
                                        break
                                else:
                                    go_dat[k] = temp
                    else:
                        go_dat[k] = temp
                else:
                    go_dat[k] = temp
            else:
                go_dat[k] = temp
        else:
            go_dat[k] = temp

    # find any kinases where multiple gene/product IDs were returned
    fix_list = []

    for x in go_dat.keys():
        if(go_dat[x].shape[0] > 1):
            fix_list.append(x)

    for x in fix_list:
        temp = go_dat[x].iloc[0]
        temp['GO Class Labels']='|'.join(list(set(go_dat[x]['GO Class Labels'].iloc[0].split('|')) | set(go_dat[x]['GO Class Labels'].iloc[1].split('|'))))
        go_dat[x] = temp

    def helper(x):
        try:
            temp = go_dat[x]['GO Class Labels'].values[0].split('|')
        except:
            temp = go_dat[x]['GO Class Labels'].split('|')
        return temp

    just_labels = {x:helper(x) for x in go_dat.keys()}
    agg_labels = [x for y in just_labels.values() for x in y]

    stopwords = pd.read_csv(path_to_stopwords).iloc[:,0].tolist()
    processes = pd.read_csv(path_to_process_list, sep='\t', header=None).set_index(0)
    
    # apply subset filter if present
    if use_go_subset:
        go_subset = set(pd.read_csv(path_to_go_subset, sep ='\t', header=None)[0])
        indices_to_keep = pd.Series(processes.index.to_list()).isin(go_subset).to_list()
        processes = processes[indices_to_keep]
         
    processes = set(processes[1].tolist())

    def stophelper_plus_is_process(x):
        try:
            temp = go_dat[x]['GO Class Labels'].values[0].split('|')
        except:
            temp = go_dat[x]['GO Class Labels'].split('|')
        return list(filter(lambda x: x in processes, filter(lambda x: x not in stopwords, temp)))

    kinase_labels = {x:stophelper_plus_is_process(x) for x in go_dat.keys()}
    labeled_kinases = pd.Series(kinase_labels)

    labeled_kinases.to_csv(out_path)
                           
    return None

def get_go_data_dict(go_file='../data/goData/kinase_go_processes.csv', go_dir='./'):
    # get the go annotations
    go_annotations = get_go_annotations_series(go_dir+go_file)

    # filter out any empty kinases
    go_annotations = go_annotations[map(lambda x: len(x) > 0, [x for x in go_annotations])]

    # generate binarized go_annotations
    go_lab_binner, bin_go_annotations = convert_go_annotations_to_one_hot(go_annotations,
                                                       return_binarizer= True)

    # create a frequency chart for every label in the go_annotations 
    # this allows us to screen out labels with only one occurence (to help reduce noise)
    frequency = defaultdict(int)
    for doc in go_annotations:
        for word in doc:
            frequency[word] += 1
    texts = [[word for word in doc if frequency[word] > 1]
             for doc in go_annotations]

    # generate binarized go_annotations for the intermediate, filtered go_annotations
    filter_go_annotations = pd.Series(texts)
    filter_go_annotations.index = go_annotations.index
    filter_go_annotations.name = go_annotations.name

    # filter out any empty filter_kinases
    filter_go_annotations = filter_go_annotations[map(lambda x: len(x) > 0, [x for x in filter_go_annotations])]

    freq_go_lab_binner, freq_go_annotations = convert_go_annotations_to_one_hot(filter_go_annotations,
                                                       return_binarizer= True)

    # create a corpus and dictionary (reference) objects for a tfidf model for the go_annotations
    dictionary = corpora.Dictionary(texts)
    corpus = [dictionary.doc2bow(doc) for doc in texts]

    # create a tfidf series for the go_annotations
    tfidf = TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]

    # we have to parse gensim's sparse format back into pd.Series format
    hold = []
    rej = []
    for i in range(len(corpus_tfidf)):
        if(len(corpus_tfidf[i]) > 1):
            a, b = map(np.asarray, zip(*corpus_tfidf[i]))
            row = np.repeat(0, len(a))
            hold += [sp.csc_matrix((b, (row, a)), shape=(1,len(dictionary)))]
        else:
            rej += [go_annotations.index[i]]


    tfidf_go_annotations = pd.Series([x for x in sp.vstack(hold).toarray()])
    tfidf_go_annotations.index = [x for x in go_annotations.index.tolist() if x not in rej]
    tfidf_go_annotations.name = go_annotations.name

    # create a dictionary with our datasets
    go_data_dict = {'go':bin_go_annotations,'freq':freq_go_annotations,'tfidf':tfidf_go_annotations}
    transform_dict = {'freq':freq_go_lab_binner}
    
    return go_data_dict, transform_dict