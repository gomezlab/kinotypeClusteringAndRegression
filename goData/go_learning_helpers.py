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

def get_tts(dat, X_col_name='GO Labels', test_size=0.2, random_state=None):
    ''' Gets the train_test_split on dat
    inputs: dat -- DataFrame of kinases with cluster number and GO Labels attached, as output by add_cluster_labels
            X_col_name -- name of the column(s) to use for the predictor variables
            test_size -- portion of data to place in the test set
            random_state -- random seed used for train_test_split
    '''
    # split into training & testing
    X_train, X_test, y_train, y_test = tts(dat[X_col_name], dat['cluster'], test_size=test_size, random_state=random_state)

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
            kwargs -- a kwarg dict needed/for customizability for the scorers
    returns: out -- either a scalar score or an iterable/array with detailed scoring information
    '''

    # this stores the scorers
    score_dict = {'accuracy':skm.accuracy_score, 'matthews':skm.matthews_corrcoef,
                'cohen':skm.cohen_kappa_score, 'confusion':skm.confusion_matrix,
                'report':skm.classification_report, 'f1':skm.f1_score,
                'hinge':skm.hinge_loss}

    scorer = score_dict[metric]

    score = scorer(y_true, y_predict, **kwargs)

    return score

def validate_learnability(n_run, dat_df, clf, X_col_name='GO Labels', test_size = 0.3, parallel=None, metrics=['accuracy'], scorer_kwargs={}):
    ''' takes the output of an add_cluster_labels and calls clf.fit(), get_tts, and score_model n_run times
    inputs: dat_df -- the input dat with 'cluster' response and X_col_name predictor
            clf -- a classifier with .fit() and .predict() methods
            X_col_name -- a string (or list of strings) corresponding to the predictors in dat_df
            test_size -- proportion of dat_df to use in the holdout sets
            parallel -- how many cores to use, default None is sequential operation on a single core
            metric -- the metric to call for score_model
    output: out -- an array of length n_run of scores from score_model on each of the holdout sets
    '''
    num_metrics = len(metrics)

    if parallel is None:
        clf_local = clone(clf) # clone original clf for parallel read/write safety
        output = []
        for i in range(n_run):
            # split, fit, score, append
            X_train, X_test, y_train, y_test = get_tts(dat_df, X_col_name=X_col_name, test_size=test_size)
            clf.fit(X_train, y_train)

            # loop through the metrics
            scores = []
            for i in range(num_metrics):
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
        kwarg_dict = {'dat_df':dat_df, 'clf':clf, 'X_col_name':X_col_name, 'test_size':test_size, 'parallel':None, 'metrics':metrics, 'scorer_kwargs':scorer_kwargs}
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
    idx_locs = [x for x, y in enumerate(np.array(list(it.combinations(range(num_classes), 2)))) if cluster_num in y]
    return svm.coef_[idx_locs,:]