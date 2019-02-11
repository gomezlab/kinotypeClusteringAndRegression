import pandas as pd
import numpy as np
from copy import copy

def fetch_hgnc_mapper(path_to_hgnc='../data/hgnc_alias_list.txt',
                     path_to_kmast='../data/KINASESmasterlist_w_Aliases.xlsx'):
    hgnc = pd.read_csv(path_to_hgnc, sep='\t')

    hgnc = hgnc[hgnc['Approved symbol'].apply(lambda x: 'withdrawn' not in x)]

    # get the original keys
    hgnc_original_keys = hgnc['Approved symbol'].unique()

    # drop columns for efficiency
    hgnc = hgnc[list(hgnc.columns)[0:6]].drop('Status', axis=1)

    # filter out Nan synonyms (not helpful)
    hgnc_syn_list = hgnc[~ hgnc.Synonyms.isna()]
    hgnc_prev_symb_list = hgnc[~ hgnc['Previous symbols'].isna()]

    #convert the synonyms column to a list
    # convert these lists to pd.Series
    # merge with original dataframe
    #drop old synonyms column 
    # melt the new columns into rows

    current_syn_list = hgnc_syn_list.Synonyms.apply(lambda x: x.split(',')) \
        .apply(pd.Series) \
        .merge(hgnc, left_index = True, right_index = True) \
        .drop(["Synonyms"], axis = 1) \
        .melt(id_vars = ['HGNC ID', 'Approved symbol', 'Approved name', 'Previous symbols'], value_name = "synonym") 

    current_syn_list = current_syn_list[~ current_syn_list.synonym.isna()]
    current_syn_list.synonym = current_syn_list.synonym.apply(lambda x: x.replace(' ',''))

    prev_symb_list = hgnc_prev_symb_list['Previous symbols'].apply(lambda x: x.split(',')) \
        .apply(pd.Series) \
        .merge(hgnc, left_index = True, right_index = True) \
        .drop(['Previous symbols'], axis = 1) \
        .melt(id_vars = ['HGNC ID', 'Approved symbol', 'Approved name', 'Synonyms'], value_name = "synonym") 

    prev_symb_list = prev_symb_list[~ prev_symb_list.synonym.isna()]
    prev_symb_list.synonym = prev_symb_list.synonym.apply(lambda x: x.replace(' ',''))
    
    hgnc_mapper = dict(zip(current_syn_list['synonym'], current_syn_list['Approved symbol']))
    # add in HGNC ID mapper
    hgnc_mapper.update(dict(zip(current_syn_list['HGNC ID'], current_syn_list['Approved symbol'])))
    hgnc_mapper_previous = dict(zip(prev_symb_list['synonym'], prev_symb_list['Approved symbol']))

    trouble_list = list(filter(lambda x: hgnc_mapper[x] != hgnc_mapper_previous[x], set(hgnc_mapper.keys())&set(hgnc_mapper_previous.keys())))

    hand_coded = {'RAGE':'MOK', 'SGK2':'SGK2', 'SGK196':'SGK196', 'MAPK3':'MAPK3'}

    hgnc_mapper_previous.update(hgnc_mapper) #overwrite the previous symbol conflicts

    hgnc_mapper = hgnc_mapper_previous
    hgnc_mapper.update({x:x for x in hgnc_original_keys}) #keep the identify maps
    hgnc_mapper.update(hand_coded) # overwrite the trouble list
    
    ### leverage additional known relations from kmast
    kmast = pd.read_excel(path_to_kmast)
    for uni,ms,rna in zip(kmast['Uniprot Protein'], kmast['MS Gene'], kmast['RNAseq Gene']):
        if uni not in hgnc_mapper.keys():
            if ms in hgnc_mapper.keys():
                hgnc_mapper[uni] = hgnc_mapper[ms]
            elif rna in hgnc_mapper.keys():
                hgnc_mapper[uni] = hgnc_mapper[rna]
        if ms not in hgnc_mapper.keys():
            if uni in hgnc_mapper.keys():
                hgnc_mapper[ms] = hgnc_mapper[uni]
            elif rna in hgnc_mapper.keys():
                hgnc_mapper[ms] = hgnc_mapper[rna]
        if rna not in hgnc_mapper.keys():
            if ms in hgnc_mapper.keys():
                hgnc_mapper[rna] = hgnc_mapper[ms]
            elif uni in hgnc_mapper.keys():
                hgnc_mapper[rna] = hgnc_mapper[uni]

    hgnc_mapper = {x.upper():y.upper() for x,y in hgnc_mapper.items() if x is not np.nan and y is not np.nan}
    
    hgnc_mapper = {x.upper():y.upper() for x,y in hgnc_mapper.items() if x is not np.nan and y is not np.nan}
    
    return hgnc_mapper