import pandas as pd
import numpy as np
import itertools as it
import functools as ft
import multiprocessing as mp
from time import time

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.ensemble import BaggingClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model.coordinate_descent import LinearModelCV
from sklearn.linear_model import ElasticNetCV
from sklearn.base import RegressorMixin
from sklearn.linear_model.base import LinearModel
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import cross_validate
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import *

# control parameter for demos
SUBSET=False
CV_JOBS = 24
N_JOBS = 1

# import the data file
agg_data = pd.read_csv('~/Github/subnetGRcurves/data/agg_growth.csv.gz', index_col=0, compression='gzip')

### Slice relevant data ###
conc = agg_data.Concentration
drugs = agg_data.Drug
targ = agg_data.Target
y = agg_data.Growth
X = agg_data[agg_data.columns[4:]]

# set parameters
clf = None
clf_loss = ['modified_huber']
threshes=[5e-2]
l1_ratios=[0, 0.1, 0.5, 0.9, 0.99]
clf_alphas = [3e-4]
clf_class_weights = [None, 'balanced', {0:1, 1:2}, {0:1, 1:6}, {0:1, 1:12}]

if SUBSET:
    clf_alphas = [1e-4]
    threshes=[1e-12]
    l1_ratios=[0]
    clf_class_weights = [None, 'balanced']

def thresholder(thresh, arr):
    return (arr > 1- thresh) * 1.

clf_transformer = ft.partial(thresholder, 1e-12)

# cross validator
logo = LeaveOneGroupOut()

# we can't grid search across threshes, so store a dict of them
thresh_grids = {}
param_grid = {
    'base_estimator__loss' : clf_loss,
    'base_estimator__class_weight' : clf_class_weights,
    'base_estimator__l1_ratio' : l1_ratios,
    'base_estimator__alpha' : clf_alphas,
}

scorer_dict = {'ap':average_precision_score,
              'roc_auc':roc_auc_score,
              'accuracy':accuracy_score,
              'f1':f1_score,
              'pearsonr':lambda x,y: pearsonr(x,y)[0]}

# convert to scorers
scorer_dict = {x:make_scorer(y) for x,y in scorer_dict.items()}

# define the threshes
for thresh in threshes:
    base_clf = BaggingClassifier(base_estimator=SGDClassifier(penalty='elasticnet',
                                                              learning_rate='optimal',
                                                              random_state=1920,
                                                              tol=1e-3),
                                 n_estimators=30,
                                 max_samples=0.632,
                                 n_jobs=N_JOBS)

    thresh_grids[thresh] = GridSearchCV(base_clf,
                                        cv=logo.split(X, y, groups=drugs),
                                        param_grid=param_grid,
                                        pre_dispatch='2*n_jobs',
                                        n_jobs=CV_JOBS,
                                        iid=False,
                                        refit=False,
                                        return_train_score=True)

# fit the models
for thresh, gCV in thresh_grids.items():
    clf_transformer = ft.partial(thresholder, thresh)
    y_bin = clf_transformer(y)

    gCV.fit(X, y_bin)

# save results
for thresh, gCV in thresh_grids.items():
    dat = pd.DataFrame(gCV.cv_results_)
    dat.to_csv('~/Github/subnetGRcurves/results/gCV_' + str(thresh) + 'b.csv')
