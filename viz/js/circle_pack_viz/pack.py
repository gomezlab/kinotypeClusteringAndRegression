import sys
import pandas as pd
import numpy as np
import functools as ft
import itertools as it
import scipy.sparse as sp
import itertools
from collections import defaultdict
from collections import Counter
import json

##### ##### ##### ##### #####
# Read in data
##### ##### ##### ##### #####

dat = pd.read_csv(sys.argv[1], index_col=0, sep='\t')

##### ##### ##### ##### #####
# Parse the dataframe for individual items
##### ##### ##### ##### #####

# get a list of the kinases
kin_list = list(dat.index)

# define the sizes
sizes = np.array(np.int(np.log(s + 1)) for s in dat['size'])
sizes = dict(zip(kin_list, sizes.tolist()))

# sizes = {n:int(10) for n in louvain_small_clusters['names'].unique()} # example

# handy to have a kin_arr
kin_arr = np.array(kin_list)

kin_parents = np.array(dat['super'].tolist())
kin_children = np.array(dat['sub'].tolist())

parents = {par:set(kin_children[kin_parents==par].tolist()) for par in np.unique(kin_parents).tolist()}

# get parents plus children for parents != children
parents_plus_children = {x:{z:kin_arr[kin_children==z].tolist() for z in y if(list(y)[0] != x)} for x,y in parents.items()}

# identify parents == children
# missing_parents = [x for x,y in parents_plus_children.items() if len(y)<1]

# update to add kinases to remove parents == children
# parents_plus_children = {x:y for x,y in parents_plus_children.items() if x not in missing_parents}

# place parents == children into separate dict
# missing_parents = {x:kin_arr[kin_parents==x].tolist() for x in missing_parents}

opacity = False

# check if we are tracking opacity
if len(np.unique(dat['opacity'])) > 1:
    opacity = True
    opacity_dict = dict(zip(kin_list, dat['opacity'].tolist()))

dark_dict = dict(zip(kin_list, dat['color']))

##### ##### ##### ##### #####
# Write the data to a json
##### ##### ##### ##### #####

# first get the kinases of depth 3
if not opacity:
    json_out = {"name":"viz",
                "children":[
                    {"name":str(parent), 
                    "children":[
                        {"name":str(child), 
                        "children":[
                            {"name":str(k),
                            "size":int(sizes[k]),
                            "studied":dark_dict[k],
                            "opacity":1}
                            for k in kinases]}
                         for child, kinases in children.items()]} 
                    for parent, children in parents_plus_children.items()]}
# opacity support 
else:
    json_out = {"name":"viz",
                "children":[
                    {"name":str(parent), 
                    "children":[
                        {"name":str(child), 
                        "children":[
                            {"name":str(k),
                            "size":int(sizes[k]),
                            "studied":dark_dict[k],
                            "opacity":opacity_dict[k]}
                            for k in kinases]}
                         for child, kinases in children.items()]} 
                    for parent, children in parents_plus_children.items()]}

# extra json
# missing_json = [
#     {"name":str(parent),
#     "children":[
#         {"name":str(k),
#         "size":int(sizes[k])}
#         for k in kinases]}
#     for parent, kinases in missing_parents.items()]

# json_out["children"] = json_out["children"] + missing_json

str_out = json.dumps(json_out).replace("},", "}, \n").replace('[{', '[\n{\n').replace(']},', ']}\n,\n')

with open('dist/viz.json', 'w') as f:
    f.write(str_out)