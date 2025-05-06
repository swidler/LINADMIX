#!/usr/bin/python3

import numpy as np
import itertools
import pandas as pd

def read_ind_file(ind_file):
    samp_ids = []
    samp_sex = {}
    groups = []
    grouped_samples = {}
    samp_group = {}
    with open(ind_file, "rt") as indfile:
        for line in indfile:
            line = line.rstrip("\n")  # remove trailing line feeds
            line = line.lstrip()  # remove leading whitespace
            if not line:
                continue
            fields = line.split()  # split by whitespace
            samp_ids.append(fields[0])
            samp_sex[fields[0]] = fields[1]
            samp_group[fields[0]] = fields[2]
            groups.append(fields[2])
    no_samples = len(samp_ids)
    types = dict(zip(samp_ids, groups))
    #get number of groups
    no_groups = len(set(x for y in types for x in types.values()))  
    group_sizes = [(x, len(list(y))) for x,y in itertools.groupby(sorted(types.values()))]
    group_names = [x[0] for x in group_sizes]
    for grp in group_names:
        grouped_samples[grp] = [x for x in samp_ids if samp_group[x] == grp]  # dict of samples per group
    group_nums = [group_names.index(x)+1 for x in groups]
    groups_by_name = {group_names[x-1]:x for x in group_nums}
    groups_by_num = {x:group_names[x-1] for x in group_nums}
    print(f"individuals: {no_samples}")
    print(f"populations: {no_groups}")
    return {"samp_ids":samp_ids, "samp_sex":samp_sex, "samp_group":samp_group, "no_groups":no_groups, "grouped_samples": grouped_samples, "groups_sizes":group_sizes, "group_names":group_names, "group_nums":group_nums, "groups_by_name":groups_by_name, "groups_by_num":groups_by_num} 

def cross_ind_file(ind, fam_file):
    samples = []
    with open(fam_file, "rt") as famfile:
        for line in famfile:
            line = line.rstrip("\n")  # remove trailing line feeds
            line = line.lstrip()  # remove leading whitespace
            fields = line.split()  # split by whitespace
            samples.append(fields[1])
            
    #idx2 = [samples.index(x) for x in ind["samp_ids"]]  # this is the index of each ind file sample in the fam file
    idx = [ind["samp_ids"].index(x) if x in ind["samp_ids"] else np.nan for x in samples]  # this is the index of each fam file sample in the ind file
    print(f"Number of samples in dataset: {len(samples)}")
    print(f"Number of samples in ADMIXTURE: {len(ind['samp_ids'])}")
    if np.sum(np.isnan(idx)) > 0:
        print("ADMIXTURE uses samples that are not in the full list")
        
    #reorder ind file samples based on fam file order
    #ind["samp_ids"] = [samples[x] for x in idx if ~np.isnan(x) and x < len(samples)]  # exclude nans and out of range indices
    ind["samp_ids"] = samples[:]
    return ind

def read_pfile(p_file):
        samp_vals = []
        with open(p_file, "rt") as pfile:
            for line in pfile:
                line = line.rstrip("\n")  # remove trailing line feeds
                line = line.lstrip()  # remove leading whitespace
                fields = line.split()  # split by whitespace
                fields = [float(x) for x in fields]
                samp_vals.append(fields)
        return samp_vals
    
def read_gfile(g_file):
    gfile = pd.read_table(g_file,sep=' ')
    ids = gfile["IID"]  # get all sample ids
    gdata = gfile.iloc[:,6:]  # get all data excluding first 6 cols
    g = gdata.to_numpy()  # convert both data structures to numpy arrays
    g = np.where(g==1,0,g)  # replace 1 values with 0
    i = ids.to_numpy()
    return (i, g)

def read_gfile2(g_file):
    with open(g_file, "rt") as gfile:
        header = gfile.readline()
        header = header.rstrip("\n")  # remove trailing line feeds
        header = header.lstrip()  # remove leading whitespace
        fields = header.split()  # split by whitespace
        id_ind = fields.index("IID")
        gdata = []
        ids = []
        for line in gfile:
            line = line.rstrip("\n")  # remove trailing line feeds
            line = line.lstrip()  # remove leading whitespace
            fields = line.split()  # split by whitespace
            ids.append(fields[id_ind])
            f = np.array(fields[6:])
            f = np.array(['nan' if x == 'NA' else x for x in fields[6:]], dtype='float64')
            gdata.append(f)
    i = np.array(ids)
    g = np.asarray(gdata)
    #g = np.where(g==1,0,g)  # replace 1 values with 0
    return (i, g)
    
    
    
    