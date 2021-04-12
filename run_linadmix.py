#!/usr/bin/python3

import linadmix_tools as l
import admix as ad
import numpy as np
from config import *


# read in ind file 
ind = l.read_ind_file(indfile_main)  # ind is a dict
ind = l.cross_ind_file(ind, indfile_admix)

# create admix object
ad = ad.Admix(samples=ind, no_samples=len(ind["samp_ids"]))
ad.read_qfile(Qfile)

# read in p and g files
P = l.read_pfile(Pfile)
(ids, G) = l.read_gfile(Gfile)

list_of_populations = target_pops + source_pops_var + source_pops_const
randQ = {}
for pop in list_of_populations:
    randQ[pop] = ad.calc_randQ(pop, P, G, num_reps, ids)
fid = open(outfile, "w")
for source in source_pops_var:
    sources = source_pops_const + [source]
    for target in target_pops:
        mix_coeff = ad.linadmix(target, sources, fid=fid)
        mix_coeff_rand = ad.linadmixerr(target, sources, randQ, list_of_populations)
        std = mix_coeff_rand["std"]
        (mix_coeff_pval, pval) = ad.linadmix(target, sources, q_val=qval)
        qvecs = ad.get_qvecs(sources)
        (target_indices, is_K) = ad.comp_target_sources(qvecs, sources, mix_coeff_pval, qval)
        randQ_sources = {}
        for pop in sources:
            randQ_sources[f"{pop}"] = ad.calc_randQ(pop, P, G, num_reps_pval, ids, minval=qval)
        randQ_target = ad.calc_randQ_target(target, P, G, num_reps_pval, ids, is_K, randQ_sources, mix_coeff, minval=qval, target_indices=target_indices);
        wvalstat = np.zeros(num_reps_pval)
        randexpq = np.zeros((num_reps_pval, ad.K))
        
        # calculate distribution
        for rep in range(num_reps_pval):
            q_target = randQ_target[rep]
            num_sources = len(sources)
            q_sources = np.zeros((num_sources, ad.K))
            i = 0
            for pop in sources:
                q_sources[i] = (randQ_sources[f"{pop}"][rep])
                i += 1
            randexpq[rep] = np.dot(mix_coeff, q_sources)
            q_expected = np.dot(mix_coeff, q_sources)
            idx = np.where(q_expected>=qval)
            wvalstat[rep] = np.sum((q_target[idx]-q_expected[idx])**2/q_expected[idx])
            
        emp_pval = (len(np.where(wvalstat>=pval)[0])+1)/(num_reps_pval+1)
        print(f"Mixing coefficients: {mix_coeff}\nStandard errors: {std}\nP-value: {emp_pval}")    
        fid.write(f"Mixing coefficients: {mix_coeff}\nStandard errors: {std}\nP-value: {emp_pval}\n\n")
fid.close()