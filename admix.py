#!/usr/bin/python3

import copy
import numpy as np
from qpsolvers import solve_qp


class Admix:  
    def __init__(self, name="unnamed", description="", source="", K=0, samples={}, groups=[], Q={}, no_samples= 0):
        self.name = name
        self.description = description
        self.source = source
        self.K = copy.deepcopy(K)  # if these are just assigned (not copied), when creating 2 objects, these
        self.groups = copy.deepcopy(groups)  # elements will be identical and changing one will change the other
        self.Q = copy.deepcopy(Q)
        self.samples = copy.deepcopy(samples)
        self.no_samples = no_samples
        #self.no_K = len(K)
        

    def __repr__(self): #defines print of object
        return "name: %s\ndescription: %s\nsource: %s\nK: %s\nsamples: %s\ngroups: %s\nQ: %s\nno_samples: %s" % (self.name, self.description, self.source, self.K, self.samples, self.groups, self.Q, self.no_samples)

    def read_qfile(self, q_file):
        samp_vals = []
        with open(q_file, "rt") as qfile:
            for line in qfile:
                line = line.rstrip("\n")  # remove trailing line feeds
                line = line.lstrip()  # remove leading whitespace
                fields = line.split()  # split by whitespace
                fields = [float(x) for x in fields]
                samp_vals.append(fields)
        samps = len(samp_vals)
        if samps != self.no_samples:
            print(f"Number of individuals in ind file ({self.no_samples}) doesn't match number in q file ({samps})")
        K = len(samp_vals[0])
        Q = {self.samples["samp_ids"][x]:samp_vals[x] for x in range(len(samp_vals))}
        self.Q = Q
        self.K = K
        
    def get_error(self, no_samples, idx, V, no_joint_sourcepop, idx_sourcepop, sig):
        # draw error of each individual
        sample_err = np.zeros((no_samples,self.K))
        for samp in range(no_samples):
            zerovec = np.zeros(len(idx[samp]))
            sample_err[samp,idx[samp]]= np.random.multivariate_normal(zerovec,V[samp])
        
        # draw error due to population covariance
        pop_err = np.zeros(self.K)
        zerovec = np.zeros(no_joint_sourcepop)
        pop_err[idx_sourcepop] = np.random.multivariate_normal(zerovec,sig)
        return(sample_err, pop_err)
    
    def calc_randQ(self, pop, P, G, num_reps, ids, proxy_pop=None, minval=0.01):
        if proxy_pop == None:
            samples = self.samples["grouped_samples"][pop]
        else:
            samples = self.samples["grouped_samples"][proxy_pop]
        no_samples = len(samples)
        indices = [np.where(ids==x)[0][0] for x in samples]
        # compute variance
        [V, idx, is_K] = self.var(samples, P, G, indices, min_Qval=minval)
        # merge relevant source populations
        idx_sourcepop = set(np.concatenate(idx).tolist())
        idx_sourcepop = list(idx_sourcepop)
        no_joint_sourcepop = len(idx_sourcepop)
        # compute covariance
        sig = []
        for samp in samples:
            sig.append(np.array(self.Q[samp])[idx_sourcepop])
        sig = np.cov(sig, rowvar=False)
        meanQ = np.zeros((num_reps,self.K))
        for rep in range(num_reps):
            (sample_err, pop_err) = self.get_error(no_samples, idx, V, no_joint_sourcepop, idx_sourcepop, sig)
            
            # generate random q-vector per individual
            expq = np.zeros((no_samples,self.K))
            j = 0
            for samp in samples:
                expq[j,idx_sourcepop] = np.array(self.Q[samp])[idx_sourcepop]
                j += 1
            #randQ = bsxfun(@plus,expq+sample_err,pop_err);
            randQ = expq + sample_err + pop_err  # ??
    
            # generate mean q-vector
            meanQ[rep,:] = randQ.mean(axis=0)  # right?
            meanQ[rep,:] = np.maximum(meanQ[rep,:],0)
    
            # correct for last source population and normalize
            if any(is_K): # ??
                meanQ[rep,self.K-1] = 1 - np.sum(meanQ[rep,idx_sourcepop])
                meanQ[meanQ<0] = 0
            meanQ[meanQ<0] = 0  # this line seems redundant
            meanQ[rep,:] = meanQ[rep,:] / np.sum(meanQ[rep,:])
        return meanQ
    
    def calc_randQ_target(self, pop, P, G, num_reps, ids, is_Ks, randQ_sources, mix_coeff_vec, minval, proxy_pop=None, target_indices=None):
        if proxy_pop == None:
            samples = self.samples["grouped_samples"][pop]
        else:
            samples = self.samples["grouped_samples"][proxy_pop]
        no_samples = len(samples)
        indices = [np.where(ids==x)[0][0] for x in samples]
        
        # compute variance
        [V, idx, is_K] = self.var(samples, P, G, indices, min_Qval=minval, target_indices=target_indices)
        
        meanQ = np.zeros((num_reps,self.K))
        for rep in range(num_reps):
            source_num = len(randQ_sources.keys())
            source_qmat = {}
            for source in randQ_sources.keys():
                source_qmat[source] = randQ_sources[source][rep]
            h0qvec = np.dot(mix_coeff_vec, list(source_qmat.values()))
            h0qmat = np.zeros((no_samples,self.K))
            for samp in range(no_samples):
                h0qmat[samp] = h0qvec
            
            # compute covariance
            idx_satall = list(set([x for y in idx for x in y]))
            no_joint_sourcepop = len(idx_satall)
            sig = []
            for samp in samples:
                sig.append(np.array(self.Q[samp])[idx_satall])
            sig = np.cov(sig, rowvar=False)
            
            (sample_err, pop_err) = self.get_error(no_samples, idx, V, no_joint_sourcepop, idx_satall, sig)
            
            # calculate expq
            expq = h0qmat.copy()
            expq[:,-1] = 0

            # generate random q-vector per individual
            randQ = expq + sample_err + pop_err  # ??
    
            # generate mean q-vector
            meanQ[rep,:] = randQ.mean(axis=0)  # right?
            meanQ[rep,:] = np.maximum(meanQ[rep,:],0)
    
            # correct for last source population and normalize
            if is_Ks or any(is_K):
                meanQ[rep,self.K-1] = 1 - np.sum(meanQ[rep,:]);
                meanQ[meanQ<0] = 0
            meanQ[meanQ<0] = 0  # this line seems redundant
            meanQ[rep,:] = meanQ[rep,:] / np.sum(meanQ[rep,:])
        return meanQ

    def var(self, samps, P, G, indices, min_Qval=0.01, target_indices=[]):
        # hard-coded arguments
        no_snps = len(P)
        no_inds = len(samps)
        V = []
        idx_sourcepop = []
        is_K = np.zeros(no_inds)
        i = 0
        # process individuals
        for samp in samps:
            # indices of source populations contributing to individual
            idx_sourcepop.append(np.where(np.array(self.Q[samp])>=min_Qval)[0])
            # special treatment of last (Kth) population
            if idx_sourcepop[i][-1] == self.K-1:  # what is this?!
                is_K[i] = 1
                idx_sourcepop[i] = idx_sourcepop[i][:-1]  # delete last element
            if target_indices != []:
                idx_sourcepop[i] = np.union1d(idx_sourcepop[i], target_indices) ##TEST
            # num of relevant source populations
            L = len(idx_sourcepop[i]);
            # initialize information matrix
            I = np.zeros((L,L))
            # compute the information matrix
            Ij = np.zeros((no_snps,L))
            counter = 0
            for pop in idx_sourcepop[i]:
                for snp in range(no_snps):
                     if ~np.isnan(G[indices[i]][snp]):
                         denom1 = (np.dot(np.array(self.Q[samp])[idx_sourcepop[i]],np.array(P[snp])[idx_sourcepop[i]]))**2
                         denom2 = np.dot(np.array(self.Q[samp])[idx_sourcepop[i]],(1 - np.array(P[snp])[idx_sourcepop[i]]))**2
                         Ij[snp] = (P[snp][pop] - P[snp][self.K-1]) * (np.array(P[snp])[idx_sourcepop[i]] - P[snp][self.K-1]) * ( (2 - G[indices[i]][snp]) / denom1 + G[indices[i]][snp] / denom2 )
                I[counter,:] = Ij.sum(axis=0)
                counter += 1
            V.append(np.linalg.inv(I))
            i += 1
        return(V, idx_sourcepop, is_K)
        
    def linadmix(self, target, sources, q_val=0, fid=None):
        samples = self.samples["grouped_samples"][target]
        qvals = []
        for samp in samples:
            qvals.append(np.array(self.Q[samp]))
        q_target = np.mean(qvals, axis=0)
        if fid:
            print(f"Target:\n\t{target} ({len(samples)} individuals)")
            fid.write(f"Target:\n\t{target} ({len(samples)} individuals)\n")
        q_sources = np.zeros((self.K, len(sources)))
        i = 0
        if fid:
            print("Sources: ")
            fid.write("Sources: \n")
        for source in sources:
            q_per_source = []
            samples = self.samples["grouped_samples"][source]
            for samp in samples:
                q_per_source.append(np.array(self.Q[samp]))
            for j in range(self.K):
                q_sources[j,i] = np.mean(np.array(q_per_source)[:,j])
            if fid:
                print(f"({i+1}) {source} ({len(samples)} individuals)")
                fid.write(f"({i+1}) {source} ({len(samples)} individuals)\n")
            i+=1
        P = 2*(np.dot(q_sources.T,q_sources))
        q = -2*(np.dot(q_sources.T,q_target))
        A = -1*np.eye(len(q_sources[0]))
        b = np.zeros(len(q_sources[0]))
        Aeq = np.ones((1,len(q_sources[0])))
        beq = np.array([1])
        mix_coeff = solve_qp(P, q, A, b, Aeq, beq)
        #resid = np.dot(q_sources, max_coeff) - q_target  # these aren't needed anymore
        #resnorm = np.dot(resid.T, resid)
        if q_val != 0:
            q_exp = np.dot(q_sources,mix_coeff)
            idx = np.where(q_exp>=q_val)
            wval = np.sum((q_target[idx]-q_exp[idx])**2/q_exp[idx])
            return (mix_coeff, wval)
        return mix_coeff
    
    def linadmixerr(self, target, sources, randQ, pops):
        num_reps = len(randQ[pops[0]])
        num_sources = len(sources)
        mixing = np.zeros((num_sources, num_reps))
        norms = np.zeros(num_reps)
        randQ_sources = {}
        # identify randQ of target
        if target in pops:
            randQ_target = randQ[target]
        else:
            raise Exception(f"Target {target} is missing in population list")
        
        # identify randQ of sources
        for source in sources:
            if source in pops:
                randQ_sources[source] = randQ[source]
            else:
               raise Exception(f"Source {source} is missing in population list") 
        
        # prepare linear regression
        A = -1*np.eye(num_sources)
        b = np.zeros(num_sources)
        Aeq = np.ones((1,num_sources))
        beq = np.array([1])
        
        # compute for all repetitions
        q_sources = np.zeros((len(randQ[pops[0]][0]), num_sources))
        for rep in range(num_reps):
            # build Q-matrix of sources
            for source in range(num_sources):
                q_sources[:,source] = randQ_sources[sources[source]][rep]

            # build q-matrix of target
            q_target = randQ_target[rep]
            
            P = 2*(np.dot(q_sources.T,q_sources))
            q = -2*(np.dot(q_sources.T,q_target))
            
            # run linear model
            mixing[:,rep] = solve_qp(P, q, A, b, Aeq, beq)
        
            
        
        mix_coeff = {}
        mix_coeff["mean"] = np.mean(mixing, axis=1)
        mix_coeff["std"] = np.std(mixing, axis=1, ddof=1)
        
        return mix_coeff
    
    def get_qvecs(self, sources):
        qvecs = {}
        #print("Populations:")
        for source in sources:
            q = []
            num_inds = len(self.samples["grouped_samples"][source])
            for ind in range(num_inds):
                q.append(self.Q[self.samples["grouped_samples"][source][ind]])
            qvecs[source] = q
            #print(f"{source} ({num_inds} individuals)")
        return qvecs
    
    def comp_target_sources(self, qvecs, sources, mix, qval):
        qmat_sources = np.zeros((len(sources),self.K))
        is_K = 0
        for source in range(len(sources)):
            qmat_sources[source,:] = np.mean(qvecs[sources[source]], axis=0)
        #vals = list(mix.values())
        qest_target = np.dot(mix,qmat_sources)
        idx_est_target = np.where(qest_target>=qval)[0]
        if idx_est_target[-1] == self.K + 1:
            is_K = 1
            idx_est_target = idx_est_target[:-1]  # delete last element
        return (idx_est_target, is_K)
        
    
        
        
                
            
  
if __name__ == "__main__":
    ad = Admix()
    print(ad)