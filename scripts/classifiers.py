"""
script to define classifier models for WINDEX emissions 
references: hannah snell, dr. lauren sugden, scott mccallum
"""

import numpy as np
import math
import os
import pickle

SQRT_2PI = math.sqrt(2 * math.pi)

class Stats():
    def __init__(self,path):
        if path != '' and path[-1] != '/':
            path += '/'
        file = open(path+'component_stats.txt')
        f = file.read()
        file.close()
        f = f.strip().splitlines()
        self.stats = [line.split('\t')[0] for line in f]
        self.stat2score = {s:np.nan for s in self.stats}
    def set_stat(self,stat,value):
        self.stat2score[stat] = value
        return self.stat2score

class Classifiers():
    def __init__(self,path2trained, odecompensation = True):

        self.classes = ['neutral','sweep', 'linked']
        file = open(path2trained+'component_stats.txt','r')
        f = file.read()
        file.close()
        self.statlist = f.strip().splitlines()
        self.num2stat = {i:self.statlist[i] for i in range(len(self.statlist))}
        self.stat2num = {y:x for x,y in self.num2stat.items()}
        self.path2AODE = os.path.join(path2trained,'AODE_params/')
        self.MARGINALS = [[] for _ in self.stat2num.keys()]
        for stat in self.statlist:
            statnum = self.stat2num[stat]
            for C in self.classes:
                self.MARGINALS[statnum].append(self._extract_marginal_params(self.gmm_fit_1D(stat,C)))
        self.odecompensation = odecompensation

        self.JOINTS = [[[] for _ in self.stat2num.keys()] for _ in self.stat2num.keys()]

        for stat1 in self.statlist:
            for stat2 in [x for x in self.statlist if x!= stat1]:
                statnum1 = self.stat2num[stat1]
                statnum2 = self.stat2num[stat2]
                if statnum1 < statnum2:
                    for C in self.classes:
                        self.JOINTS[statnum1][statnum2].append(self._extract_joint_params(self.gmm_fit(stat1,stat2,C)))

    def classify(self,pivec):
        numerators = []
        for i in range(len(self.classes)):
            C = self.classes[i]
            pi = pivec[i]
            liklhd = self.nb_likelihood(C)
            numerators.append(pi*liklhd)
        denom = sum(numerators)
        posteriors = [float(x)/denom for x in numerators]
        return posteriors
    
    def nb_likelihood(self,C,S):
        product = 1
        vallist = []
        for stat in self.statlist:
            score = S.stat2score[stat]
            if score != -998:
                MARGS = self.MARGINALS[self.stat2num[stat]] # marginals for stat
                M = MARGS[self.classes.index(C)] # for specific class
                value = self.GMM_pdf(M,score)
                product = product*value
                vallist.append(value)

        if len(vallist) == 0: 
            #print("Error: one or more input rows have no defined statistics. Please remove from your data or recalculate.")
            return None
            
        if self.odecompensation:
            n = len(self.statlist) - len(vallist)
            xbar = float(sum(vallist))/len(vallist)
            for _ in range(n):
                product = product*xbar
        return product

    def _extract_marginal_params(self, G): # one-time conversion of fitted 1D GMM into plain floats (call once)
        params = []
        for i in range(len(G.weights_)): 
            w = float(G.weights_[i])
            mu = float(G.means_[i].item())
            sigma = math.sqrt(float(G.covariances_[i].item()))
            params.append((w, mu, sigma))
        return params

    def _extract_joint_params(self, G): # one-time conversion of fitted 1D GMM into plain floats (call once)
        params = []
        for i in range(len(G.weights_)): 
            w = float(G.weights_[i])
            mu1 = float(G.means_[i][0])
            mu2 = float(G.means_[i][1])
            sigma1 = math.sqrt(float(G.covariances_[i][0][0]))
            sigma2 = math.sqrt(float(G.covariances_[i][1][1]))
            ro = float(G.covariances_[i][0][1]) / (sigma1 * sigma2)
            cond_sigma_1given2 = sigma1 * math.sqrt(1 - ro**2)
            cond_sigma_2given1 = sigma2 * math.sqrt(1 - ro**2)
            params.append((w, mu1, mu2, sigma1, sigma2, ro, cond_sigma_1given2, cond_sigma_2given1))
        return params

    def aode_likelihood(self, C, S): # from ode function in SWIFr.py; return aode likelihood for a single classification scenario
        ode_likelihoods = []
        scenarionum = self.classes.index(C)
        for keystat in self.statlist:
            score = S.stat2score[keystat]
            if score != -998: # changed to check if stat value is -998 
                Likelihood_list = []
                MARGS = self.MARGINALS[self.stat2num[keystat]]
                M = MARGS[self.classes.index(C)]
                value = self.GMM_pdf(M, score)
                Likelihood = value # marginal likelihood for keystat
                stats_undefined = 0 # keep track of how many comparison stats are undefined
                for stat in self.statlist:
                    if stat != keystat:
                        score2 = S.stat2score[stat]
                        if score2 == -998:
                            stats_undefined += 1
                        else:
                            if self.stat2num[keystat] < self.stat2num[stat]:
                                H = self.conditional_GMM(score, 2, self.JOINTS[self.stat2num[keystat]][self.stat2num[stat]][scenarionum])
                                value = self.GMM_pdf(H, score2)
                                Likelihood = Likelihood * value
                                Likelihood_list.append(value)
                            else:
                                H = self.conditional_GMM(score, 1, self.JOINTS[self.stat2num[stat]][self.stat2num[keystat]][scenarionum])
                                value = self.GMM_pdf(H, score2)
                                Likelihood = Likelihood * value
                                Likelihood_list.append(value)
                
                if self.odecompensation: # do this 
                    n = stats_undefined
                    if len(Likelihood_list) > 0: # else: no other stats were defined so Likelihood is single-stat marginal
                        xbar = float(sum(Likelihood_list))/len(Likelihood_list) # average conditional density for defined comparison stats
                        for i in range(n):
                            Likelihood = Likelihood * xbar

                ode_likelihoods.append(Likelihood)

        if len(ode_likelihoods) == 0: # check to make sure ode_likelihoods is non-empty
            return None
        
        # return average of ode likelihoods for aode emission of specific class
        return sum(ode_likelihoods) / len(ode_likelihoods)

    def gmm_fit_1D(self,stat,scenario):
        G = pickle.load(open(os.path.join(self.path2AODE,stat+'_'+scenario+'_1D_GMMparams.p'),'rb'))
        return G

    def gmm_fit(self,stat1,stat2,scenario):
        G = pickle.load(open(self.path2AODE+stat1+'_'+stat2+'_'+scenario+'_GMMparams.p','rb'))
        return G

    def GMM_pdf(self, params, x):
        pdf = 0.0
        for w, mu, sigma in params: 
            u = (x - mu) / sigma 
            pdf += w * math.exp(-u * u/2) / (SQRT_2PI * sigma)
        return pdf

    def normpdf(self,x,mu,sigma):
        u = float(x-mu)/sigma
        y = (1/(math.sqrt(2*math.pi)*abs(sigma)))*math.exp(-u*u/2)
        return y
	    
    def conditional_GMM(self, condval, keystat, joint_params):
        H = [] 
        if  keystat == 1: 
            for w, mu1, mu2, sigma1, sigma2, ro, cs1g2, _ in joint_params: 
                cond_mu = mu1 + sigma1 * ro * (condval - mu2) / sigma2
                H.append((w, cond_mu, cs1g2))
        else: 
            for w, mu1, mu2, sigma1, sigma2, ro, _, cs2g1 in joint_params: 
                cond_mu = mu2 + sigma2 * ro * (condval - mu1) / sigma1
                H.append((w, cond_mu, cs2g1))
        return H
        
class Mix1D():
    def __init__(self):
        self.means_ = []
        self.weights_ = []
        self.covariances_ = []
