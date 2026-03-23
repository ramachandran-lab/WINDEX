"""
script to define classifier models for WINDEX emissions 
references: hannah snell, dr. lauren sugden, scott mccallum
"""
import numpy as np
import math
import os
import pickle

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
    def __init__(self,path2trained, odecompensation=True):

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
                self.MARGINALS[statnum].append(self.gmm_fit_1D(stat,C))
        self.odecompensation = odecompensation

        self.JOINTS = [[[] for _ in self.stat2num.keys()] for _ in self.stat2num.keys()]

        for stat1 in self.statlist:
            for stat2 in [x for x in self.statlist if x!= stat1]:
                statnum1 = self.stat2num[stat1]
                statnum2 = self.stat2num[stat2]
                if statnum1 < statnum2:
                    for C in self.classes:
                        self.JOINTS[statnum1][statnum2].append(self.gmm_fit(stat1,stat2,C))

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

    def ode_likelihood(self,C, S): # adapted from ode function in SWIFr.py to only return likelihood for a single classification scenario
        ode_likelihoods = []
        for keystat in self.statlist:
            score = S.stat2score[keystat]
            if np.isnan(score):
                return 'n/a'
            else:
                Likelihood_list = []
                MARGS = self.MARGINALS[self.stat2num[keystat]]
                scenarionum = self.classes.index(C)
                M = MARGS[self.classes.index(C)]
                value = self.GMM_pdf(M,score)
                Likelihood = value
                Likelihood_list.append(value)
                stats_undefined = 0 # keep track of how many comparison stats are undefined
                for stat in self.statlist:
                    if stat != keystat:
                        score2 = S.stat2score[stat]
                        if np.isnan(score2):
                            stats_undefined += 1
                        else:
                            if self.stat2num[keystat] < self.stat2num[stat]:
                                H = self.conditional_GMM(score,2,self.JOINTS[self.stat2num[keystat]][self.stat2num[stat]][scenarionum])
                                value = self.GMM_pdf(H,score2)
                                Likelihood = Likelihood*value
                                Likelihood_list.append(value)
                            else:
                                H = self.conditional_GMM(score,1,self.JOINTS[self.stat2num[stat]][self.stat2num[keystat]][scenarionum])
                                value = self.GMM_pdf(H,score2)
                                Likelihood = Likelihood*value
                                Likelihood_list.append(value)
                if self.odecompensation:
                    n = stats_undefined
                    xbar = float(sum(Likelihood_list))/len(Likelihood_list) # average conditional density for defined comparison stats
                    for i in range(n):
                        Likelihood = Likelihood*xbar
    
            ode_likelihoods.append(Likelihood)
        
        # return average of ode likelihoods for aode emission
        return sum(ode_likelihoods) / len(ode_likelihoods)

    def gmm_fit_1D(self,stat,scenario):
        G = pickle.load(open(os.path.join(self.path2AODE,stat+'_'+scenario+'_1D_GMMparams.p'),'rb'))
        return G

    def gmm_fit(self,stat1,stat2,scenario):
        G = pickle.load(open(self.path2AODE+stat1+'_'+stat2+'_'+scenario+'_GMMparams.p','rb'))
        return G

    def GMM_pdf(self,G,x):
        w = G.weights_
        mu = G.means_
        C = G.covariances_
        pdf = 0
        for i in range(len(w)):
            pdf += w[i]*self.normpdf(x,mu[i][0],math.sqrt(C[i][0]))
        return pdf

    def normpdf(self,x,mu,sigma):
        u = float(x-mu)/sigma
        y = (1/(math.sqrt(2*math.pi)*abs(sigma)))*math.exp(-u*u/2)
        return y
	    
    def conditional_GMM(self,condval,keystat,G):
        # keystat = 1 if want stat1|stat2, keystat = 2 if want stat2|stat1
        H = Mix1D()

        for i in range(len(G.weights_)):
            sigma1 = math.sqrt(G.covariances_[i][0][0])
            sigma2 = math.sqrt(G.covariances_[i][1][1])
            ro = float(G.covariances_[i][0][1])/(sigma1*sigma2)
            mu1 = G.means_[i][0]
            mu2 = G.means_[i][1]
            if keystat == 1:
                H.weights_.append(G.weights_[i])
                H.means_.append([mu1 + float(sigma1*ro*(condval-mu2))/sigma2])
                H.covariances_.append([(1-ro**2)*sigma1**2])

            elif keystat == 2:
                H.weights_.append(G.weights_[i])
                H.means_.append([mu2 + float(sigma2*ro*(condval-mu1))/sigma1])
                H.covariances_.append([(1-ro**2)*sigma2**2])
        return H

class Mix1D():
    def __init__(self):
        self.means_ = []
        self.weights_ = []
        self.covariances_ = []
