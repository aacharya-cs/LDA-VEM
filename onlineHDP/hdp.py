# the hdp class in python
# implements the truncated hdp model
# by chongw@cs.princeton.edu

import numpy as np
import scipy.special as sp
import scipy
import scipy.linalg
import os, sys, math, time
import utils
from utils import *
from corpus import document, corpus, parse_line
from itertools import izip
import random
sys.path.insert(0, "/lusr/u/ayan/Documents/DSLDA_SDM/DSLDA/liblinear-1.92/python/")
from liblinearutil import *

meanchangethresh = 0.00001
random_seed = 999931111
np.random.seed(random_seed)
random.seed(random_seed)

def save_beta_ss(filename, m_var_beta_ss):
        f = file(filename, "w") 
        for beta in m_var_beta_ss:
            line = ' '.join([str(x) for x in beta])  
            f.write(line + '\n')
        f.close()

class hdp_hyperparameter:
    def __init__(self, alpha_a, alpha_b, gamma_a, gamma_b, hyper_opt=False):
        self.m_alpha_a = alpha_a
        self.m_alpha_b = alpha_b
        self.m_gamma_a = gamma_a
        self.m_gamma_b = gamma_b
        self.m_hyper_opt = hyper_opt 

class suff_stats:
    def __init__(self, T, size_vocab, K, nsz):
        self.m_var_sticks_ss = np.zeros(T) 
        self.m_var_beta_ss = np.zeros((T, size_vocab))
        self.m_var_zeta = np.zeros((nsz, K))
    
    def set_zero(self):
        self.m_var_sticks_ss.fill(0.0)
        self.m_var_beta_ss.fill(0.0)
        self.m_var_zeta.fill(0.0)


class hdp:
    ''' hdp model using john's new stick breaking'''
    def save_topics(self, filename):
        f = file(filename, "w") 
        for beta in self.m_beta:
            line = ' '.join([str(x) for x in beta])  
            f.write(line + '\n')
        f.close()
    

    def write_local_sticks(self, v):
        filename = "localsticks_a.txt"; 
        f = file(filename, "a") 
        for x in v[0]:  
            f.write(str(x) + '\t')
        f.write('\n')
	f.close()
        filename = "localsticks_b.txt"; 
        f = file(filename, "a")
        for y in v[1]:  
            f.write(str(y) + '\t')
        f.write('\n')
	f.close()

    def __init__(self, T, K, D,  size_vocab, eta, trsz, hdp_hyperparam):
        ''' this follows the convention of the HDP paper'''
        ''' gamma, first level concentration ''' 
        ''' alpha, second level concentration '''
        ''' eta, the topic Dirichlet '''
        ''' T, top level truncation level '''
        ''' K, second level truncation level '''
        ''' size_vocab, size of vocab'''
        ''' hdp_hyperparam, the hyperparameter of hdp '''
    
        self.m_hdp_hyperparam = hdp_hyperparam

        self.m_T = T # higher level truncation
        self.m_K = K # for now, we assume all the same for the second level truncation
        self.m_size_vocab = size_vocab

        # print "%d %d %d" %(T, size_vocab, D)

        self.m_beta = np.random.gamma(1.0, 1.0, (T, size_vocab)) * D*100/(T*size_vocab)
        (log_m_beta, log_norm) = utils.log_normalize(self.m_beta)
        self.m_beta = np.exp(log_m_beta)
	self.save_topics("lambda.txt");
        self.m_eta = eta  

        self.m_alpha = hdp_hyperparam.m_alpha_a/hdp_hyperparam.m_alpha_b
        self.m_gamma = hdp_hyperparam.m_gamma_a/hdp_hyperparam.m_gamma_b
        self.m_var_sticks = np.zeros((2, T-1))
        self.m_var_sticks[0] = 1.0
        self.m_var_sticks[1] = self.m_gamma
        self.r = np.zeros((6, self.m_K))
        self.dmu = np.zeros((trsz, 6)) 

        # variational posterior parameters for hdp
        self.m_var_gamma_a = hdp_hyperparam.m_gamma_a
        self.m_var_gamma_b = hdp_hyperparam.m_gamma_b

    # E-step per document
    def doc_e_step(self, doc, ss, trlabel, docnum, Elogbeta, Elogsticks_1st, Elogsticks_2nd, var_converge, fresh=False):

        Elogbeta_doc = Elogbeta[:, doc.words] 
        v = np.zeros((2, self.m_K-1))

        phi = np.ones((doc.length, self.m_K)) * 1.0/self.m_K  # should be zeta

        likelihood = 0.0
        old_likelihood = -1e1000
        converge = 1.0 
        eps = 1e-100
        
        iter = 0
        max_iter = 10
        #(TODO): support second level optimization in the future
        while iter < max_iter: #and (converge < 0.0 or converge > var_converge):
            ### update variational parameters
            # smallphi
            
            var_phi = np.dot(phi.T, (Elogbeta_doc * doc.counts).T) + Elogsticks_1st
            (log_var_phi, log_norm) = utils.log_normalize(var_phi)
            var_phi = np.exp(log_var_phi)

 
            # phi  #zeta
            sval   = np.zeros((1, self.m_K))
            nwords = np.sum(doc.counts)
            tmp    = (self.r[trlabel,:] - self.r)
            sval   = np.dot(self.dmu[docnum,:],tmp) 
            sval   = sval/nwords
	    sval   = 0;

            phi = np.dot(var_phi, Elogbeta_doc).T + Elogsticks_2nd + sval
            (log_phi, log_norm) = utils.log_normalize(phi)
            phi = np.exp(log_phi)
            phi_all = phi * np.array(doc.counts)[:,np.newaxis]


            # local sticks
            v[0] = 1.0 + np.sum(phi_all[:,:self.m_K-1], 0)      #a_{nt}
            phi_cum = np.flipud(np.sum(phi_all[:,1:], 0))
            v[1] = self.m_alpha + np.flipud(np.cumsum(phi_cum)) #b_{nt}
            Elogsticks_2nd = expect_log_sticks(v)
            
            if iter==max_iter-1: 
              self.write_local_sticks(v);

            likelihood = 0.0
            # compute likelihood
            # var_phi part/ C in john's notation
            likelihood += np.sum((Elogsticks_1st - log_var_phi) * var_phi)

            # v part/ v in john's notation, john's beta is alpha here
            log_alpha = np.log(self.m_alpha)
            likelihood += (self.m_K-1) * log_alpha
            dig_sum = sp.psi(np.sum(v, 0))
            likelihood += np.sum((np.array([1.0, self.m_alpha])[:,np.newaxis]-v) * (sp.psi(v)-dig_sum))
            likelihood -= np.sum(sp.gammaln(np.sum(v, 0))) - np.sum(sp.gammaln(v))

            # Z part 
            likelihood += np.sum((Elogsticks_2nd - log_phi) * phi)

            # X part, the data part
            likelihood += np.sum(phi.T * np.dot(var_phi, Elogbeta_doc * doc.counts))

            converge = (likelihood - old_likelihood)/abs(old_likelihood)
            old_likelihood = likelihood
            
            if converge < 0:
                print "warning, likelihood is decreasing!"
            
            iter = iter + 1
            
        # update the suff_stat ss 
        ss.m_var_sticks_ss += np.sum(var_phi, 0)   
        ss.m_var_beta_ss[:, doc.words] += np.dot(var_phi.T, phi.T * doc.counts)
        ss.m_var_zeta[docnum,:] = np.sum((phi.T * doc.counts).T,0)

        return(likelihood)

    def optimal_ordering(self, ss):
        s = [(a, b) for (a,b) in izip(ss.m_var_sticks_ss, range(self.m_T))]
        x = sorted(s, key=lambda y: y[0], reverse=True)
        idx = [y[1] for y in x]
        print idx
        filename = "idx.txt"; 
        f = file(filename, "w")
        for y in idx:  
            f.write(str(y) + '\t')
        f.write('\n')
	f.close()

        ss.m_var_sticks_ss[:] = ss.m_var_sticks_ss[idx]
        ss.m_var_beta_ss[:] = ss.m_var_beta_ss[idx,:]

    def do_m_step(self, ss, trlabels):

        print "M-step1"       
        self.optimal_ordering(ss)

        ## update global sticks 
        self.m_var_sticks[0] = ss.m_var_sticks_ss[:self.m_T-1] + 1.0  #u_{k}
        var_phi_sum = np.flipud(ss.m_var_sticks_ss[1:])
        self.m_var_sticks[1] = np.flipud(np.cumsum(var_phi_sum)) + self.m_gamma #v_{k}

        ## update topic parameters and normaize them
        self.m_beta = self.m_eta + ss.m_var_beta_ss #lambda_{kv}
        beta_sum = np.sum(self.m_beta, axis=1)
        self.m_beta = self.m_beta/beta_sum[:, np.newaxis]

        self.save_topics("lambda_updated.txt");
        save_beta_ss("beta_ss.txt", ss.m_var_beta_ss);

        raw_input()

        if self.m_hdp_hyperparam.m_hyper_opt:     
            print "M-step2"
            m = 2;  
            ## running SVM
            """x = utils.converttoliblinear(ss.m_var_zeta)
	    y = trlabels    
            m = train(y, x, '-c 1 -s 4')
            save_model('zzzzzz',m)
            self.dmu = utils.read_2D_array('alphaval.txt')
            self.dmu = scipy.array(self.dmu)
            self.r = utils.read_2D_array('wval.txt')
            self.r = (scipy.array(self.r)).T
            p_label, p_acc, p_val = predict(y,x,m)
            ACC, MSE, SCC = evaluations(y, p_label)
            print "accuracy: %f" %ACC"""

            """self.m_var_gamma_a = self.m_hdp_hyperparam.m_gamma_a + self.m_T - 1
            dig_sum = sp.psi(np.sum(self.m_var_sticks, 0))
            Elog1_W = sp.psi(self.m_var_sticks[1]) - dig_sum
            self.m_var_gamma_b = self.m_hdp_hyperparam.m_gamma_b - np.sum(Elog1_W)
            self.m_gamma = self.m_hdp_hyperparam.m_gamma_a/self.m_hdp_hyperparam.m_gamma_b"""

        return m
    
    ## one iteration of the em
    def em(self, c, trlabels, var_converge, nsz, fresh):

        # initializing sufficient statistics
        ss = suff_stats(self.m_T, self.m_size_vocab, self.m_K, nsz)
        ss.set_zero()
        
        # prepare all needs for a single doc
        Elogbeta = dirichlet_expectation(self.m_beta) # the topics
        Elogsticks_1st = expect_log_sticks(self.m_var_sticks) # global sticks 
        Elogsticks_2nd = expect_log_sticks(self.m_var_sticks[:,:self.m_K-1]) # local sticks; have to change later 
        likelihood = 0.0
        docnum = 0
        for doc in c.docs:
            # print docnum
            likelihood += self.doc_e_step(doc, ss, trlabels[docnum], docnum, Elogbeta, Elogsticks_1st, Elogsticks_2nd, var_converge, fresh=fresh)
            docnum += 1
 
        # collect the likelihood from other parts
        # the prior for gamma
        if self.m_hdp_hyperparam.m_hyper_opt:
            log_gamma = sp.psi(self.m_var_gamma_a) -  np.log(self.m_var_gamma_b)
            likelihood += self.m_hdp_hyperparam.m_gamma_a * np.log(self.m_hdp_hyperparam.m_gamma_b) \
                    - sp.gammaln(self.m_hdp_hyperparam.m_gamma_a) 

            likelihood -= self.m_var_gamma_a * np.log(self.m_var_gamma_b) \
                    - sp.gammaln(self.m_var_gamma_a) 

            likelihood += (self.m_hdp_hyperparam.m_gamma_a - self.m_var_gamma_a) * log_gamma \
                    - (self.m_hdp_hyperparam.m_gamma_b - self.m_var_gamma_b) * self.m_gamma
        else:
            log_gamma = np.log(self.m_gamma)

        # the W/sticks part 
        likelihood += (self.m_T-1) * log_gamma
        dig_sum = sp.psi(np.sum(self.m_var_sticks, 0))
        likelihood += np.sum((np.array([1.0, self.m_gamma])[:,np.newaxis] - self.m_var_sticks) * (sp.psi(self.m_var_sticks) - dig_sum))
        likelihood -= np.sum(sp.gammaln(np.sum(self.m_var_sticks, 0))) - np.sum(sp.gammaln(self.m_var_sticks))
        
        # the beta part    
        likelihood += np.sum((self.m_eta - self.m_beta) * Elogbeta)
        likelihood += np.sum(sp.gammaln(self.m_beta) - sp.gammaln(self.m_eta))
        likelihood += np.sum(sp.gammaln(self.m_eta*self.m_size_vocab) - sp.gammaln(np.sum(self.m_beta, 1)))

        # run m step
        model = self.do_m_step(ss, trlabels) 
        return (likelihood, model)

    # inference of local variational parameters only in test phase
    def doc_inference(self, doc, docnum, Elogbeta, Elogsticks_1st, var_converge, m_var_zeta):

        Elogbeta_doc = Elogbeta[:, doc.words] 
        v = np.zeros((2, self.m_K-1))         

        phi = np.ones((doc.length, self.m_K)) * 1.0/self.m_K  # should be zeta

        # the following line is of no use
        Elogsticks_2nd = expect_log_sticks(v)

        likelihood = 0.0
        old_likelihood = -1e1000
        converge = 1.0 
        eps = 1e-100
        
        iter = 0
        max_iter = 100
        #(TODO): support second level optimization in the future
        while iter < 20: #and (converge < 0.0 or converge > var_converge):
            ### update variational parameters
            # var_phi

            var_phi = np.dot(phi.T, (Elogbeta_doc * doc.counts).T) + Elogsticks_1st
            (log_var_phi, log_norm) = utils.log_normalize(var_phi)
            var_phi = np.exp(log_var_phi)

            # phi  #zeta
            phi = np.dot(var_phi, Elogbeta_doc).T + Elogsticks_2nd 
            (log_phi, log_norm) = utils.log_normalize(phi)
            phi = np.exp(log_phi)
            phi_all = phi * np.array(doc.counts)[:,np.newaxis]
              
            # local sticks
            v[0] = 1.0 + np.sum(phi_all[:,:self.m_K-1], 0)  #a_{jt}
            phi_cum = np.flipud(np.sum(phi_all[:,1:], 0))
            v[1] = self.m_alpha + np.flipud(np.cumsum(phi_cum)) #b_{jt}
            Elogsticks_2nd = expect_log_sticks(v)

            likelihood = 0.0
            # compute likelihood
            # var_phi part/ C in john's notation
            likelihood += np.sum((Elogsticks_1st - log_var_phi) * var_phi)

            # v part/ v in john's notation, john's beta is alpha here
            log_alpha = np.log(self.m_alpha)
            likelihood += (self.m_K-1) * log_alpha
            dig_sum = sp.psi(np.sum(v, 0))
            likelihood += np.sum((np.array([1.0, self.m_alpha])[:,np.newaxis]-v) * (sp.psi(v)-dig_sum))
            likelihood -= np.sum(sp.gammaln(np.sum(v, 0))) - np.sum(sp.gammaln(v))

            # Z part 
            likelihood += np.sum((Elogsticks_2nd - log_phi) * phi)

            # X part, the data part
            likelihood += np.sum(phi.T * np.dot(var_phi, Elogbeta_doc * doc.counts))

            converge = (likelihood - old_likelihood)/abs(old_likelihood)
            old_likelihood = likelihood
            
            """if converge < 0:
                print "warning, likelihood is decreasing!" """
            
            iter = iter + 1
        
        m_var_zeta[docnum,:] = np.sum((phi.T * doc.counts).T,0)

        return(likelihood, m_var_zeta)

