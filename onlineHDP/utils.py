#! /usr/bin/python

''' several useful functions '''
import numpy as np
import scipy.special as sp

def dirichlet_expectation(alpha):
    if (len(alpha.shape) == 1):
        return(sp.psi(alpha) - sp.psi(np.sum(alpha)))
    return(sp.psi(alpha) - sp.psi(np.sum(alpha, 1))[:, np.newaxis])

def expect_log_sticks(sticks):
    dig_sum = sp.psi(np.sum(sticks, 0))
    ElogW = sp.psi(sticks[0]) - dig_sum
    Elog1_W = sp.psi(sticks[1]) - dig_sum

    n = len(sticks[0]) + 1
    Elogsticks = np.zeros(n)
    Elogsticks[0:n-1] = ElogW
    Elogsticks[1:] = Elogsticks[1:] + np.cumsum(Elog1_W)
    return Elogsticks 

def log_normalize(v):
    ''' return log(sum(exp(v)))'''
    log_max = 100.0
    if len(v.shape) == 1:
        max_val = np.max(v)
        log_shift = log_max - np.log(len(v)+1.0) - max_val
        tot = np.sum(np.exp(v + log_shift))
        log_norm = np.log(tot) - log_shift
        v = v - log_norm
    else:
        max_val = np.max(v, 1)
        log_shift = log_max - np.log(v.shape[1]+1.0) - max_val
        tot = np.sum(np.exp(v + log_shift[:,np.newaxis]), 1)

        log_norm = np.log(tot) - log_shift
        v = v - log_norm[:,np.newaxis]
    return (v, log_norm)


def log_sum(log_a, log_b):
	''' we know log(a) and log(b), compute log(a+b) '''
	v = 0.0;
	if (log_a < log_b):
		v = log_b+np.log(1 + np.exp(log_a-log_b))
	else:
		v = log_a+np.log(1 + np.exp(log_b-log_a))
	return v


def argmax(x):
	''' find the index of maximum value '''
	n = len(x)
	val_max = x[0]
	idx_max = 0
	for i in range(1, n):
		if x[i]>val_max:
			val_max = x[i]
			idx_max = i		
	return idx_max	


def converttoliblinear(xx):
        prob_x = []
	for s in xx:
            xi = {}
            i = 1
            for val in s:
               xi[i] = float(val)
               i = i + 1 
            prob_x += [xi]
	return (prob_x)

def read_1D_array(data_file_name):	
	prob_x = []
	for line in open(data_file_name):
		for e in line.split():
			prob_x += [int(e)]
	return (prob_x)


def read_2D_array(data_file_name):		
        prob_y = [] 
	for line in open(data_file_name):
                prob_x = []
                line = line.rstrip()
		for e in line.split(" "):
                        prob_x += [(float(e))]
                prob_y += [prob_x]
	return (prob_y)	

        """f = file("idx.txt", "w")
        for x in idx:   
               f.write(str(x) + '\n')
        f.close()"""
def writefiles(x):
        f = file("Elogbeta.txt", "w")
        for xx in Elogbeta:
               line = ' '.join(str(x) for x in xx)   
               f.write(line + '\n')
        f.close()

        """f = file("lambda.txt", "w") 
        for beta in self.m_beta:
            line = ' '.join([str(x) for x in beta])  
            f.write(line + '\n')
        f.close()"""
	
