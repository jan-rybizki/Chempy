import numpy as np
from cem_function import cem
import multiprocessing as mp
import os
from parameter import ModelParameters

def one_chain(args):
	'''
        This function is testing the startpoint of an MCMC chain and tries to find parameter values for which Chempy does not return -inf 
        '''
        seed,a,startpoint = args
	np.random.seed(seed[0])
	chain = np.zeros(shape=(a.ndim))
	backup = np.array(a.to_optimize)
	result = -np.inf
	while result == -np.inf:
		for j,name in enumerate(a.to_optimize):
			#(mean, std) = a.priors.get(name)
			mean = startpoint[j]
			std = np.abs(0.02 * startpoint[j])
			if std == 0:
				std = 0.02
			(lower, upper) = a.constraints.get(name)
			value = mean + np.random.normal(0,std)
			while value < lower and lower is not None or value > upper and upper is not None:
				value = mean + np.random.normal(0,std)
			chain[j] = value 
		a.to_optimize = backup
		result = posterior_probability(chain,a)
	return chain

def creating_chain(a,startpoint):
	'''
        This function creates the initial parameter values for an MCMC chain.

        INPUT:
        a = default parameter values from parameter.py
        startpoint = from where the pointcloud of walkers start in a small sphere

        OUTPUT:
        returns the array of the initial startpoints
        '''
        args = [(np.random.random_integers(low = 0, high = 1e9,size = 1),a,startpoint) for i in range(a.nwalkers)]
	p = mp.Pool(a.nwalkers)
	t = p.map(one_chain,args)
	p.close()
	p.join()
	return np.array(t)


def gaussian_log(x,x0,xsig):
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)

def posterior_probability(x,a):
        '''
        Just returning the posterior probability of Chempy and the list of blobs
        '''
	s,t = cem(x,a)
	return s,t
