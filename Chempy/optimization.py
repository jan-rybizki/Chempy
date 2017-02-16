import numpy as np
from cem_function import cem
import multiprocessing as mp
import os
from parameter import ModelParameters

def one_chain(args):
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

def prior_sampling_single(args):
	seed,a = args
	np.random.seed(seed)
	chain = np.zeros(shape=(a.ndim))
	result = -np.inf
	while result == -np.inf:
		for j,name in enumerate(a.to_optimize):
			(mean, std) = a.priors.get(name)
			assert std > 0
			(lower, upper) = a.constraints.get(name)
			value = mean + np.random.normal(0,std)
			while value < lower and lower is not None or value > upper and upper is not None:
				value = mean + np.random.normal(0,std)
			chain[j] = value 
		result = posterior_probability(chain,a)
	return chain

def prior_sampling(args):
	seed,a = args
	np.random.seed(seed[0])
	chain = np.zeros(shape=(a.ndim))
	result = -np.inf
	while result == -np.inf:
		for j,name in enumerate(a.to_optimize):
			(mean, std) = a.priors.get(name)
			assert std > 0
			(lower, upper) = a.constraints.get(name)
			value = mean + np.random.normal(0,std)
			while value < lower and lower is not None or value > upper and upper is not None:
				value = mean + np.random.normal(0,std)
			chain[j] = value 
		result = posterior_probability(chain,a)
	return chain

def creating_chain(a,startpoint):
	args = [(np.random.random_integers(low = 0, high = 1e9,size = 1),a,startpoint) for i in range(a.nwalkers)]
	p = mp.Pool(a.nwalkers)
	t = p.map(one_chain,args)
	p.close()
	p.join()
	return np.array(t)


def gaussian_log(x,x0,xsig):
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)

def posterior_probability(x,a):
	s,t = cem(x,a)
	return s,t