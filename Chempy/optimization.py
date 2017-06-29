import numpy as np

import multiprocessing as mp
import os
from .parameter import ModelParameters
from .cem_function import posterior_function

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
			# new formulation due to problem in python 3
			if lower is None and upper is not None:
				while value > upper:
					value = mean + np.random.normal(0,std)
			elif upper is None and lower is not None:
				while value < lower:
					value = mean + np.random.normal(0,std)
			elif lower is not None and upper is not None:
				while value < lower or value > upper:
					value = mean + np.random.normal(0,std)
			### alternative in python 3
			#while value < lower and lower is not None or value > upper and upper is not None:
			#	value = mean + np.random.normal(0,std)
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
	s,t = posterior_function(x,a)#cem(x,a)#posterior_function(x,a)
	return s,t

def minimizer_initial(identifier):
	'''
	This is a function that is minimizing the posterior of Chempy from initial conditions (and can be called in multiprocesses)

	INPUT:

	   a = model parameters

	OUTPUT:

	   res.x = the free Chempy parameter for which the posterior was minimal (log posterior is maximized)
	'''
	from scipy.optimize import minimize
	from .cem_function import posterior_function_for_minimization
	from .parameter import ModelParameters

	a = ModelParameters()
	a.stellar_identifier = identifier

	res = minimize(fun = posterior_function_for_minimization,
		x0 = a.p0,
		args = (a),
		method = 'Nelder-Mead',
		tol = a.tol_minimization,
		options = {'maxiter':a.maxiter_minimization})
	if a.verbose:
		print(res.message)
	return res.x

def minimizer_local(args):
	from scipy.optimize import minimize
	from .cem_function import posterior_function_local_for_minimization
	from .parameter import ModelParameters

	a = ModelParameters()

	changing_parameter, identifier, global_parameters, errors, elements = args
	
	res = minimize(fun = posterior_function_local_for_minimization,
		x0 = changing_parameter,
		args = (identifier , global_parameters, errors, elements),
		method = 'Nelder-Mead',
		tol = a.tol_minimization,
		options = {'maxiter':a.maxiter_minimization})
	if a.verbose:
		print(res.message)
	return res.x

def minimizer_global(changing_parameter, tol, maxiter, verbose, result):
	'''
	This is a function that minimizes the posterior coming from global optimization

	INPUT:

	   changing_parameter = the global SSP parameters (parameters that all stars share)

	   tol = at which change in posterior the minimization should stop

	   maxiter = maximum number of iteration

	   verbose = print or print not result (bool)

	   result = the complete parameter set is handed over as an array of shape(len(stars),len(all parameters)). From those the local ISM parameters are taken

	OUTPUT:

	   rex.x = for which global parameters the minimization returned the best posterior
	'''
	from scipy.optimize import minimize
	from .cem_function import global_optimization
	
	res = minimize(fun = global_optimization,
		x0 = changing_parameter,
		args = (result),
		method = 'Nelder-Mead',
		tol = tol,
		options = {'maxiter':maxiter})
	if verbose:
		print(res.message)
	return res.x