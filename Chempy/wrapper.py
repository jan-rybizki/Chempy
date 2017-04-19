import numpy as np 
from .weighted_yield import SSP, lifetime_Argast, lifetime_Raiteri
from .imf import IMF
from .yields import SN2_feedback, AGB_feedback, SN1a_feedback, Hypernova_feedback


class SSP_wrap():
	'''
	This is the wrapper around the SSP function. It preloads the needed classes and calls all nucleosynthetic enrichment processes when the enrichment is calculated.
	'''
	def __init__(self, a):
		'''
		Upon initialization the default IMF, CC-SN yields, SN Ia yields and AGB yields is loaded.

		INPUT:
		
		   a = Modelparameter class. So the default IMF etc are loaded. If we want other yield sets etc. loaded we need to specify that in paramter.py
		'''

		## loading the IMF and the yieldsets prescribed in a (containing all the model parameters)
		basic_imf = IMF(a.mmin,a.mmax,a.mass_steps)
		getattr(basic_imf, a.imf_type_name)(a.imf_parameter)
		basic_sn2 = SN2_feedback()
		getattr(basic_sn2, a.yield_table_name_sn2)()
		basic_1a = SN1a_feedback()
		getattr(basic_1a, a.yield_table_name_1a)()
		basic_agb = AGB_feedback()
		getattr(basic_agb, a.yield_table_name_agb)()
		### mixing of Nomoto CC-SN and HN yields
		if a.yield_table_name_sn2 == 'Nomoto2013':
			basic_hn = Hypernova_feedback()
			getattr(basic_hn, a.yield_table_name_hn)()
			for item in basic_sn2.metallicities:
				x = np.copy(basic_sn2.table[item])
				y = np.copy(basic_hn.table[item])
				for jtem in basic_hn.masses:
					basic_sn2.table[item]['mass_in_remnants'][np.where(basic_sn2.table[item]['Mass']==jtem)] = a.sn2_to_hn * (x['mass_in_remnants'][np.where(x['Mass']==jtem)]) + (1-a.sn2_to_hn) * (y['mass_in_remnants'][np.where(y['Mass']==jtem)])
					basic_sn2.table[item]['unprocessed_mass_in_winds'][np.where(basic_sn2.table[item]['Mass']==jtem)] = a.sn2_to_hn * (x['unprocessed_mass_in_winds'][np.where(x['Mass']==jtem)]) + (1-a.sn2_to_hn) * (y['unprocessed_mass_in_winds'][np.where(y['Mass']==jtem)])
					hn_mass = []
					sn_mass = []
					for stem in basic_sn2.elements:
						sn_mass.append(x[stem][np.where(x['Mass']==jtem)])
						hn_mass.append(y[stem][np.where(y['Mass']==jtem)])
						basic_sn2.table[item][stem][np.where(basic_sn2.table[item]['Mass']==jtem)]= a.sn2_to_hn * (x[stem][np.where(x['Mass']==jtem)]) + (1-a.sn2_to_hn) * (y[stem][np.where(y['Mass']==jtem)])
		## to pass the information on to the feedback calculation
		self.a = a
		self.imf = basic_imf
		self.sn2 = basic_sn2
		self.sn1a = basic_1a
		self.agb = basic_agb

	def calculate_feedback(self, z, elements, element_fractions, time_steps):
		'''
		The feedback is calculated for the initializes SSP.

		INPUT:
		
		   z = metallicity of the SSP in mass fraction (not normed to solar!)
		
		   elements = which elements to follow
		
		   element_fractions = the birth material of the SSP in the same order as 'elements'
		
		   time_steps = the time-steps for which the enrichment of the SSP should be calculated (usually the time-steps until the end of the chempy simulation)
		'''
		basic_ssp = SSP(False, float(z), np.copy(self.imf.x), np.copy(self.imf.dm), np.copy(self.imf.dn), np.copy(time_steps), list(elements), str(self.a.stellar_lifetimes), str(self.a.interpolation_scheme), bool(self.a.only_net_yields_in_process_tables))
		basic_ssp.sn2_feedback(list(self.sn2.elements), dict(self.sn2.table), np.copy(self.sn2.metallicities), float(self.a.sn2mmin), float(self.a.sn2mmax),list(element_fractions))
		basic_ssp.agb_feedback(list(self.agb.elements), dict(self.agb.table), list(self.agb.metallicities), float(self.a.agbmmin), float(self.a.agbmmax),np.hstack(element_fractions))
		basic_ssp.sn1a_feedback(list(self.sn1a.elements), list(self.sn1a.metallicities), dict(self.sn1a.table), str(self.a.time_delay_functional_form), float(self.a.sn1ammin), float(self.a.sn1ammax), self.a.sn1a_parameter, float(self.a.total_mass), bool(self.a.stochastic_IMF))
		# exposing these tables to the outside wrapper
		self.table = basic_ssp.table
		self.sn2_table = basic_ssp.sn2_table
		self.agb_table = basic_ssp.agb_table
		self.sn1a_table = basic_ssp.sn1a_table
		self.inverse_imf = basic_ssp.inverse_imf

def initialise_stuff(a):
	'''
	Convenience function initialising the solar abundance, SFR and infall with the default values provided in parameter.py as a
	'''
	from .solar_abundance import solar_abundances
	from .sfr import SFR 
	from .infall import INFALL

	basic_solar = solar_abundances()
	getattr(basic_solar, a.solar_abundance_name)()

	basic_sfr = SFR(a.start,a.end,a.time_steps)
	if a.basic_sfr_name == 'gamma_function':
		getattr(basic_sfr, a.basic_sfr_name)(S0 = a.S_0 * a.mass_factor,a_parameter = a.a_parameter, loc = a.sfr_beginning, scale = a.sfr_scale)
	elif a.basic_sfr_name == 'model_A':
		basic_sfr.model_A(a.mass_factor*a.S_0,a.t_0,a.t_1)
	elif a.basic_sfr_name == 'prescribed':
		basic_sfr.prescribed(a.mass_factor, a.name_of_file)
	elif a.basic_sfr_name == 'doubly_peaked':
		basic_sfr.doubly_peaked(S0 = a.mass_factor*a.S_0, peak_ratio = a.peak_ratio, decay = a.sfr_decay, t0 = a.sfr_t0, peak1t0 = a.peak1t0, peak1sigma = a.peak1sigma)
	basic_sfr.sfr = a.total_mass * np.divide(basic_sfr.sfr,sum(basic_sfr.sfr))

	basic_infall = INFALL(np.copy(basic_sfr.t),np.copy(basic_sfr.sfr))
	if a.basic_infall_name == 'exponential':
		getattr(basic_infall, a.basic_infall_name)((a.infall_amplitude,a.tau_infall,a.infall_time_offset,a.c_infall,a.norm_infall))
	elif a.basic_infall_name == 'gamma_function':
		getattr(basic_infall, a.basic_infall_name)(mass_factor = a.norm_infall, a_parameter = a.infall_a_parameter, loc = a.infall_beginning, scale = a.infall_scale)
	elif a.basic_infall_name == 'sfr_related':
		getattr(basic_infall, a.basic_infall_name)()


	return basic_solar, basic_sfr, basic_infall

def Chempy(a):
	'''
	Chemical evolution run with the default parameters using the net yields.

	INPUT: 
	
	   a = ModelParameters() from parameter.py

	OUTPUT:
	
	   cube = The ISM evolution class
	
	   abundances = The abundances of the ISM
	'''
	from .infall import PRIMORDIAL_INFALL
	from .time_integration import ABUNDANCE_MATRIX
	from .making_abundances import mass_fraction_to_abundances
	from numpy.lib.recfunctions import append_fields	
	basic_solar, basic_sfr, basic_infall = initialise_stuff(a)
	elements_to_trace = a.elements_to_trace
	basic_primordial = PRIMORDIAL_INFALL(list(elements_to_trace),np.copy(basic_solar.table))
	basic_primordial.primordial()
	cube = ABUNDANCE_MATRIX(np.copy(basic_sfr.t),np.copy(basic_sfr.sfr),np.copy(basic_infall.infall),list(elements_to_trace),list(basic_primordial.symbols),list(basic_primordial.fractions),float(a.gas_at_start),list(basic_primordial.symbols),list(basic_primordial.fractions),float(a.gas_reservoir_mass_factor),float(a.outflow_feedback_fraction),bool(a.check_processes),float(a.starformation_efficiency),float(a.gas_power), float(a.sfr_factor_for_cosmic_accretion), list(basic_primordial.symbols), list(basic_primordial.fractions))
	basic_ssp = SSP_wrap(a)
	for i in range(len(basic_sfr.t)-1):
		j = len(basic_sfr.t)-i
		element_fractions = []
		for item in elements_to_trace:
			element_fractions.append(float(np.copy(cube.cube[item][max(i-1,0)]/cube.cube['gas'][max(i-1,0)])))## gas element fractions from one time step before	
		metallicity = float(cube.cube['Z'][i])
		time_steps = np.copy(basic_sfr.t[:j])
		basic_ssp.calculate_feedback(float(metallicity), list(elements_to_trace), list(element_fractions), np.copy(time_steps))
		cube.advance_one_step(i+1,np.copy(basic_ssp.table),np.copy(basic_ssp.sn2_table),np.copy(basic_ssp.agb_table),np.copy(basic_ssp.sn1a_table))
		if cube.cube['gas'][i] < 0:
			print(i, basic_sfr.t[i])
			print('gas became negative. returning -inf')
			return -np.inf, [0]
		if cube.gas_reservoir['gas'][i] < 0:
			print('gas_reservoir became negative. returning -inf')
			return -np.inf, [0]

	abundances,elements,numbers = mass_fraction_to_abundances(np.copy(cube.cube),np.copy(basic_solar.table))
	weights = cube.cube['sfr']
	abundances = append_fields(abundances,'weights',weights)
	abundances = append_fields(abundances,'time', cube.cube['time'])
	abundances = np.array(abundances)

	return cube, abundances




def Chempy_gross(a):
	'''
	Chemical evolution run with the default parameters but now using solar scaled material (testing the worse case when total yields provided).

	INPUT: 
	
	   a = ModelParameters() from parameter.py

	OUTPUT:
	
	   cube = The ISM evolution class
	
	   abundances = The abundances of the ISM
	'''
	from infall import PRIMORDIAL_INFALL
	from time_integration import ABUNDANCE_MATRIX
	from making_abundances import mass_fraction_to_abundances
	from numpy.lib.recfunctions import append_fields	
	basic_solar, basic_sfr, basic_infall = initialise_stuff(a)
	elements_to_trace = a.elements_to_trace
	basic_primordial = PRIMORDIAL_INFALL(list(elements_to_trace),np.copy(basic_solar.table))
	basic_primordial.primordial(0)
	cube = ABUNDANCE_MATRIX(np.copy(basic_sfr.t),np.copy(basic_sfr.sfr),np.copy(basic_infall.infall),list(elements_to_trace),list(basic_primordial.symbols),list(basic_primordial.fractions),float(a.gas_at_start),list(basic_primordial.symbols),list(basic_primordial.fractions),float(a.gas_reservoir_mass_factor),float(a.outflow_feedback_fraction),bool(a.check_processes),float(a.starformation_efficiency),float(a.gas_power), float(a.sfr_factor_for_cosmic_accretion), list(basic_primordial.symbols), list(basic_primordial.fractions))
	basic_ssp = SSP_wrap(a)
	for i in range(len(basic_sfr.t)-1):
		j = len(basic_sfr.t)-i
		metallicity = float(cube.cube['Z'][i])
		solar_scaled_material = PRIMORDIAL_INFALL(list(elements_to_trace),np.copy(basic_solar.table))
		solar_scaled_material.solar(np.log10(metallicity/basic_solar.z))
		element_fractions = list(solar_scaled_material.fractions)
		for item in elements_to_trace:
			element_fractions.append(float(np.copy(cube.cube[item][max(i-1,0)]/cube.cube['gas'][max(i-1,0)])))## gas element fractions from one time step before	
		time_steps = np.copy(basic_sfr.t[:j])
		basic_ssp.calculate_feedback(float(metallicity), list(elements_to_trace), list(element_fractions), np.copy(time_steps))
		cube.advance_one_step(i+1,np.copy(basic_ssp.table),np.copy(basic_ssp.sn2_table),np.copy(basic_ssp.agb_table),np.copy(basic_ssp.sn1a_table))
	abundances,elements,numbers = mass_fraction_to_abundances(np.copy(cube.cube),np.copy(basic_solar.table))
	weights = cube.cube['sfr']
	abundances = append_fields(abundances,'weights',weights)
	abundances = np.array(abundances)

	return cube, abundances


def multi_star_optimization():
	'''
	This function will optimize the parameters of all stars in a hierachical manner (similar to gibbs sampling)

	INPUT: 

	   a = will be loaded from parameter.py (prepare all variables there)

	OUTPUT:

	   log_list = a list of intermediate results (so far only for debugging)
	'''
	import time
	import multiprocessing as mp
	from .optimization import minimizer_initial, minimizer_global, minimizer_local
	from .cem_function import global_optimization_error_returned
	from .parameter import ModelParameters
	
	a = ModelParameters()

	start_time = time.time()

	log_list = []
	# I: Minimization for each star seperately
	# 1: for each star make initial conditions (each star needs other model parameters)	
	parameter_list = []
	for item in a.stellar_identifier_list:
		parameter_list.append(item)
	# 2: call posterior_function_for_minimization with scipy.optimize.minimize in multiprocess for each star and recover the found parameters
	p = mp.Pool(len(parameter_list))
	t = p.map(minimizer_initial, parameter_list)
	p.close()
	p.join()
	result = np.vstack(t)

	log_list.append(np.copy(result))
	log_list.append('initial minimization')
	initial = time.time()
	print('first minimization for each star separately took: %2.f seconds' %(initial - start_time))

	# IV: repeat II and III until posterior does not change much
	result[:,:len(a.SSP_parameters)] = np.mean(result[:,:len(a.SSP_parameters)], axis = 0)
	posteriors = []
	while True:
		if len(posteriors) > 1:
			if np.abs(posteriors[-1] - posteriors[-2]) < a.gibbs_sampler_tolerance:
				break

		initial = time.time()
		# II: Global parameter minimization:
		# 1: only SSP parameters free. Use mean SSP parameter values and individual (but fixed ISM parameter values)
		changing_parameter = result[0,:len(a.SSP_parameters)]
		# 2: Call each star in multiprocess but only return the predictions
		# 3: Calculate the likelihood for each star and optimize the common model error (is all done within minimizer global, which is calling 'global optimization')
		x = minimizer_global(changing_parameter,  a.tol_minimization, a.maxiter_minimization, a.verbose, result)

		# 4: return global SSP parameters and common model error
		posterior, error_list, elements = global_optimization_error_returned(x, result)
		posteriors.append(posterior)
		print(posteriors)

		global_iteration1 = time.time()
		print('first global minimization took: %2.f seconds' %(global_iteration1 - initial))	

		# III: Local parameter minimization:
		# 1: Use fixed global parameters and fixed common errors make initial conditions
		result[:,:len(a.SSP_parameters)] = x

		log_list.append((np.copy(x),posterior, error_list,elements))
		log_list.append('global minimization')

		p0_list = []
		parameter_list = []
		x_list = []
		error_list_mp = []
		element_list_mp = []

		for i,item in enumerate(a.stellar_identifier_list):
			parameter_list.append(item)
			p0_list.append(result[i,len(a.SSP_parameters):])
			x_list.append(x)
			error_list_mp.append(error_list)
			element_list_mp.append(elements)

		args = zip(p0_list,parameter_list,x_list,error_list_mp,element_list_mp)

		# 2: Minimize each star ISM parameters in multiprocess
		p = mp.Pool(len(parameter_list))
		t = p.map(minimizer_local, args)
		p.close()
		p.join()
		local_parameters = np.vstack(t)
		result[:,len(a.SSP_parameters):] = local_parameters

		log_list.append(np.copy(result))
		log_list.append('local minimization')
		local_iteration1 = time.time()
		print('first local minimization took: %2.f seconds' %(local_iteration1 - global_iteration1))	

	log_list.append(posteriors)
	print(log_list)

	# V: MCMC run
	## reshape the result to have global parameters in the front and the local parameters following
	changing_parameter = list(result[0,:len(a.SSP_parameters)])
	for i in range(result.shape[0]):
		changing_parameter.append(list(result[i,len(a.SSP_parameters):]))
	changing_parameter = np.hstack(changing_parameter)
	## jitter the parameters to initialise the chain (add a validation later, i.e. testing that the particular parameters yield a result)
	mcmc_multi(changing_parameter, error_list, elements)
	# 1: Free all parameters and optimize common error (SSP should be the same for all stars)
	# 2: Plug everything into emcee and sample the posterior
	return log_list

def mcmc(a):
	'''
	Convenience function to use the MCMC. A subdirectory mcmc/ will be created in the current directory and intermediate chains will be stored there.
	
	The MCMC will sample the volume of best posterior for the likelihood functions that are declared in parameter.py. Default is ['sol_norm','gas_reservoir','sn_ratio'] which corresponds to 'Sun+' from the paper.
	'''
	import time
	import os
	import multiprocessing as mp
	from .optimization import creating_chain, posterior_probability
	import emcee

	start1 = time.time()
	directory = 'mcmc/'
	if os.path.exists(directory):
		if a.verbose:
			print('%s already existed. Content might be overwritten' %(directory))
	else:
		os.makedirs(directory)
	
	a.check_processes = False
	a.number_of_models_overplotted = 1
	a.only_net_yields_in_process_tables = False
	a.testing_output = False
	a.summary_pdf = False
	a.nthreads = mp.cpu_count()
	if a.nthreads == 4:
		a.nthreads = 2
	
	chain = creating_chain(a,np.copy(a.p0))
	sampler = emcee.EnsembleSampler(a.nwalkers,a.ndim,posterior_probability,threads=a.nthreads, args = [a])
	pos,prob,state,blobs = sampler.run_mcmc(chain,a.mburn)
	
	mean_prob = mean_prob_beginning = np.zeros((a.m))
	posterior_list = []
	posterior_std_list = []
	for i in range(a.m):
		print('step ', i+1 , 'of ',a.m)
		pos, prob, state, blobs = sampler.run_mcmc(pos, a.save_state_every, rstate0=state, lnprob0=prob, blobs0 = blobs, storechain = True)
		np.save('%s/flatchain' %(directory),sampler.chain)
		np.save('%s/flatlnprobability' %(directory),sampler.lnprobability)
		np.save('%s/flatblobs' %(directory),sampler.blobs)
		posterior = np.load('%s/flatlnprobability.npy' %(directory))
		posterior_list.append(np.mean(posterior, axis = 0)[-1])
		posterior_std_list.append(np.std(posterior, axis = 0)[-1])
		np.save('%s/flatmeanposterior' %(directory), posterior_list)
		np.save('%s/flatstdposterior' %(directory), posterior_std_list)
		print(np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[-1])
		
		if i>202:
			print('posterior -1, -100, -200',np.mean(posterior, axis = 0)[-1], np.mean(posterior, axis = 0)[-100], np.mean(posterior, axis = 0)[-200])
			print('posterior 0, 100, 200',np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[100], np.mean(posterior, axis = 0)[200])
		#print("Mean acceptance fraction:", sampler.acceptance_fraction)
		elapsed1 = (time.time() - start1)
		print('calculation so far took', elapsed1, ' seconds')
		if i>300 and np.abs(np.mean(posterior, axis = 0)[-1] - np.mean(posterior, axis = 0)[-100]) < 0.5 and np.abs(np.mean(posterior, axis = 0)[-1] - np.mean(posterior, axis = 0)[-200]) < 0.5:
			break


def mcmc_multi(changing_parameter, error_list, elements):
	'''
	Convenience function to use the MCMC for multiple zones (and therefore multiple observations). A subdirectory mcmc/ will be created in the current directory and intermediate chains will be stored there.
	The MCMC will sample the volume of best posterior for the likelihood functions that are declared in parameter.py. 
	Default is a list of Proto-sun, Arcturus and B-stars. The MCMC uses many walkers and can use multiple threads. Each walker will evaluate a series of Chempy zones and add their posterior together which then will be returned.
	
	INPUT:

	   changing_parameter = the parameter vector for initialization (will usually be found from minimization before). The initial chain will be created by jittering slightly the initial parameter guess

	   error_list = the vector of element errors

	   elements = the corresponding element symbols

	OUTPUT:

	   The function will create a folder and store the chain as well as the predicted element values

	The MCMC stops when the convergence criteria is met, which is when the median posterior of all walkers does not change much inbetween 200 steps anymore.
	'''
	import time
	import os
	import multiprocessing as mp
	from .cem_function import  posterior_function_many_stars
	from .parameter import ModelParameters
	import emcee

	a = ModelParameters()
	start1 = time.time()
	directory = 'mcmc/'
	if os.path.exists(directory):
		if a.verbose:
			print('%s already existed. Content might be overwritten' %(directory))
	else:
		os.makedirs(directory)
	
	nthreads = mp.cpu_count()
	if nthreads == 4:
		nthreads = 2
	ndim = len(changing_parameter)
	a.nwalkers = max(a.nwalkers, int(ndim*2))
	chain = np.empty(shape = (a.nwalkers,ndim))
	for i in range(a.nwalkers):
		result = -np.inf
		while result == -np.inf:
			jitter = np.random.normal(loc = 0, scale = 0.001, size = ndim)
			result, dummy = posterior_function_many_stars(changing_parameter + jitter,error_list,elements)
		chain[i] = changing_parameter + jitter

	sampler = emcee.EnsembleSampler(a.nwalkers,ndim,posterior_function_many_stars,threads=nthreads, args = [error_list,elements])
	pos,prob,state,blobs = sampler.run_mcmc(chain,a.mburn)
	
	mean_prob = mean_prob_beginning = np.zeros((a.m))
	posterior_list = []
	posterior_std_list = []
	for i in range(a.m):
		print('step ', i+1 , 'of ',a.m)
		pos, prob, state, blobs = sampler.run_mcmc(pos, a.save_state_every, rstate0=state, lnprob0=prob, blobs0 = blobs, storechain = True)
		np.save('%s/flatchain' %(directory),sampler.chain)
		np.save('%s/flatlnprobability' %(directory),sampler.lnprobability)
		np.save('%s/flatblobs' %(directory),sampler.blobs)
		posterior = np.load('%s/flatlnprobability.npy' %(directory))
		posterior_list.append(np.mean(posterior, axis = 0)[-1])
		posterior_std_list.append(np.std(posterior, axis = 0)[-1])
		np.save('%s/flatmeanposterior' %(directory), posterior_list)
		np.save('%s/flatstdposterior' %(directory), posterior_std_list)
		print(np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[-1])
		
		if i>202:
			print('posterior -1, -100, -200',np.mean(posterior, axis = 0)[-1], np.mean(posterior, axis = 0)[-100], np.mean(posterior, axis = 0)[-200])
			print('posterior 0, 100, 200',np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[100], np.mean(posterior, axis = 0)[200])
		#print("Mean acceptance fraction:", sampler.acceptance_fraction)
		elapsed1 = (time.time() - start1)
		print('calculation so far took', elapsed1, ' seconds')
		if i>a.min_mcmc_iterations and np.abs(np.mean(posterior, axis = 0)[-1] - np.mean(posterior, axis = 0)[-100]) < 0.5 and np.abs(np.mean(posterior, axis = 0)[-1] - np.mean(posterior, axis = 0)[-200]) < 0.5:
			break