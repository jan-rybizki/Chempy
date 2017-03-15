import numpy as np
import os
from .sfr import SFR
from .solar_abundance import solar_abundances
import time
from .data_to_test import likelihood_function, wildcard_likelihood_function, elements_plot, arcturus, sol_norm, plot_processes, save_abundances,  cosmic_abundance_standard, ratio_function, star_function, gas_reservoir_metallicity
import multiprocessing as mp
from .wrapper import initialise_stuff, Chempy

def gaussian_log(x,x0,xsig):
	'''
	function to calculate the gaussian probability (its normed to Pmax and given in log)
	
	INPUT:
	
	   x = where is the data point or parameter value
	
	   x0 = mu
	
	   xsig = sigma
	'''
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)

def lognorm_log(x,mu,factor):
	'''
	this function provides Prior probability distribution where the factor away from the mean behaves like the sigma deviation in normal_log 
	
	for example if mu = 1 and factor = 2 
	
	for	1 it returns 0
	
	for 0,5 and 2 it returns -0.5
	
	for 0.25 and 4 it returns -2.0
	
	and so forth
	
	Can be used to specify the prior on the yield factors
	'''
	y = np.log(np.divide(x,mu))
	y = np.divide(y,np.log(factor))
	y = gaussian_log(y,0.,1.)
	return y

def gaussian(x,x0,xsig):
	'''
	function to calculate the gaussian probability (its normed to Pmax and given in log)
	
	INPUT:
	
	   x = where is the data point or parameter value
	
	   x0 = mu
	
	   xsig = sigma
	'''
	factor = 1. / (np.sqrt(xsig * xsig * 2. * np.pi))
	exponent = -np.divide((x - x0) * (x - x0),2 * xsig * xsig)
	return factor * np.exp(exponent)

def lognorm(x,mu,factor):
	'''
	this function provides Prior probability distribution where the factor away from the mean behaves like the sigma deviation in normal_log 
	BEWARE: this function is not a properly normalized probability distribution. It only provides relative values.
	
	INPUT:

	   x = where to evaluate the function, can be an array

	   mu = peak of the distribution

	   factor = the factor at which the probability decreases to 1 sigma
	
	Can be used to specify the prior on the yield factors
	'''
	y = np.log(np.divide(x,mu))
	y = np.divide(y,np.log(factor))
	y = gaussian(y,0.,1.)
	return y

def shorten_sfr(a):
	'''
	This function crops the SFR to the length of the age of the star and ensures that enough stars are formed at the stellar birth epoch

	INPUT:

	   a = Modelparameters

	OUTPUT:
	
	   the function will update the modelparameters, such that the simulation will end when the star is born and it will also check whether there is enough sfr left at that epoch
	'''
	try:
		star = np.load('%s.npy' %(a.stellar_identifier))
	except Exception as ex:
		from . import localpath
		star = np.load(localpath + 'input/stars/' + a.stellar_identifier + '.npy')
	age_of_star = star['age'][0]
	assert (age_of_star <= 13.0), "Age of the star must be below 13Gyr"

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
	mass_normalisation = a.total_mass
	mean_sfr = sum(basic_sfr.sfr) / a.end
	
	# at which time in the simulation is the star born
	star_time = basic_sfr.t[-1] - age_of_star
	cut = [np.where(np.abs(basic_sfr.t - star_time) == np.min(np.abs(basic_sfr.t - star_time)))]
	if len(cut[0][0]) != 1:
		cut = cut[0][0][0]
	
	# updating the end time and the model steps and rescale the total mass
	time_model = float(basic_sfr.t[cut])
	a.end = time_model
	a.time_steps = int(cut[0][0]) + 1
	a.total_mass = sum(basic_sfr.sfr[0:a.time_steps])

	# check whether the sfr is enough at end to produce reasonable number of stars (which is necessary in order to have a probability to observe a star at all)
	sfr_at_end = float(basic_sfr.sfr[cut] / basic_sfr.dt)
	fraction_of_mean_sfr = sfr_at_end / mean_sfr 	
	assert fraction_of_mean_sfr > 0.05, 'The total SFR of the last age bin is below 5% of the mean SFR'
	return a

def cem(changing_parameter,a):
	'''
	This is the function calculating the chemical evolution for a specific parameter set (changing_parameter) and for a specific observational constraint specified in a (e.g. 'solar_norm' calculates the likelihood of solar abundances coming out of the model). It returns the posterior and a list of blobs. It can be used by an MCMC.
	This function actually encapsulates the real cem function in order to capture exceptions and in that case return -inf. This makes the MCMC runs much more stable

	INPUT: 
	
	   changing_parameter = parameter values of the free parameters as an array
	
	   a = model parameters specified in parameter.py. There are also the names of free parameters specified here

	OUTPUT:
	
	   log posterior, array of blobs
	
	   the blobs contain the prior values, the likelihoods and the actual values of each predicted data point (e.g. elemental abundance value)
	'''
	try:
		posterior, blobs = cem_real(changing_parameter,a)
		return posterior, blobs
	except Exception as ex:
		import traceback; traceback.print_exc()
	return -np.inf, [0]

def cem_real(changing_parameter,a):
	'''
	real chempy function. description can be found in cem
	'''
	for i,item in enumerate(a.to_optimize):
		setattr(a, item, changing_parameter[i])
		val = getattr(a, item)

	start_time = time.time()
	### PRIOR calculation, values are stored in parameter.py
	prior_names = []
	prior = []
	for name in a.to_optimize:
		(mean, std, functional_form) = a.priors.get(name)
		val = getattr(a, name)
		prior_names.append(name)
		if functional_form == 0:
			prior.append(gaussian_log(val, mean, std))
		elif functional_form == 1:
			prior.append(lognorm_log(val, mean, std))
	a.prior = prior

	for name in a.to_optimize:
		(lower, upper) = a.constraints.get(name)
		val = getattr(a, name)
		if lower is not None and val<lower:
			print('%s lower border is violated' %(name))
			return -np.inf, [0]
		if upper is not None and val>upper:
			print('%s upper border is violated' %(name))
			return -np.inf, [0]

	if not a.testing_output:
		print(changing_parameter,mp.current_process()._identity[0])#,a.observational_constraints_index
	else:
		print(changing_parameter)
	
	### So that the parameter can be plotted in linear space
	if 'log10_N_0' in a.to_optimize:
		a.N_0 = np.power(10,a.log10_N_0)
	if 'log10_sn1a_time_delay' in a.to_optimize:
		a.sn1a_time_delay = np.power(10,a.log10_sn1a_time_delay)
	if 'log10_starformation_efficiency' in a.to_optimize:
		a.starformation_efficiency = np.power(10,a.log10_starformation_efficiency)
	if 'log10_gas_reservoir_mass_factor' in a.to_optimize:
		a.gas_reservoir_mass_factor = np.power(10,a.log10_gas_reservoir_mass_factor)
	if 'log10_sfr_scale' in a.to_optimize:
		a.sfr_scale = np.power(10,a.log10_sfr_scale)

	if a.imf_type_name == 'salpeter':
		a.imf_parameter = (a.high_mass_slope)
	elif a.imf_type_name == 'Chabrier_2':
		a.imf_parameter = (a.chabrier_para1, a.chabrier_para2, a.chabrier_para3,a.high_mass_slope)
	elif a.imf_type_name == 'normed_3slope':	
		a.imf_parameter = (a.imf_slope_1,a.imf_slope_2,a.high_mass_slope,a.imf_break_1,a.imf_break_2)
	if a.time_delay_functional_form == 'maoz':
		a.sn1a_parameter = [a.N_0,a.sn1a_time_delay,a.sn1a_exponent,a.dummy]
	elif a.time_delay_functional_form == 'normal':
		a.sn1a_parameter = [a.number_of_pn_exlopding,a.sn1a_time_delay,a.sn1a_timescale,a.sn1a_gauss_beginning]
	elif a.time_delay_functional_form == 'gamma_function':
		a.sn1a_parameter = [a.sn1a_norm,a.sn1a_a_parameter,a.sn1a_beginning,a.sn1a_scale]

	basic_solar = solar_abundances()
	getattr(basic_solar, a.solar_abundance_name)()
	elements_to_trace = a.elements_to_trace
        
	directory = 'model_temp/'
	### Model is calculated
	if a.calculate_model:
		cube, abundances = Chempy(a)
		cube1 = cube.cube
		gas_reservoir = cube.gas_reservoir
		if a.testing_output:
			if os.path.exists(directory):
				print(directory, ' already exists. Content might be overwritten')
			else:
				os.makedirs(directory)
			np.save(directory + '%s_elements_to_trace' %(a.name_string), elements_to_trace)
			np.save(directory + '%s_gas_reservoir' %(a.name_string),gas_reservoir)
			np.save(directory + '%s_cube' %(a.name_string),cube1)
			np.save(directory + '%s_abundances' %(a.name_string),abundances)
	else:
		cube1 = np.load(directory + '%s_cube.npy' %(a.name_string))
		abundances = np.load(directory + '%s_abundances.npy' %(a.name_string))
		gas_reservoir = np.load(directory + '%s_gas_reservoir.npy' %(a.name_string))
		elements_to_trace = np.load(directory + '%s_elements_to_trace.npy' %(a.name_string))


	### LIKELIHOOD is being calculated
	a.probability = []
	a.abundance_list = []
	a.names = []
	### these functions need to return lists of the probabilities, likelihoods and elements names. The latter ones are important for the blobs so that the MCMC result can be enriched with elemental likelihoods and the like.
	if 'gas_reservoir' in a.observational_constraints_index:
		probabilities, result, names = gas_reservoir_metallicity(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,np.copy(gas_reservoir),a.number_of_models_overplotted,a.produce_mock_data,a.use_mock_data,a.error_inflation, np.copy(basic_solar.z))
		a.probability.append(probabilities)
		a.abundance_list.append(result)
		a.names.append(names)
	if 'sn_ratio' in a.observational_constraints_index:
		probabilities, result, names = ratio_function(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,np.copy(gas_reservoir),a.number_of_models_overplotted,a.produce_mock_data,a.use_mock_data,a.error_inflation)
		a.probability.append(probabilities)
		a.abundance_list.append(result)
		a.names.append(names)
	if 'cas' in a.observational_constraints_index:
		probabilities, abundance_list, element_names = cosmic_abundance_standard(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,np.copy(basic_solar.table),a.number_of_models_overplotted,a.produce_mock_data,a.use_mock_data,a.error_inflation)
		a.probability.append(probabilities)
		a.abundance_list.append(abundance_list)
		a.names.append(element_names)
	if 'sol_norm' in a.observational_constraints_index:
		probabilities, abundance_list, element_names = sol_norm(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,a.element_names,np.copy(basic_solar.table),a.number_of_models_overplotted,a.produce_mock_data,a.use_mock_data,a.error_inflation)
		a.probability.append(probabilities)
		a.abundance_list.append(abundance_list)
		a.names.append(element_names)
	if 'arcturus' in a.observational_constraints_index:
		probabilities, abundance_list, element_names = arcturus(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,a.element_names,np.copy(basic_solar.table),a.number_of_models_overplotted,a.arcturus_age,a.produce_mock_data,a.use_mock_data,a.error_inflation)
		a.probability.append(probabilities)
		a.abundance_list.append(abundance_list)
		a.names.append(element_names)
	if 'wildcard' in a.observational_constraints_index:
		probabilities, abundance_list, element_names = wildcard_likelihood_function(a.summary_pdf,a.stellar_identifier, np.copy(abundances))
		a.probability.append(probabilities)
		a.abundance_list.append(abundance_list)
		a.names.append(element_names)
	### These functions are for plotting and saving purposes they just return 0 probability not interfering with the likelihood. But they will make the MCMC blobs crash. Therefore take them out when running the MCMC
	if 'stars_at_end' in a.observational_constraints_index:
		a.probability.append(star_function(a.summary_pdf,a.name_string,np.copy(abundances),np.copy(cube1),elements_to_trace,np.copy(gas_reservoir),a.number_of_models_overplotted))
	if 'save_abundances' in a.observational_constraints_index:
		a.probability.append(save_abundances(a.summary_pdf,a.name_string,np.copy(abundances)))
	if 'plot_processes' in a.observational_constraints_index:
		a.probability.append(plot_processes(a.summary_pdf,a.name_string,cube.sn2_cube,cube.sn1a_cube,cube.agb_cube,a.element_names,np.copy(cube1),a.number_of_models_overplotted))
	if 'elements' in a.observational_constraints_index:
		a.probability.append(elements_plot(a.name_string,basic_ssp.agb.elements, basic_ssp.sn2.elements, basic_ssp.sn1a.elements,elements_to_trace, basic_solar.table,60))

	### to flatten the sublists so that the likelihood can be calculated and the blobs are in a flattened format
	a.names =  [item for sublist in a.names for item in sublist]
	a.names += ['m-%s' %(item) for item in a.names]
	a.abundance_list = [item for sublist in a.abundance_list for item in sublist]
	a.probability = [item for sublist in a.probability for item in sublist]
	a.names += prior_names
	if a.testing_output:
		#print a.names
		np.save("model_temp/blobs_name_list", a.names)
	if np.isnan(sum(a.probability)):
		return -np.inf, [0]
	if a.testing_output:
		print('l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability))#, mp.current_process()._identity[0]
	else:
		print('l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability),'|', mp.current_process()._identity[0])
	return sum(a.probability) + sum(a.prior), np.hstack((a.probability,a.abundance_list,a.prior))

def cem2(a):
	'''
	This is the function calculating the chemical evolution for a specific parameter set (changing_parameter) and for a specific observational constraint specified in a (e.g. 'solar_norm' calculates the likelihood of solar abundances coming out of the model). It returns the posterior and a list of blobs. It can be used by an MCMC.
	This function actually encapsulates the real cem function in order to capture exceptions and in that case return -inf. This makes the MCMC runs much more stable

	INPUT: 
	
	   a = model parameters specified in parameter.py and alteres by posterior_function

	OUTPUT:
	
	   predictions, name_of_prediction
	
	   the predicted element abundances for the time of the birth of the star (specified in a) are given back, as well as the corona metallicity at that time and the SN-ratio at that time.
	'''
	try:
		posterior, blobs = cem_real2(a)
		return posterior, blobs
	except Exception as ex:
		import traceback; traceback.print_exc()
	return -np.inf, [0]

def cem_real2(a):
	'''
	real chempy function. description can be found in cem2
	'''
	## The time until which Chempy is calculated is cropped to the stellar birth time. Also the SFR should not be below 1/20th of the mean SFR
	a = shorten_sfr(a)
	basic_solar = solar_abundances()
	getattr(basic_solar, a.solar_abundance_name)()
	elements_to_trace = list(a.elements_to_trace)
        
	directory = 'model_temp/'
	### Model is calculated
	if a.calculate_model:
		cube, abundances = Chempy(a)
		cube1 = cube.cube
		gas_reservoir = cube.gas_reservoir
		if a.testing_output:
			if os.path.exists(directory):
				print(directory, ' already exists. Content might be overwritten')
			else:
				os.makedirs(directory)
			np.save(directory + '%s_elements_to_trace' %(a.name_string), elements_to_trace)
			np.save(directory + '%s_gas_reservoir' %(a.name_string),gas_reservoir)
			np.save(directory + '%s_cube' %(a.name_string),cube1)
			np.save(directory + '%s_abundances' %(a.name_string),abundances)
	else:
		cube1 = np.load(directory + '%s_cube.npy' %(a.name_string))
		abundances = np.load(directory + '%s_abundances.npy' %(a.name_string))
		gas_reservoir = np.load(directory + '%s_gas_reservoir.npy' %(a.name_string))
		elements_to_trace = np.load(directory + '%s_elements_to_trace.npy' %(a.name_string))

	# predicted values are written out and returned together with corona metallicity and SN-ratio
	abundance_list = []
	for item in elements_to_trace:
		abundance_list.append(abundances[item][-1])
	
	abundance_list.append(gas_reservoir['Z'][-1])
	elements_to_trace.append('Zcorona')

	abundance_list.append(cube1['sn2'][-1]/cube1['sn1a'][-1])
	elements_to_trace.append('SNratio')

	return(abundance_list,elements_to_trace)


def posterior_function(changing_parameter,a):
	'''
	The posterior function is the interface between the optimizing function and Chempy. Usually the likelihood will be calculated with respect to a so called 'stellar wildcard'.
	Wildcards can be created according to the tutorial 6. A few wildcards are already stored in the input folder. Chempy will try the current folder first. If no wildcard npy file with the name a.stellar_identifier is found it will look into the Chempy/input/stars folder.

	INPUT: 
	
	   changing_parameter = parameter values of the free parameters as an array
	
	   a = model parameters specified in parameter.py. There are also the names of free parameters specified here

	OUTPUT:
	
	   log posterior, array of blobs
	
	   the blobs contain the likelihoods and the actual values of each predicted data point (e.g. elemental abundance value)
	'''
	try:
		posterior, blobs = posterior_function_real(changing_parameter,a)
		return posterior, blobs
	except Exception as ex:
		import traceback; traceback.print_exc()
	return -np.inf, [0]

def extract_parameters_and_priors(changing_parameter, a):
	'''
	This function extracts the parameters from changing parameters and writes them into the ModelParamaters (a), so that Chempy can evaluate the changed parameter settings
	'''
	for i,item in enumerate(a.to_optimize):
		setattr(a, item, changing_parameter[i])
		val = getattr(a, item)

	start_time = time.time()
	### PRIOR calculation, values are stored in parameter.py
	prior_names = []
	prior = []
	for name in a.to_optimize:
		(mean, std, functional_form) = a.priors.get(name)
		val = getattr(a, name)
		prior_names.append(name)
		if functional_form == 0:
			prior.append(gaussian(val, mean, std))
		elif functional_form == 1:
			prior.append(lognorm(val, mean, std))
	a.prior = prior

	# check the borders of the free parameters
	for name in a.to_optimize:
		(lower, upper) = a.constraints.get(name)
		val = getattr(a, name)
		if lower is not None and val<lower:
			assert False, '%s lower border is violated' %(name)
		if upper is not None and val>upper:
			assert False, '%s upper border is violated' %(name)

	if not a.testing_output:
		print(changing_parameter,mp.current_process()._identity[0])#,a.observational_constraints_index
	else:
		print(changing_parameter)
	
	### So that the parameter can be plotted in linear space
	if 'log10_N_0' in a.to_optimize:
		a.N_0 = np.power(10,a.log10_N_0)
	if 'log10_sn1a_time_delay' in a.to_optimize:
		a.sn1a_time_delay = np.power(10,a.log10_sn1a_time_delay)
	if 'log10_starformation_efficiency' in a.to_optimize:
		a.starformation_efficiency = np.power(10,a.log10_starformation_efficiency)
	if 'log10_gas_reservoir_mass_factor' in a.to_optimize:
		a.gas_reservoir_mass_factor = np.power(10,a.log10_gas_reservoir_mass_factor)
	if 'log10_sfr_scale' in a.to_optimize:
		a.sfr_scale = np.power(10,a.log10_sfr_scale)

	if a.imf_type_name == 'salpeter':
		a.imf_parameter = (a.high_mass_slope)
	elif a.imf_type_name == 'Chabrier_2':
		a.imf_parameter = (a.chabrier_para1, a.chabrier_para2, a.chabrier_para3,a.high_mass_slope)
	elif a.imf_type_name == 'normed_3slope':	
		a.imf_parameter = (a.imf_slope_1,a.imf_slope_2,a.high_mass_slope,a.imf_break_1,a.imf_break_2)
	if a.time_delay_functional_form == 'maoz':
		a.sn1a_parameter = [a.N_0,a.sn1a_time_delay,a.sn1a_exponent,a.dummy]
	elif a.time_delay_functional_form == 'normal':
		a.sn1a_parameter = [a.number_of_pn_exlopding,a.sn1a_time_delay,a.sn1a_timescale,a.sn1a_gauss_beginning]
	elif a.time_delay_functional_form == 'gamma_function':
		a.sn1a_parameter = [a.sn1a_norm,a.sn1a_a_parameter,a.sn1a_beginning,a.sn1a_scale]
	return(a)

def posterior_function_real(changing_parameter,a):
	'''
	This is the actual posterior function. But the functionality is explained in posterior_function.
	'''
	
	start_time = time.time()
	# the values in a are updated according to changing_parameters and the prior list is appended
	a = extract_parameters_and_priors(changing_parameter, a)
	

	# the log prior is calculated
	prior = sum(np.log(a.prior))

	
	precalculation = time.time()
	#print('precalculation: ', start_time - precalculation)

	## The endtime is changed for the actual calculation but restored to default afterwards
	backup = a.end ,a.time_steps, a.total_mass
	
	# call Chempy and return the abundances at the end of the simulation = time of star's birth and the corresponding element names as a list
	abundance_list,elements_to_trace = cem_real2(a)
	a.end ,a.time_steps, a.total_mass = backup
	
	# The last two entries of the abundance list are the Corona metallicity and the SN-ratio
	abundance_list = abundance_list[:-2]
	elements_to_trace = elements_to_trace[:-2]

	model = time.time()
	#print('model: ', precalculation - model)

	# a likelihood is calculated where the model error is optimized analytically
	likelihood, element_list, model_error, star_error_list, abundance_list, star_abundance_list = likelihood_function(a.stellar_identifier, abundance_list, elements_to_trace)


	error_optimization = time.time()
	#print('error optimization: ', model - error_optimization)

	if not a.testing_output:
		print('prior = ', prior, 'likelihood = ', likelihood, mp.current_process()._identity[0])
	else:
		print('prior = ', prior, 'likelihood = ', likelihood)

	return(prior+likelihood,[0])


def posterior_function_for_minimization(changing_parameter,a):
	'''
	calls the posterior function but just returns the negative log posterior instead of posterior and blobs
	'''
	posterior, blobs = posterior_function(changing_parameter,a)
	return -posterior