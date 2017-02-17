import numpy as np
import os
from .sfr import SFR
from .infall import INFALL, PRIMORDIAL_INFALL
from .time_integration import ABUNDANCE_MATRIX
from .solar_abundance import solar_abundances
from .making_abundances import mass_fraction_to_abundances
from numpy.lib.recfunctions import append_fields
import time
from .data_to_test import elements_plot, arcturus, sol_norm, plot_processes, save_abundances,  cosmic_abundance_standard, ratio_function, star_function, gas_reservoir_metallicity
import multiprocessing as mp
from .wrapper import SSP_wrap, initialise_stuff, Chempy

def gaussian_log(x,x0,xsig):
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)

def lognorm_log(x,mu,factor):
	'''
	this function provides Prior probability distribution where the factor away from the mean behaves like the sigma deviation in normal_log 
	for example if mu = 1 and factor = 2 for
	1 it returns 0
	for 0,5 and 2 it returns -0.5
	for 0.25 and 4 it returns -2.0
	and so forth
	Can be used to specify the prior on the yield factors
	'''
	y = np.log(np.divide(x,mu))
	y = np.divide(y,np.log(factor))
	y = gaussian_log(y,0.,1.)
	#y = np.log(y)
	return y

def cem(changing_parameter,a):
	try:
		posterior, blobs = cem_real(changing_parameter,a)
		return posterior, blobs
	except Exception as ex:
		import traceback; traceback.print_exc()
		return -np.inf, [0]

def cem_real(changing_parameter,a):
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

	if a.imf_type_name == 'salpeter':
		a.imf_parameter = (a.imf_slope)
	elif a.imf_type_name == 'Chabrier_2':
		a.imf_parameter = (a.chabrier_para1, a.chabrier_para2, a.chabrier_para3,a.high_mass_slope)
	elif a.imf_type_name == 'normed_3slope':	
		a.imf_parameter = (a.imf_slope_1,a.imf_slope_2,a.imf_slope_3,a.imf_break_1,a.imf_break_2)
	if a.time_delay_functional_form == 'maoz':
		a.sn1a_parameter = [a.N_0,a.sn1a_time_delay,a.sn1a_exponent,a.dummy]
	elif a.time_delay_functional_form == 'normal':
		a.sn1a_parameter = [a.number_of_pn_exlopding,a.sn1a_time_delay,a.sn1a_timescale,a.sn1a_gauss_beginning]
	elif a.time_delay_functional_form == 'gamma_function':
		a.sn1a_parameter = [a.sn1a_norm,a.sn1a_a_parameter,a.sn1a_beginning,a.sn1a_scale]

	basic_solar = solar_abundances()
	getattr(basic_solar, a.solar_abundance_name)()
	elements_to_trace = a.elements_to_trace

	### Model is calculated
	if a.calculate_model:
		cube, abundances= Chempy(a)
		cube1 = cube.cube
		gas_reservoir = cube.gas_reservoir
		directory = 'model_temp/'
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
		np.save("blobs_name_list", a.names)
	if np.isnan(sum(a.probability)):
		return -np.inf, [0]
	if a.testing_output:
		print('l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability))#, mp.current_process()._identity[0]
	else:
		print('l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability),'|', mp.current_process()._identity[0])
	return sum(a.probability) + sum(a.prior), np.hstack((a.probability,a.abundance_list,a.prior))
