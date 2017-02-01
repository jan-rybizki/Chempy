import numpy as np
from sfr import SFR
from infall import INFALL, PRIMORDIAL_INFALL
from time_integration import ABUNDANCE_MATRIX
from solar_abundance import solar_abundances
from making_abundances import mass_fraction_to_abundances
from numpy.lib.recfunctions import append_fields
import time
from data_to_test import elements_plot, arcturus, sol_norm, plot_processes, save_abundances,  cosmic_abundance_standard, ratio_function, star_function, gas_reservoir_metallicity
import multiprocessing as mp
from wrapper import SSP_wrap

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
	except Exception, ex:
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
			print '%s lower border is violated' %(name)
			return -np.inf, [0]
		if upper is not None and val>upper:
			print '%s upper border is violated' %(name)
			return -np.inf, [0]

	if not a.testing_output:
		print changing_parameter,mp.current_process()._identity[0]#,a.observational_constraints_index
	else:
		print changing_parameter
	
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

	if a.calculate_model:

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

		elements_to_trace = a.elements_to_trace

		basic_pri_gas_at_start = PRIMORDIAL_INFALL(list(elements_to_trace),np.copy(basic_solar.table))
		basic_pri_gas_at_start.primordial(0)
		
		basic_primordial = PRIMORDIAL_INFALL(list(elements_to_trace),np.copy(basic_solar.table))
		basic_primordial.primordial(0)

		cube = ABUNDANCE_MATRIX(np.copy(basic_sfr.t),np.copy(basic_sfr.sfr),np.copy(basic_infall.infall),list(elements_to_trace),list(basic_primordial.symbols),list(basic_primordial.fractions),float(a.gas_at_start),list(basic_pri_gas_at_start.symbols),list(basic_pri_gas_at_start.fractions),float(a.gas_reservoir_mass_factor),float(a.outflow_feedback_fraction),a.check_processes,a.starformation_efficiency,a.gas_power, a.sfr_factor_for_cosmic_accretion, a.cosmic_accretion_elements,a.cosmic_accretion_element_fractions)
		basic_ssp = SSP_wrap(a)
		time_preload = (time.time() - start_time)
		### Here the real chemical evolution is started and iterated over the SFH
		for i in range(len(basic_sfr.t)-1):
			j = len(basic_sfr.t)-i
			element_fractions = []
			for item in elements_to_trace:
				element_fractions.append(np.copy(cube.cube[item][max(i-1,0)]/cube.cube['gas'][max(i-1,0)]))## gas element fractions from one time step before	
			metallicity = float(cube.cube['Z'][i])
			time_steps = np.copy(basic_sfr.t[:j])
			basic_ssp.calculate_feedback(float(metallicity), list(elements_to_trace), list(element_fractions), np.copy(time_steps))
			cube.advance_one_step(i+1,np.copy(basic_ssp.table),np.copy(basic_ssp.sn2_table),np.copy(basic_ssp.agb_table),np.copy(basic_ssp.sn1a_table))
			if cube.cube['gas'][i] < 0:
				print i, basic_sfr.t[i]
				print 'gas became negative. returning -inf'
				return -np.inf, [0]
			if cube.gas_reservoir['gas'][i] < 0:
				print 'gas_reservoir became negative. returning -inf'
				return -np.inf, [0]
		## Model is safed temporarily
		time_calculation = (time.time() - start_time - time_preload)
		abundances,elements,numbers = mass_fraction_to_abundances(np.copy(cube.cube),np.copy(basic_solar.table))
		weights = cube.cube['sfr']
		abundances = append_fields(abundances,'weights',weights)
		abundances = np.array(abundances)
		np.save('model_temp/%s_elements_to_trace' %(a.name_string), elements_to_trace)
		np.save('model_temp/%s_gas_reservoir' %(a.name_string),cube.gas_reservoir)
		np.save('model_temp/%s_cube' %(a.name_string),cube.cube)
		np.save('model_temp/%s_abundances' %(a.name_string),abundances)
		cube1 = cube.cube
		gas_reservoir = cube.gas_reservoir
	else:
		cube1 = np.load('model_temp/%s_cube.npy' %(a.name_string))
		abundances = np.load('model_temp/%s_abundances.npy' %(a.name_string))
		gas_reservoir = np.load('model_temp/%s_gas_reservoir.npy' %(a.name_string))
		elements_to_trace = np.load('model_temp/%s_elements_to_trace.npy' %(a.name_string))


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
		print a.names
		np.save("name_list", a.names)
	if np.isnan(sum(a.probability)):
		return -np.inf, [0]
	if a.testing_output:
		print 'l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability)#, mp.current_process()._identity[0]
	else:
		print 'l: ', sum(a.probability), 'pr: ', sum(a.prior), 'po: ', sum(a.prior) + sum(a.probability),'|', mp.current_process()._identity[0]
	return sum(a.probability) + sum(a.prior), np.hstack((a.probability,a.abundance_list,a.prior))