import numpy as np


class ModelParameters(object):
	'''
	In this class the model parameters are specified. It contains a lot of information which is (not always) necessary to run Chempy.
	The individual definitions are given as comments.
	'''
	
	# Which zero point of abundances shall be used. Asplund 2005 is corrected to VESTA abundances
	solar_abundance_name_list = ['Lodders09','Asplund09','Asplund05_pure_solar','Asplund05_apogee_correction']
	solar_abundance_name_index = 1
	solar_abundance_name = solar_abundance_name_list[solar_abundance_name_index]

	# Observational constraints
	#stellar_identifier_list = ['Proto-sun', 'Arcturus', 'B-stars']
	#stellar_identifier_list = ['2M01233744+3414451', '2M02484368+3106550', '2M05510326+1129561', '2M09031459+0648573', '2M09422500+4846338', '2M02011031+2426397', '2M09055837+0505324', '2M20092234+5601366']
	#indices = [78,130,122,156,113,34, 128,167] # low alpha sequence
	#indices = [0, 163, 27,  98,  95,  17,  71,  79] # random
	#indices = [158, 24, 152, 56, 100, 21, 17, 126] # This is the list for middle alpha sequence
	#indices = [147, 0, 3, 128, 1, 156, 113, 110] # extremes in alpha over iron space
	#indices = [85, 94, 15, 110, 30, 11, 7, 3] # high alpha sequence
	#indices = [78,130,122,156,113,34, 128,167,85, 94, 15, 110, 30, 11, 7, 3] #low alpha + high alpha
	#indices = [0, 163, 27,  98,  95,  17,  71,  79, 78,130,122,156,113,34, 128,167,85, 94, 15, 110, 30, 11, 7, 3] #low alpha + high alpha + random
	#stellar_identifier_list = []
	#for item in indices:
	#	stellar_identifier_list.append("Rob_%d" %item)
	#stellar_identifier_list = ['Proto-sun', 'Arcturus', 'B-stars']
	# 'prior' can be used as stellar_identifier, then the prior will be sampled with Chempy.wrapper.mcmc() routine
	stellar_identifier_list = ['Proto-sun']
	stellar_identifier = 'Proto-sun'

	# Convergense parameters of minimization and MCMC
	maxiter_minimization = 500
	min_mcmc_iterations = 300
	mcmc_tolerance = 0.5
	gibbs_sampler_tolerance = 1e-1
	gibbs_sampler_maxiter = 10
	tol_minimization = 1e-1
	nwalkers = 64
	mburn = 1
	save_state_every = 1
	m = 1000 # For 7 free parameters 300 iterations are usually enough. The mcmc routine is stopping after 300 if the posterior mean is converged for more than 200 iterations.
	error_marginalization = True # Marginalizing over the model error or using the best model error value
	flat_model_error_prior = [0.,1.,51] # Flat prior for the error marginalization [begin, end, number of evaluations inbetween]
	beta_error_distribution = [True, 1, 3] # Instead of a flat prior for the error marginalization we use a beta distribution with a = 1 and b = 3 as default (wikipedia and scipy have the same parametrization) putting more weight to small model errors
	zero_model_error = False # a boolean that can be used to restore the old Chempy behaviour of 0 model error, will only work if error_marginalization is set to False
	send_email = False
	
	verbose = 0
	# Time discretization, so far only linear time-steps are implemented
	start = 0 # birth of disc, always set to 0
	end = 13.5
	time_steps = 28#541#241#35#1401
	total_mass = 1#45.07
	stochastic_IMF = False
	number_of_models_overplotted = 1 ### with the positions from an mcmc run
	testing_output = False
	summary_pdf = False
	name_string = 'Chempy_default'
	parameter_names = [r'$\alpha_\mathrm{IMF}$',r'$\log_{10}\left(\mathrm{N}_\mathrm{Ia}\right)$',r'$\log_{10}\left(\tau_\mathrm{Ia}\right)$',r'$\log_{10}\left(\mathrm{SFE}\right)$',r'$\log_{10}\left(\mathrm{SFR}_\mathrm{peak}\right)$',r'$\mathrm{x}_\mathrm{out}$']
	# SFR still model A from Just&Jahreiss 2010 should be changed
	# arbitrary function can be implemented here
	basic_sfr_name_list = ['model_A','gamma_function','prescribed', 'doubly_peaked']
	basic_sfr_index = 1
	basic_sfr_name = basic_sfr_name_list[basic_sfr_index]
	if basic_sfr_name == 'model_A':
		mass_factor = 1.
		S_0 = 45.07488
		t_0 = 5.6
		t_1 = 8.2
	elif basic_sfr_name == 'gamma_function':
		mass_factor = 1.
		S_0 = 1#45.07488
		a_parameter = 2
		sfr_beginning = 0
		sfr_scale = 3.5 # SFR peak in Gyr for a = 2
	elif basic_sfr_name == 'prescribed':
		mass_factor = 1.
		name_of_file = 'input/Daniel_Weisz/ic1613.lcid.final.sfh'
	elif basic_sfr_name == 'doubly_peaked':
		mass_factor = 1.
		S_0 = 45.07488
		peak_ratio = 0.8
		sfr_decay = 3.5
		sfr_t0 = 2.
		peak1t0 = 0.8
		peak1sigma = 0.8

	basic_infall_name_list = ["exponential","constant","sfr_related","peaked_sfr","gamma_function"]
	basic_infall_index = 2
	basic_infall_name = basic_infall_name_list[basic_infall_index]
	starformation_efficiency = 0.
	gas_power = 0.
	if basic_infall_name == 'sfr_related':
		starformation_efficiency = np.power(10,-0.3)
		gas_power = 1.0
	if basic_infall_name == 'exponential':
		infall_amplitude = 10 # not needed just a dummy
		tau_infall = -0.15
		infall_time_offset = 0
		c_infall = -1.
		norm_infall = 0.9
	if basic_infall_name == 'gamma_function':
		norm_infall = 1.0 # not needed just a dummy
		infall_a_parameter = 2
		infall_beginning = 0
		infall_scale = 3.3

	yield_table_name_sn2_list = ['chieffi04','Nugrid','Nomoto2013','Portinari', 'chieffi04_net', 'Nomoto2013_net']
	yield_table_name_sn2_index = 2
	yield_table_name_sn2 = yield_table_name_sn2_list[yield_table_name_sn2_index]

	yield_table_name_hn_list = ['Nomoto2013']
	yield_table_name_hn_index = 0
	yield_table_name_hn = yield_table_name_hn_list[yield_table_name_hn_index]

	##### Karakas2016 needs much more calculational resources (order of magnitude) using 2010 net yields from Karakas are faster and only N is significantly underproduced
	yield_table_name_agb_list = ['Karakas','Nugrid','Karakas_net_yield','Ventura','Karakas16_net']
	yield_table_name_agb_index = 2
	yield_table_name_agb = yield_table_name_agb_list[yield_table_name_agb_index]

	yield_table_name_1a_list = ['Iwamoto','Thielemann','Seitenzahl']
	yield_table_name_1a_index = 2
	yield_table_name_1a = yield_table_name_1a_list[yield_table_name_1a_index]

	mmin = 0.1
	mmax = 100
	mass_steps = 5000 #2000 # 200000

	imf_type_name_list = ['normed_3slope','Chabrier_1','Chabrier_2','salpeter','BrokenPowerLaw']
	imf_type_index = 1
	imf_type_name = imf_type_name_list[imf_type_index]
	if imf_type_name == 'Chabrier_2':
		chabrier_para1 = 22.8978
		chabrier_para2 = 716.4
		chabrier_para3 = 0.25
		high_mass_slope = -2.3
		imf_parameter = (22.8978, 716.4, 0.25,-2.29)
	if imf_type_name == 'Chabrier_1':
		chabrier_para1 = 0.69
		chabrier_para2 = 0.079
		high_mass_slope = -2.29
		imf_parameter = (0.69, 0.079, -2.29)
	if imf_type_name == 'salpeter':
		imf_slope = 2.35
		imf_parameter = (2.35)
	if imf_type_name == 'BrokenPowerLaw':
		imf_break_1 = 0.5
		imf_break_2 = 1.39
		imf_break_3 = 6
		imf_slope_1 = -1.26
		imf_slope_2 = -1.49
		imf_slope_3 = -3.02
		imf_slope_4 = -2.3
		imf_parameter = ((0.5,1.39,6),(-1.26,-1.49,-3.02,-2.3))
	if imf_type_name == 'normed_3slope':	
		imf_break_1 = 0.5
		imf_break_2 = 1.0
		imf_slope_1 = -1.3
		imf_slope_2 = -2.3
		imf_slope_3 = -2.29
		imf_parameter = (imf_slope_1,imf_slope_2,imf_slope_3,imf_break_1,imf_break_2)
	name_infall_list = ['primordial','solar','simple','alpha']
	name_infall_index = 1
	name_infall = name_infall_list[name_infall_index]

	interpolation_list = ['linear','logarithmic']
	interpolation_index = 1
	interpolation_scheme = interpolation_list[interpolation_index] ## could be a variant to change the interpolation scheme
	stellar_lifetimes_list = ['Argast_2000','Raiteri_1996']
	stellar_lifetimes_index = 0
	stellar_lifetimes = stellar_lifetimes_list[stellar_lifetimes_index] ## which stellar lifetime approximation to use

	sn2_to_hn = 1.

	sn2mmin = 8.
	sn2mmax = 100.

	bhmmin = float(sn2mmax) ## maximum of hypernova
	bhmmax = float(mmax) ## maximum of the IMF
	percentage_of_bh_mass = 0.25 # the rest 75% will be given back to the ISM with the abundances from the step before

	agbmmin = 0.5
	agbmmax = 8

	sagbmmin = float(agbmmax)
	sagbmmax = float(sn2mmin)
	percentage_to_remnant = 0.13 # see Kobayashi 2011 the remnant mass is about 13%

	time_delay_functional_form_list = ['normal','maoz','gamma_function']
	time_delay_index = 1
	time_delay_functional_form = time_delay_functional_form_list[time_delay_index]
	if time_delay_functional_form == 'maoz':
		N_0 = np.power(10,-2.75)
		sn1a_time_delay = np.power(10,-0.8)
		sn1a_exponent = 1.12
		dummy = 0.0
		sn1a_parameter = [N_0,sn1a_time_delay,sn1a_exponent,dummy]
	if time_delay_functional_form == 'normal':
		number_of_pn_exlopding = 0.003
		sn1a_time_delay = 1.
		sn1a_timescale = 3.2
		sn1a_gauss_beginning = 0.25
		sn1a_parameter = [number_of_pn_exlopding,sn1a_time_delay,sn1a_timescale,sn1a_gauss_beginning]
	if time_delay_functional_form == 'gamma_function':
		sn1a_norm = 0.0024 #number of sn1a exploding within end of simulation time per 1Msun
		sn1a_a_parameter = 1.3
		sn1a_beginning = 0
		sn1a_scale = 3
		sn1a_parameter = [sn1a_norm,sn1a_a_parameter,sn1a_beginning,sn1a_scale]
	sn1ammin = 1#float(agbmmin) #Maoz Timedelay should be independent of sn1a_mmin and sn1a_mmax. N_0 just determines the number of SN1a exploding per 1Msun over the time of 15Gyr
	sn1ammax = 8#float(sagbmmax)
	gas_at_start = 0. #*dt yields the Msun/pc^2 value

	gas_reservoir_mass_factor = np.power(10,0.0)#3.0
	sfr_factor_for_cosmic_accretion = 1.
	cosmic_accretion_elements = ['H','He']
	cosmic_accretion_element_fractions = [0.76,0.24]
	outflow_feedback_fraction = 0.5
	## various output modes
	check_processes = False
	only_net_yields_in_process_tables = True
	calculate_model = True #just loading the outcome of the last ssp if False


	####### Evaluate model
	element_names = ['He','C', 'N', 'O', 'F','Ne','Na', 'Mg', 'Al', 'Si', 'P','S', 'Ar','K', 'Ca','Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni']#, 'Zn','Y', 'Ba']# Runs with sun
	elements_to_trace = ['Al', 'Ar', 'B', 'Be', 'C', 'Ca', 'Cl', 'Co', 'Cr', 'Cu', 'F', 'Fe', 'Ga', 'Ge', 'H', 'He', 'K', 'Li', 'Mg', 'Mn', 'N', 'Na', 'Ne', 'Ni', 'O', 'P', 'S', 'Sc', 'Si', 'Ti', 'V', 'Zn']

	#observational_constraints_index = ['sol_norm']#['gas_reservoir','sn_ratio','sol_norm']#,'wildcard ','cas','arcturus','stars_at_end', 'plot_processes', 'save_abundances', 'elements']
	arcturus_age = 7.1# 7.1 +1.5 -1.2

	produce_mock_data = False
	use_mock_data = False
	error_inflation = 1.


	# If some parameter is in to optimise there needs to be a prior and constraints defined
	if True:
		#prior
		SSP_parameters =  [-2.29 ,-2.75 ,	-0.8 ]#,0.2]#, 0.7, 0.3, 0.0]
		SSP_parameters_to_optimize = ['high_mass_slope', 'log10_N_0', 'log10_sn1a_time_delay']#,'log10_sfr_factor_for_cosmic_accretion']#,'log10_gas_reservoir_mass_factor','log10_a_parameter','log10_gas_power']
	else:
		SSP_parameters = []
		SSP_parameters_to_optimize = []
	assert len(SSP_parameters) == len(SSP_parameters_to_optimize)
	if True:
		#prior
		ISM_parameters =  [-0.3, 0.55,	0.5]#, 0.3]#,0.2]#, 0.7, 0.3, 0.0]
		ISM_parameters_to_optimize = ['log10_starformation_efficiency', 'log10_sfr_scale', 'outflow_feedback_fraction']#,'log10_gas_reservoir_mass_factor']#,'log10_sfr_factor_for_cosmic_accretion']#,'log10_gas_reservoir_mass_factor','log10_a_parameter','log10_gas_power']
	else:
		ISM_parameters = []
		ISM_parameters_to_optimize = []
	assert len(ISM_parameters) == len(ISM_parameters_to_optimize)

	p0 = np.hstack((SSP_parameters,ISM_parameters))
	to_optimize = np.array(SSP_parameters_to_optimize + ISM_parameters_to_optimize)
	ndim = len(to_optimize)

	constraints = {
	'high_mass_slope' : (-4.,-1.),
	'log10_N_0' : (-5,-1), 
	'log10_sn1a_time_delay' : (-3,1.),
	'log10_starformation_efficiency' : (-3,2),
	'log10_sfr_scale' : (-1,1),
	'sfr_scale' : (0.0,None),
	'outflow_feedback_fraction' : (0.,1.),
	'log10_gas_reservoir_mass_factor': (None,None),
	

	'N_0' : (0.,1.),
	'sn1a_time_delay' : (0.,end),
	'a_parameter' : (0.,None),
	'starformation_efficiency' : (0.,None), 
	'gas_power': (1.,2.),
	'log10_a_parameter' : (None,None),
	'log10_gas_power' : (None,None),
	'log10_gas_reservoir_mass_factor': (None,None),
	'log10_sfr_factor_for_cosmic_accretion': (None,None),

	'mass_factor' : (0,None),
	'norm_infall' : (0.,2.),
	'tau_infall' : (None,None),
	'c_infall' : (None,None),
	'gas_at_start' : (0.,2.),
	'gas_reservoir_mass_factor' : (0.,20.),
	'infall_scale' : (0.0,end),
	'sn1a_norm' : (0.,None),
	'sn1a_scale' : (0.,None),
	}
	# the prior entry is (mean,std,0)
	# functional form 0 is a gaussian with log values and 1 is for fractions where the sigma distances are in factors from the mean (see cem_function.py)
	# for functional form 1 read (mean,factor,1)
	priors = {
	## gaussian priors
	'high_mass_slope' : (-2.3,0.3,0),	
	'log10_N_0' : (-2.75,0.3,0),
	'log10_sn1a_time_delay' : (-0.8,0.3,0),
	'log10_starformation_efficiency' : (-0.3,0.3,0),
	'log10_sfr_scale' : (0.55,0.1,0),
	'sfr_scale' : (3.5,1.5,0),
	'outflow_feedback_fraction' : (0.5,0.1,0),
	'log10_gas_reservoir_mass_factor' : (0.3,0.3,0),

	'a_parameter' : (3.,3.,0),
	'infall_scale' : (3.3,0.5,0),
	'gas_power': (1.5,0.2,0),
	'log10_sfr_factor_for_cosmic_accretion': (0.2,0.3,0),
	'log10_a_parameter' : (0.3,0.2,0),
	'log10_gas_power' : (0,0.15,0),
	
	## Priors on factors
	'starformation_efficiency' : (0.5,3.,1), 
	'mass_factor' : (1.,1.2,1),
	'norm_infall' : (1.,1.2,1),
	'sn1a_time_delay' : (0.3,3.,1),
	'N_0' : (0.001,3.,1),
	'gas_at_start' : (0.1,2.,1),
	'gas_reservoir_mass_factor' : (3.,2.,1),
	}
