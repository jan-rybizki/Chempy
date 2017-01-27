import numpy as np


def gaussian_log(x,x0,xsig):
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)


class ModelParameters(object):
	# Which zero point of abundances shall be used. Asplund 2005 is corrected to VESTA abundances
	solar_abundance_name_list = ['Lodders09','Asplund09','Asplund05_pure_solar','Asplund05_apogee_correction']
	solar_abundance_name_index = 1
	solar_abundance_name = solar_abundance_name_list[solar_abundance_name_index]

	# Time discretization, so far only linear time-steps are implemented
	start = 0 # birth of disc, always set to 0
	end = 13.5
	time_steps = 28#541#241#35#1401
	total_mass = 1#45.07
	stochastic_IMF = False
	number_of_models_overplotted = 1 ### with the positions from an mcmc run
	testing_output = False

	only_one_SSP = False
	if only_one_SSP:
		metallicity_of_SSP = 0.0134
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

	basic_infall_name_list = ["exponential","constant","sfr_related","peaked_sfr","just","simon",'sarah','gamma_function']
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

	yield_table_name_sn2_list = ['chieffi04','Nugrid','Nomoto2013','Portinari']
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

	mmin = 0.08
	mmax = 100
	mass_steps = 5000 #2000 # 200000

	imf_type_name_list = ['normed_3slope','Chabrier_1','Chabrier_2','salpeter','BrokenPowerLaw']
	imf_type_index = 2
	imf_type_name = imf_type_name_list[imf_type_index]
	if imf_type_name == 'Chabrier_2':
		chabrier_para1 = 22.8978
		chabrier_para2 = 716.4
		chabrier_para3 = 0.25
		high_mass_slope = -2.3

	if imf_type_name == 'Chabrier_1':
		chabrier_para1 = 0.852464
		chabrier_para2 = 0.237912
		chabrier_para3 = 0.69
		chabrier_para4 = 0.079

	if imf_type_name == 'salpeter':
		imf_slope = 2.35

	if imf_type_name == 'BrokenPowerLaw':
		imf_break_1 = 0.5
		imf_break_2 = 1.39
		imf_break_3 = 6
		imf_slope_1 = -1.26
		imf_slope_2 = -1.49
		imf_slope_3 = -3.02
		imf_slope_4 = -2.3

	if imf_type_name == 'normed_3slope':	
		imf_break_1 = 0.5
		imf_break_2 = 1.0
		imf_slope_1 = -1.3
		imf_slope_2 = -2.2
		imf_slope_3 = -2.7

	name_infall_list = ['primordial','solar','simple','alpha']
	name_infall_index = 1
	name_infall = name_infall_list[name_infall_index]

	if name_infall == 'solar':
		begin_metallicity_dex = -0.35

	interpolation_list = ['linear','logarithmic']
	interpolation_index = 1
	interpolation_scheme = interpolation_list[interpolation_index] ## could be a variant to change the interpolation scheme
	stellar_lifetimes_list = ['Argast_2000','Raiteri_1996']
	stellar_lifetimes_index = 0
	stellar_lifetimes = stellar_lifetimes_list[stellar_lifetimes_index] ## which stellar lifetime approximation to use
	with_extra_output = False # For the example SSP to get

	sn2_to_hn = 0.5

	sn2_method_list = ['IRA', 'discretised']
	sn2_method_index = 1
	sn2_method_name = sn2_method_list[sn2_method_index]
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
	if time_delay_functional_form == 'normal':
		number_of_pn_exlopding = 0.003
		sn1a_time_delay = 1.
		sn1a_timescale = 3.2
		sn1a_gauss_beginning = 0.25
	if time_delay_functional_form == 'gamma_function':
		sn1a_norm = 0.0024 #number of sn1a exploding within end of simulation time per 1Msun
		sn1a_a_parameter = 1.3
		sn1a_beginning = 0
		sn1a_scale = 3

	sn1ammin = 1#float(agbmmin) #Maoz Timedelay should be independent of sn1a_mmin and sn1a_mmax. N_0 just determines the number of SN1a exploding per 1Msun over the time of 15Gyr
	sn1ammax = 8#float(sagbmmax)
	#gas_at_start_parameters = (-0.7,0.4,0,-0.03)
	iron_at_start_in_dex = -1.
	gas_at_start = 0. #/40 yields the Msun/pc^2 value

	gas_at_start_name_list = ['sn2','primordial']
	gas_at_start_index = 1
	gas_at_start_name = gas_at_start_name_list[gas_at_start_index]
	gas_reservoir_mass_factor = np.power(10,0.3)#3.0
	sfr_factor_for_cosmic_accretion = 1.
	cosmic_accretion_elements = ['H','He']
	cosmic_accretion_element_fractions = [0.76,0.24]
	outflow_feedback_fraction = 0.5
	## various output modes
	check_processes = True
	only_net_yields_in_process_tables = True
	calculate_model = True #just loading the outcome of the last ssp if False
	Verbose = False
	#sigma_astro = 0.03 # astrophysical scatter added to the reported observational error

	Q = 2.
	lnodds = -3.
	astronomical_spread = 0.05
	normalising_element = 'Fe' # to which elements the others are normalised to in the likelihood calculation
	dr13 = True #switch between Dr13 and Dr12 cannon data
	survey_name = 'APOGEE-Keith' # chose from RAVE4, APOGEE-DR12, APOGEE-DR13, APOGEE-Keith, GES-DR4, Cepheids, Bensby14
	## Rave giants: 4250,5250,1.7,2.8, Rave dwarfs: 5250,7000,3.8,5.0, APOGEE RC: 4500,5200,2.,3.1
	dist = [0.,10.] #in kpc
	rgal = [7.5,8.5]
	#rgal = [6.8,7.8]
	#rgal = [8.8,9.8]

	zgal = [-0.4,0.4]
	form_factor = 'red-clump'

	if form_factor == 'red-clump':
		logg = [2.,3.1]
		teff = [4500.,5200.]
	elif form_factor == 'rave-red-clump':
		logg = [1.7,2.8]
		teff = [4250.,5250.]
	elif form_factor == 'rave-dwarfs':
		logg = [3.8,5.0]
		teff = [5250.,7000.]
	elif form_factor == 'all':
		logg = [0,6.0]
		teff = [2250.,16000.]
	elif form_factor == 'uves-dwarfs':
		logg = [3.5,4.5]
		teff = [5500.,6500.]
	#elif form_factor == 'uves-giants':
	#	logg = [3.5,4.5]
	#	teff = [5500.,6500.]
	many_survey_plot = False
	plot_model = True

	name_string = 'single_suncasarc'#%d' %(time_steps)#'mass_analytic'#%d' %(total_mass/1e2)

	nwalkers = 64
	mburn = 1
	save_state_every = 1
	m = 1000
	####### Evaluate model
	#V,K,
	# element lists
	model_apogee_all_dr12 =  ['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni']
	model_apogee_all_dr13 =  ['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe','Ni','P','Cr','Co','Cu','Rb']
	apogee_keith_elements = ['Fe','C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Ni','P','Cr','Co','Cu','Ba'] 
	bensby_all_elements = ['O','Na','Mg','Al','Si','Ca','Ti','Cr','Fe','Ni','Zn','Y','Ba']
	element_list_arcturus = ['Fe', 'C', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Zn']
	element_list_cas = ['H','He','O','Mg','Si','Fe','Ne','C','N']
	RAVE4_elements = ['Al','Si','Fe','Ti','Ni','Mg']
	Cepheids_elements = ['Fe','Y','La','Ce','Nd','Eu','Na','Al','Mg','Si','Ca'] 
	Ges4_elements = ['Fe','Mg','O']

	element_names = ['O','Na','Mg','Al','Si','Ca','Ti','Cr','Fe','Ni','Zn','Y','Ba']## Bensby
	element_names = ['O','Mg','Al','Si','P','Ca','Cr','Mn','Fe','Ni']## DR13
	element_names = ['Fe', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Ti', 'Cr', 'Mn', 'Co', 'Ni', 'Zn']## Runs with arcturus
	element_names = ['He','C', 'N', 'O', 'F','Ne','Na', 'Mg', 'Al', 'Si', 'P','S', 'Ar','K', 'Ca','Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni']#, 'Zn','Y', 'Ba']# Runs with sun
	
	element_names_for_triangle_plot = ['Fe','Mg','Mn','C','Ca','Co','Na']
	#'apogee','gas_reservoir','gas_at_end','stars_at_end','sn_ratio','cas','sol_norm','sfr','mock_abundances','plot_processes','alpha_corner','all_alpha','save_abundances','elements'
	observational_constraints_index = ['apogee','gas_reservoir','stars_at_end','sn_ratio','gas_at_end','cas','sol_norm','arcturus']#'flexible_survey',
	observational_constraints_index = []#['plot_processes','sn_ratio','gas_reservoir','stars_at_end','cas','arcturus','sol_norm','elements']#['cas', 'sol_norm', 'arcturus','sn_ratio','gas_reservoir','stars_at_end']
	observational_constraints_index = ['gas_reservoir','sn_ratio','sol_norm','stars_at_end','cas','arcturus','elements', 'plot_processes']#'gas_at_end','sfr','stars_at_end'
	observational_constraints_index = ['gas_reservoir','stars_at_end','sol_norm','arcturus','elements', 'plot_processes','cas']#'gas_at_end','sfr','stars_at_end'
	#observational_constraints_index = ['sol_norm', 'sn_ratio', 'gas_reservoir','arcturus','cas','stars_at_end','plot_processes']
	observational_constraints_index = ['gas_reservoir','sn_ratio','cas','sol_norm','arcturus','stars_at_end']
	arcturus_age = 7.1# 7.1 +1.5 -1.2
	make_multi_zone = False
	multi_zone_list = [['sol_norm','gas_reservoir','sn_ratio'], ['arcturus','gas_reservoir','sn_ratio']]

	produce_mock_data = False
	use_mock_data = False
	error_inflation = 1.
	#observational_constraints_index = ['apogee','save_abundances','plot_processes']
	# If some parameter is in to optimise there needs to be a prior and constraints defined
	# Some names should not be changed or used for optimizing names. 'offset' and 'sigma' are such names
	if False:
		#error_parameters = [-2., 2., 0.05]#, 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. ]
		#error_parameters_to_optimize = ['lnodds','astronomical_spread','Q']#,,'sigma_O','sigma_Na','sigma_Mg','sigma_Al','sigma_Si','sigma_Ca','sigma_Ti','sigma_Cr','sigma_Fe','sigma_Ni','sigma_Zn','sigma_Y','sigma_Ba','Q_O','Q_Na','Q_Mg','Q_Al','Q_Si','Q_Ca','Q_Ti','Q_Cr','Q_Fe','Q_Ni','Q_Zn','Q_Y','Q_Ba']
		#prior
		error_parameters = [ -3., 2., 0.05]#, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]#, 0.05, 0.05]#, 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. ]
		#start_value
		#error_parameters = [ -6.7, 2., 0.14, 0.1, 0.09, 0.16, 0.13, 0.04, 0.04, 0.15, 0.13, 0.09]#, 0.05, 0.05]#, 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. ]
		error_parameters_to_optimize = ['lnodds','Q','astronomical_spread']#'sigma_O','sigma_Mg','sigma_Al','sigma_Si','sigma_P','sigma_Ca','sigma_Cr','sigma_Mn','sigma_Fe','sigma_Ni']#,'sigma_Y','sigma_Ba']#,'Q_O','Q_Na','Q_Mg','Q_Al','Q_Si','Q_Ca','Q_Ti','Q_Cr','Q_Fe','Q_Ni','Q_Zn','Q_Y','Q_Ba']
		#error_parameters = [ -2., 2., 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]#, 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. ]
		#error_parameters_to_optimize = ['lnodds','Q','sigma_O','sigma_Na','sigma_Mg','sigma_Al','sigma_Si','sigma_Ca','sigma_Ti','sigma_Cr','sigma_Fe','sigma_Ni','sigma_Zn','sigma_Y','sigma_Ba']#,'sigma_Y','sigma_Ba']#,'Q_O','Q_Na','Q_Mg','Q_Al','Q_Si','Q_Ca','Q_Ti','Q_Cr','Q_Fe','Q_Ni','Q_Zn','Q_Y','Q_Ba']
	else:
		error_parameters = []
		error_parameters_to_optimize = []			
	# Prior [-2.3 ,0.005 ,0.2 ,1.0 ,3.3,1.0, 3.5,0.5]
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
		ISM_parameters =  [-0.3, 3.5,	0.5, 0.3]#,0.2]#, 0.7, 0.3, 0.0]
		ISM_parameters_to_optimize = ['log10_starformation_efficiency', 'sfr_scale', 'outflow_feedback_fraction','log10_gas_reservoir_mass_factor']#,'log10_sfr_factor_for_cosmic_accretion']#,'log10_gas_reservoir_mass_factor','log10_a_parameter','log10_gas_power']
	else:
		ISM_parameters = []
		ISM_parameters_to_optimize = []
	assert len(ISM_parameters) == len(ISM_parameters_to_optimize)

	if False:
		yield_factors = [1.,1.,1., 1.,1.,1.,1.,1.,1.,1.,1.]
		yields_factors_to_optimize = ['SN2_Mg','SN2_Si','SN2_S','SN2_Ca','SN2_C','SN2_N','SN2_Al','SN2_Na', 'SN2_K','SN2_Mn','SN2_Ni']
		yield_factors = [1.,1.,1., 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1., 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]
		yields_factors_to_optimize = ['SN2_O','SN2_Na','SN2_Mg','SN2_Al','SN2_Si','SN2_Ca','SN2_Ti','SN2_Cr', 'SN2_Mn', 'SN2_Fe','SN2_Co','SN2_Ni','SN2_Zn','SN1a_O','SN1a_Na','SN1a_Mg','SN1a_Al','SN1a_Si','SN1a_Ca','SN1a_Ti','SN1a_Cr', 'SN1a_Mn','SN1a_Fe','SN1a_Co','SN1a_Ni','SN1a_Zn']

	else:
		yield_factors = []
		yields_factors_to_optimize = []
	assert len(yield_factors) == len(yields_factors_to_optimize)
	# Prior [0.0,0.,0., 0.,0.,0.,0.,0.,0.,0.,0.]
	if False:
		offsets = [0.0,0.0,0., 0.,0.,0.,0.,0.,0.,0.,0.]
		offsets_to_optimize = ['offset_Mg','offset_Si','offset_S','offset_Ca','offset_C','offset_N','offset_Al','offset_Na', 'offset_K','offset_Mn','offset_Ni']
	else:
		offsets = []
		offsets_to_optimize = []
	assert len(offsets) == len(offsets_to_optimize) 
	


	p0 = np.hstack((error_parameters,SSP_parameters,ISM_parameters,yield_factors,offsets))
	if make_multi_zone:
		for i in range(len(multi_zone_list)):
			p0 = np.hstack((p0,ISM_parameters))
	to_optimize = np.array(error_parameters_to_optimize + SSP_parameters_to_optimize + ISM_parameters_to_optimize+yields_factors_to_optimize+offsets_to_optimize)
	#p0 = np.load('mcmc/best_parameter_values.npy')[0]
	#assert len(p0) == len(to_optimize)
	ndim = len(to_optimize)
	#print ndim

	constraints = {
	'high_mass_slope' : (-4.,-1.),
	'log10_N_0' : (None,0), 
	'log10_sn1a_time_delay' : (None,1.),
	'log10_starformation_efficiency' : (None,None),
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
	'sigma_error_population' : (0.,10),
	'fraction_error_population' : (0.,1.),
	'gas_reservoir_mass_factor' : (0.,20.),
	'infall_scale' : (0.0,end),
	'sn1a_norm' : (0.,None),
	'sn1a_scale' : (0.,None),
	'sigma_astro' : (0.,None),
	'epsilon' : (0.,1.),



	'SN2_C'  : (0,150),
	'SN2_N'  : (0,150),
	'SN2_O'  : (0,150),
	'SN2_Na' : (0,150),
	'SN2_Mg' : (0,150),
	'SN2_Al' : (0,150),
	'SN2_Si' : (0,150),
	'SN2_S'  : (0,150),
	'SN2_P'  : (0,150),
	'SN2_K'  : (0,150),
	'SN2_Ca' : (0,150),
	'SN2_Ti' : (0,150),
	'SN2_V'  : (0,150),
	'SN2_Cr' : (0,150),
	'SN2_Co' : (0,150),
	'SN2_Mn' : (0,150),
	'SN2_Fe' : (0,150),
	'SN2_Zn' : (0,150),
	'SN2_Ni' : (0,150),


	'SN1a_C'  : (0,150),
	'SN1a_N'  : (0,150),
	'SN1a_O'  : (0,150),
	'SN1a_Na' : (0,150),
	'SN1a_Mg' : (0,150),
	'SN1a_Al' : (0,150),
	'SN1a_Si' : (0,150),
	'SN1a_S'  : (0,150),
	'SN1a_P'  : (0,150),
	'SN1a_K'  : (0,150),
	'SN1a_Ca' : (0,150),
	'SN1a_Ti' : (0,150),
	'SN1a_V'  : (0,150),
	'SN1a_Cr' : (0,150),
	'SN1a_Co' : (0,150),
	'SN1a_Mn' : (0,150),
	'SN1a_Fe' : (0,150),
	'SN1a_Zn' : (0,150),
	'SN1a_Ni' : (0,150),

	'agb_C'  : (0,150),
	'agb_N'  : (0,150),
	'agb_O'  : (0,150),
	'agb_Na' : (0,150),
	'agb_Mg' : (0,150),
	'agb_Al' : (0,150),
	'agb_Si' : (0,150),
	'agb_S'  : (0,150),
	'agb_Mn' : (0,150),
	'agb_Fe' : (0,150),
	'agb_Ni' : (0,150),

	'offset_O'	: (-1,1), 
	'offset_Mg' : (-1,1),
	'offset_Si'	: (-1,1),
	'offset_P'  : (-1,1),
	'offset_S'	: (-1,1),
	'offset_Ca'	: (-1,1),
	'offset_C'  : (-1,1),
	'offset_N'  : (-1,1),
	'offset_Al' : (-1,1),
	'offset_Na' : (-1,1),
	'offset_K'  : (-1,1),
	'offset_Cr' : (-1,1),
	'offset_Mn' : (-1,1),
	'offset_Fe' : (-1,1),
	'offset_Ni' : (-1,1),
	'offset_C+N': (-1,1),
	'Q_O' : (1,None),
	'Q_Na': (1,None),
	'Q_Mg': (1,None),
	'Q_Al': (1,None),
	'Q_Si': (1,None),
	'Q_Ca': (1,None),
	'Q_Ti': (1,None),
	'Q_Cr': (1,None),
	'Q_Fe': (1,None),
	'Q_Ni': (1,None),
	'Q_Zn': (1,None),
	'Q_Y' : (1,None),
	'Q_Ba': (1,None),
	'Q' 	: (1,None),
	'lnodds': (None,0.),
	'astronomical_spread': (0.,1.),
	'sigma_O' : (0.,1.),
	'sigma_Na': (0.,1.),
	'sigma_Mg': (0.,1.),
	'sigma_Al': (0.,1.),
	'sigma_Si': (0.,1.),
	'sigma_Ca': (0.,1.),
	'sigma_Ti': (0.,1.),
	'sigma_Cr': (0.,1.),
	'sigma_Fe': (0.,1.),
	'sigma_Ni': (0.,1.),
	'sigma_Zn': (0.,1.),
	'sigma_Y' : (0.,1.),
	'sigma_Ba': (0.,1.),
	'sigma_P': (0.,1.),
	'sigma_K' : (0.,1.),
	'sigma_Mn': (0.,1.),
	}
	# the prior entry is (mean,std,0)
	# functional form 0 is a gaussian with log values and 1 is for fractions where the sigma distances are in factors from the mean (see cem_function.py)
	# for functional form 1 read (mean,factor,1)
	priors = {
	## gaussian priors
	'high_mass_slope' : (-2.29,0.2,0),	
	'log10_N_0' : (-2.75,0.3,0),
	'log10_sn1a_time_delay' : (-0.8,0.3,0),
	'log10_starformation_efficiency' : (-0.3,0.3,0),
	'sfr_scale' : (3.5,1.5,0),
	'outflow_feedback_fraction' : (0.5,0.2,0),
	'log10_gas_reservoir_mass_factor' : (0.3,0.3,0),

	'a_parameter' : (3.,3.,0),
	'infall_scale' : (3.3,0.5,0),
	'lnodds': (-3.,0.2,0),
	'gas_power': (1.5,0.2,0),
	
	
	

	'log10_sfr_factor_for_cosmic_accretion': (0.2,0.3,0),
	'log10_a_parameter' : (0.3,0.2,0),
	'log10_gas_power' : (0,0.15,0),

	'offset_O'	: (0,0.2,0),
	'offset_Mg' : (0,0.2,0),
	'offset_Si'	: (0,0.2,0),
	'offset_P'  : (0,0.2,0),
	'offset_S'	: (0,0.2,0),
	'offset_Ca'	: (0,0.2,0),
	'offset_Cr' : (0,0.2,0),
	'offset_C'  : (0,0.2,0),
	'offset_N'  : (0,0.2,0),
	'offset_Al' : (0,0.2,0),
	'offset_Na' : (0,0.2,0),
	'offset_K'  : (0,0.2,0),
	'offset_Mn' : (0,0.2,0),
	'offset_Fe' : (0,0.2,0),
	'offset_Ni' : (0,0.2,0),
	'offset_C+N': (0,0.2,0),
	
	## Priors on factors
	'starformation_efficiency' : (0.5,3.,1), 
	'mass_factor' : (1.,1.2,1),
	'norm_infall' : (1.,1.2,1),
	'sn1a_time_delay' : (0.3,3.,1),
	'N_0' : (0.001,3.,1),
	'gas_at_start' : (0.1,2.,1),
	'gas_reservoir_mass_factor' : (3.,2.,1),
	'Q' 	: (2.,1.2,1),
	'Q_O'   : (2.,1.2,1),
	'Q_Na'  : (2.,1.2,1),
	'Q_Mg'  : (2.,1.2,1),
	'Q_Al'  : (2.,1.2,1),
	'Q_Si'  : (2.,1.2,1),
	'Q_Ca'  : (2.,1.2,1),
	'Q_Ti'  : (2.,1.2,1),
	'Q_Cr'  : (2.,1.2,1),
	'Q_Fe'  : (2.,1.2,1),
	'Q_Ni'  : (2.,1.2,1),
	'Q_Zn'  : (2.,1.2,1),
	'Q_Y'   : (2.,1.2,1),
	'Q_Ba'  : (2.,1.2,1),


	'astronomical_spread': (0.05,2.,1),
	'sigma_O'   : (0.05,2.,1),
	'sigma_Na'  : (0.05,2.,1),
	'sigma_Mg'  : (0.05,2.,1),
	'sigma_Al'  : (0.05,2.,1),
	'sigma_Si'  : (0.05,2.,1),
	'sigma_Ca'  : (0.05,2.,1),
	'sigma_Ti'  : (0.05,2.,1),
	'sigma_Cr'  : (0.05,2.,1),
	'sigma_Fe'  : (0.05,2.,1),
	'sigma_Ni'  : (0.05,2.,1),
	'sigma_Zn'  : (0.05,2.,1),
	'sigma_Y'   : (0.05,2.,1),
	'sigma_Ba'  : (0.05,2.,1),
	'sigma_P'  : (0.05,2.,1),
	'sigma_K'   : (0.05,2.,1),
	'sigma_Mn'  : (0.05,2.,1),
	
	'SN2_C'  : (1,2.,1),
	'SN2_N'  : (1,2.,1),
	'SN2_O'  : (1,2.,1),
	'SN2_Na' : (1,2.,1),
	'SN2_Mg' : (1,2.,1),
	'SN2_Al' : (1,2.,1),
	'SN2_Si' : (1,2.,1),
	'SN2_P'  : (1,2.,1),
	'SN2_S'  : (1,2.,1),
	'SN2_K'  : (1,2.,1),
	'SN2_Ca' : (1,2.,1),
	'SN2_Cr' : (1,2.,1),
	'SN2_Co' : (1,2.,1),
	'SN2_Ti' : (1,2.,1),
	'SN2_V'  : (1,2.,1),
	'SN2_Mn' : (1,2.,1),
	'SN2_Fe' : (1,2.,1),
	'SN2_Ni' : (1,2.,1),
	'SN2_Zn' : (1,2.,1),

	'SN1a_C'  : (1,2.,1),
	'SN1a_N'  : (1,2.,1),
	'SN1a_O'  : (1,2.,1),
	'SN1a_Na' : (1,2.,1),
	'SN1a_Mg' : (1,2.,1),
	'SN1a_Al' : (1,2.,1),
	'SN1a_Si' : (1,2.,1),
	'SN1a_P'  : (1,2.,1),
	'SN1a_S'  : (1,2.,1),
	'SN1a_K'  : (1,2.,1),
	'SN1a_Ca' : (1,2.,1),
	'SN1a_Cr' : (1,2.,1),
	'SN1a_Co' : (1,2.,1),
	'SN1a_Ti' : (1,2.,1),
	'SN1a_V'  : (1,2.,1),
	'SN1a_Mn' : (1,2.,1),
	'SN1a_Fe' : (1,2.,1),
	'SN1a_Ni' : (1,2.,1),
	'SN1a_Zn' : (1,2.,1),

	'agb_C'  : (1,2.,1),
	'agb_N'  : (1,2.,1),
	'agb_O'  : (1,2.,1),
	'agb_Na' : (1,2.,1),
	'agb_Mg' : (1,2.,1),
	'agb_Al' : (1,2.,1),
	'agb_Si' : (1,2.,1),
	'agb_S'  : (1,2.,1),
	'agb_Mn' : (1,2.,1),
	'agb_Fe' : (1,2.,1),
	'agb_Ni' : (1,2.,1),
	}
