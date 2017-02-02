import numpy as np 
from weighted_yield import SSP, lifetime_Argast, lifetime_Raiteri
from imf import IMF
from yields import SN2_feedback
from yields import AGB_feedback
from yields import SN1a_feedback
from yields import Hypernova_feedback

class SSP_wrap():
	def __init__(self, a):

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
		
		basic_ssp = SSP(False, float(z), np.copy(self.imf.x), np.copy(self.imf.dm), np.copy(self.imf.dn), np.copy(time_steps), list(elements), str(self.a.stellar_lifetimes), str(self.a.interpolation_scheme), bool(self.a.only_net_yields_in_process_tables))
		basic_ssp.sn2_feedback(list(self.sn2.elements), dict(self.sn2.table), np.copy(self.sn2.metallicities), float(self.a.sn2mmin), float(self.a.sn2mmax),list(element_fractions))
		basic_ssp.agb_feedback(list(self.agb.elements), dict(self.agb.table), list(self.agb.metallicities), float(self.a.agbmmin), float(self.a.agbmmax),np.hstack(element_fractions))
		basic_ssp.sn1a_feedback(list(self.sn1a.elements), list(self.sn1a.metallicities), dict(self.sn1a.table), str(self.a.time_delay_functional_form), float(self.a.sn1ammin), float(self.a.sn1ammax), self.a.sn1a_parameter, float(self.a.total_mass), bool(self.a.stochastic_IMF))
		# exposing these tables to the outside wrapper
		self.table = basic_ssp.table
		self.sn2_table = basic_ssp.sn2_table
		self.agb_table = basic_ssp.agb_table
		self.sn1a_table = basic_ssp.sn1a_table

def initialise_stuff(a):

	from solar_abundance import solar_abundances
	from sfr import SFR 
	from infall import INFALL

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
		element_fractions = []
		for item in elements_to_trace:
			element_fractions.append(float(np.copy(cube.cube[item][max(i-1,0)]/cube.cube['gas'][max(i-1,0)])))## gas element fractions from one time step before	
		metallicity = float(cube.cube['Z'][i])
		time_steps = np.copy(basic_sfr.t[:j])
		basic_ssp.calculate_feedback(float(metallicity), list(elements_to_trace), list(element_fractions), np.copy(time_steps))
		cube.advance_one_step(i+1,np.copy(basic_ssp.table),np.copy(basic_ssp.sn2_table),np.copy(basic_ssp.agb_table),np.copy(basic_ssp.sn1a_table))
	abundances,elements,numbers = mass_fraction_to_abundances(np.copy(cube.cube),np.copy(basic_solar.table))
	weights = cube.cube['sfr']
	abundances = append_fields(abundances,'weights',weights)
	abundances = np.array(abundances)

	return cube.cube, abundances, cube.gas_reservoir

def Chempy_gross(a):
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

	return cube.cube, abundances, cube.gas_reservoir