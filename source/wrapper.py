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
		
		basic_ssp = SSP(False, z, np.copy(self.imf.x), np.copy(self.imf.dm), np.copy(self.imf.dn), np.copy(time_steps), list(elements), str(self.a.stellar_lifetimes), str(self.a.interpolation_scheme), bool(self.a.only_net_yields_in_process_tables))
		basic_ssp.sn2_feedback(list(self.sn2.elements), dict(self.sn2.table), np.copy(self.sn2.metallicities), float(self.a.sn2mmin), float(self.a.sn2mmax),list(element_fractions))
		basic_ssp.agb_feedback(list(self.agb.elements), dict(self.agb.table), list(self.agb.metallicities), float(self.a.agbmmin), float(self.a.agbmmax),np.hstack(element_fractions))
		basic_ssp.sn1a_feedback(list(self.sn1a.elements), list(self.sn1a.metallicities), dict(self.sn1a.table), str(self.a.time_delay_functional_form), float(self.a.sn1ammin), float(self.a.sn1ammax), self.a.sn1a_parameter, float(self.a.total_mass), bool(self.a.stochastic_IMF))
		# exposing these table to the outside wrapper
		self.table = basic_ssp.table
		self.sn2_table = basic_ssp.sn2_table
		self.agb_table = basic_ssp.agb_table
		self.sn1a_table = basic_ssp.sn1a_table