import numpy as np
from .making_abundances import abundance_to_mass_fraction_normed_to_solar,abundance_to_mass_fraction

class PRIMORDIAL_INFALL(object):
	def __init__(self,elements,solar_table):
		'''
		This class calculates the chemical abundance fractions and can be used to provide primordial or solar scaled material for infall or gas composition at the beginning of the simulation

		INPUT upon initialisation are a list of elements and the solar table from the solar_abundance class.

		The elements are actually sorted by their element number


		'''
		self.elements = np.hstack(elements)
		element_names = list(solar_table['Symbol'])
		element_number = []
		element_masses = []
		element_abundances = []
		for item in element_names:
			element_number.append(int(solar_table['Number'][np.where(solar_table['Symbol']==item)]))
			element_masses.append(solar_table['Mass'][np.where(solar_table['Symbol']==item)])
			element_abundances.append(solar_table['photospheric'][np.where(solar_table['Symbol']==item)])

		sorted_index = np.argsort(np.array(element_number))
		element_number = [element_number[i] for i in sorted_index]
		element_masses = [element_masses[i] for i in sorted_index]
		element_names = [element_names[i] for i in sorted_index]
		element_abundances = [element_abundances[i] for i in sorted_index]
		self.all_elements = np.hstack(element_names)
		self.numbers = np.hstack(element_number)
		self.masses = np.hstack(element_masses)  
		self.all_abundances = np.hstack(element_abundances)
		self.all_fractions = abundance_to_mass_fraction(self.all_elements,self.masses,self.all_abundances,self.all_abundances,self.all_elements)
	  

	def primordial(self,):
		'''
		This returns primordial abundance fractions.
		'''
		self.symbols = np.hstack(self.elements)
		self.fractions = np.zeros(len(self.symbols))
		self.fractions[np.where(self.symbols=='H')] = 0.76
		self.fractions[np.where(self.symbols=='He')] = 0.24
		
	def solar(self, metallicity_in_dex, helium_fraction = 0.285):
		'''
		solar values scaled to a specific metallicity

		INPUT 

		   metallicity_in_dex =
		   helium_fraction = 
		'''
		self.symbols = []
		self.abundances = []
		for i, item in enumerate(self.elements):
			self.abundances.append(0.)
			self.symbols.append(item)
		self.abundances = np.hstack(self.abundances)
		self.symbols = np.hstack(self.symbols)
		self.fractions = abundance_to_mass_fraction_normed_to_solar(self.all_elements,self.masses,self.all_abundances,self.abundances,self.symbols)
		divisor = np.power(10,float(metallicity_in_dex))
		tmp = 0.
		for i,item in enumerate(self.symbols):
			if item not in ['H','He']:
				self.fractions[i] *= divisor
				tmp += self.fractions[i]
			if item in ['He']:
				self.fractions[i] = helium_fraction
				tmp += self.fractions[i]
		self.fractions[np.where(self.symbols == 'H')] = 1-tmp
		

	def sn2(self,paramet):
		'''
		This can be used to produce alpha enhanced initial abundances

		the fractions of the CC SN feedback and the iron abundance in dex needs to be specified
		'''
		sn2_fractions,iron_dex = paramet
		self.symbols = self.elements 

		solar_iron_fraction = self.all_fractions[np.where(self.all_elements == 'Fe')]
		scaled_iron_fraction = np.power(10,iron_dex)*solar_iron_fraction
		iron_fraction_in_sn2_feedback = sn2_fractions[np.where(self.elements == 'Fe')]
		self.fractions = np.divide(sn2_fractions,iron_fraction_in_sn2_feedback/scaled_iron_fraction)
		solar_helium_fraction = self.all_fractions[np.where(self.all_elements == 'He')]
		self.fractions[np.where(self.elements == 'He')] = solar_helium_fraction
		self.fractions[np.where(self.elements == 'H')] = 1 - sum(self.fractions[np.where(self.elements != 'H')])
		self.z = sum(self.fractions[np.where(np.logical_and(self.elements != 'H',self.elements != 'He'))])

class INFALL(object):
	'''
	This class provides the infall mass over time and is matched to the SFR class
	'''
	def __init__(self, t, sfr):
		"""
		Upon initialisation the timesteps and the sfr need to be provided

		Input:

		   t = timesteps
		
		   sfr = the SFR for the timesteps
		
		t is time in Gyr over which infall takes place as a numpy array
		sfr is the star formation rate
		"""
		self.t = t
		self.sfr = sfr

	def constant(self, paramet = 1):
		"""Constant gas infall of amount in Msun/pc^2/Gyr (default is 1)
		For test purposes only."""
		amount = paramet
		self.infall = np.zeros(len(self.t)) + amount

	def linear(self, paramet = (6.3, -0.5)):
		"""Linear gas infall rate (usually decreasing) in Msun/pc^2/Gyr
		with an initial infall rate of start (default 6.5)
		and a decrease/increase of slope * t from above (default -0.5)"""
		start, slope = paramet


	def polynomial(self, paramet = ([-0.003, 0.03, -0.3, 5.])):
		"""Polynomial gas infall rate in Msun/pc^2/Gyr.
		coeff: 1D array of coefficients in decreasing powers.
		The number of coeff given determines the order of the polynomial.
		Default is -0.004t^3 + 0.04t^2 - 0.4t + 6 for okay-ish results"""
		coeff = paramet
		poly = np.poly1d(coeff)

	def gamma_function(self,mass_factor = 1,a_parameter = 2, loc = 0, scale = 3):
		'''
		the gamma function for a_parameter = 2 and loc = 0 produces a peak at scale so we have a two parameter sfr.
		Later we can also release a to have a larger diversity in functional form.
		'''
		from scipy.stats import gamma
		self.infall = gamma.pdf(self.t,a_parameter,loc,scale)
		norm = sum(self.sfr)*mass_factor
		self.infall = np.divide(self.infall*norm,sum(self.infall))


	def exponential(self, paramet = ( -0.24, 0., 1.)):
		"""
		Exponential gas infall rate in Msun/pc^2/Gyr.
		The exponent is b * t + c, whole thing shifted up by d and normalised by e to the SFR.
		Default is b = -0.15 and e = 1, rest 0
		"""
		b, d, e = paramet
		sfr_norm = e
		norm = sum(self.sfr)*sfr_norm
		self.infall = np.exp(b * self.t ) + d
		self.infall = np.divide(self.infall*norm,sum(self.infall))
 
	def sfr_related(self, ):
		'''
		the infall will be calculated during the Chempy run according to the star formation efficiency usually following a Kennicut-Schmidt law
		'''
		self.infall = np.zeros_like(self.sfr)
