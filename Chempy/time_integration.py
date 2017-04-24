import numpy as np 

class ABUNDANCE_MATRIX(object):
	'''
	This class contains all information necessary to characterize the chemical evolution of the open-box one-zone Chempy model.
	
	It calculates the mass flow between the different components. And can advance the chemical evolution when the enrichment from the SSP is provided.
	'''
	def __init__(self, time, sfr, infall, list_of_elements,infall_symbols,infall_fractions,gas_at_start,gas_at_start_symbols,gas_at_start_fractions,gas_reservoir_mass_factor,outflow_feedback_fraction,check_processes,starformation_efficiency,gas_power, sfr_factor_for_cosmic_accretion, cosmic_accretion_elements, cosmic_accretion_element_fractions):
		'''
		Upon initialization the provided information is stored and initial conditions as provided by the other chempy classes are calculated.
		
		The class assigns all the information into different tables that can be queeried.
		
		Most importantly self.cube contains the ISM evolution self.gas_reservoir the corona evolution and self.sn2/sn1a/agb_feedback the enrichment from the individual nucleosynthetic processes.

		INPUT:
		
		   time = the time-steps
		
		   sfr = the corresponding star formation rate
		
		   infall = the infall at that time (can be 0 when 'sfr_related' is chosen)
		
		   list_of_elements = which elements to trace (a list of Symbols)
		
		   infall_symbols = list of element symbols for the infall
		
		   infall_fractions = the corresponding elemental fractions
		
		   gas_at_start = how much gas at start do we have (default = 0)
		
		   gas_at_start_symbols = list of elements at beginnin
		
		   gas_at_start_fractions = the corresponding fractions
		
		   gas_reservoir_mass_factor = how much more mass does the corona has compared to the integrated SFR
		
		   outflow_feedback_fraction = how much enrichment goes into the corona (in fraction, the rest goes into the ISM)
		
		   check_processes = boolean, should the individual nucleosynthetic processes be tracked (usually not necessary during the MCMC but useful when exporing a single parameter configuration)
		
		   starformation_efficiency = the SFE for a linear Kennicut-Schmidt law
		
		   gas_power = The Schmidt_exponent (default = 1, i.e. linear)
		
		   sfr_factor_for_cosmic_accretion = how much more gas should be infalling in the corona compared to the SFR
		
		   cosmic_accretion_elements = element list of this cosmic infall
		
		   cosmic_accretion_element_fractions = the corresponding fractions (all element fractions are usually primordial)

		OUTPUT:
		
		a few structured arrays will be created most notably:
		
		   .cube = ISM evolution
		
		   .gas_reservoir = corona evolution
		
		   .sn2_feedback = CC-SN feedback
		
		   .sn1a_feedback = SN Ia feedback
		
		   .agb_feedback = AGB feedback
		'''
		self.time = time
		self.dt = time[1] - time[0]
		self.sfr = sfr#np.divide(sfr,sum(sfr))
		#self.sfr[0] = 0.41
		self.infall = infall
		self.elements = list_of_elements
		self.additional = ['sfr','infall','time','feedback','mass_in_remnants','stars','gas','Z','alpha','sn1a','sn2','pn','bh','hn']
		self.names =  self.additional + self.elements
		self.base = np.zeros(len(time))
		self.infall_symbols = infall_symbols
		self.infall_fractions = infall_fractions
		## fractions of the corona gas (gas reservoir) at start
		self.gas_at_start = gas_at_start * sum(sfr)#*self.dt #now normalised to total sfr#* (1./self.dt)
		self.gas_at_start_symbols = gas_at_start_symbols
		self.gas_at_start_fractions = gas_at_start_fractions
		## fractions of the ISM at start
		self.outflow_feedback_fraction = outflow_feedback_fraction
		self.check_processes = check_processes
		self.starformation_efficiency = starformation_efficiency * self.dt
		self.gas_power = gas_power
		self.sfr_factor_for_cosmic_accretion = sfr_factor_for_cosmic_accretion
		self.cosmic_accretion_elements = cosmic_accretion_elements
		self.cosmic_accretion_element_fractions = cosmic_accretion_element_fractions
		## fractions of the cosmic inflow into the corona gas (is the same at all times)
		list_of_arrays = []
		for i in range(len(list_of_elements)+len(self.additional)):
			list_of_arrays.append(self.base)
		self.cube = np.core.records.fromarrays(list_of_arrays,names=self.names)
		if self.check_processes:
			self.process_feedback_names = ['kinetic_energy','number_of_events','mass_in_remnants'] + self.elements
			self.sn2_cube = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)
			self.sn1a_cube = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)
			self.agb_cube = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)
		self.gas_reservoir = np.core.records.fromarrays(list_of_arrays,names=self.names)
		# Setting up the table for the gas and feedback composition
		self.cube['time'] = time
		self.cube['sfr'] = sfr
		if gas_at_start >= 0.00000001:
			self.cube['gas'][0] = self.gas_at_start
			for i,item in enumerate(self.gas_at_start_symbols):
				self.cube[item][0] = self.gas_at_start_fractions[i]*self.cube['gas'][0]	
		self.cube['infall'] = infall
		if self.starformation_efficiency != 0.:
			'''
			this applies when using the Kennicut-Schmidt law (infall = 'sfr-related'):
			Infall at start is overwritten to the value required by the sfr with a specific starformation efficiency and infall should at least be as big as sfr
			'''
			gas_needed = max(self.sfr[0] * 1.0000001 ,np.power(self.sfr[0] / float(self.starformation_efficiency), 1./float(self.gas_power) )) ## the factor 1.00000001 is added because otherwise no mass for mixing will be left resulting in errors
			self.cube['infall'][0] = gas_needed
			self.infall[0] = gas_needed
		# NEW PRESCRIPTION (separating the abundance fractions from infall and sfr such that first infall occurs and then stars are formed from that material
		for i,item in enumerate(self.infall_symbols):
			self.cube[item][0] += self.infall_fractions[i]*self.infall[0]		
		self.cube['gas'][0] += self.infall[0] 
		gas_mass_temp = float(self.cube['gas'][0])
		for i,item in enumerate(self.elements):
			if gas_mass_temp == 0.:
				self.cube[item][0] = 0.
				assert self.sfr[0] == 0.
			else:
				self.cube[item][0] -= (self.cube[item][0]/gas_mass_temp)*self.sfr[0]
		self.cube['gas'][0] -= self.sfr[0]
		self.cube['stars'][0] = self.sfr[0]

		self.cube['feedback'][0] = 0.
		self.cube['mass_in_remnants'][0] = 0.		

		self.gas_reservoir['infall'] = self.infall
		self.gas_reservoir['time'] = self.time
		self.gas_reservoir['sfr'] = self.sfr
		
		cosmic_inflow = self.sfr[0] * self.sfr_factor_for_cosmic_accretion
		self.gas_reservoir['gas'][0] = cosmic_inflow
		for i,item in enumerate(self.cosmic_accretion_elements):
			self.gas_reservoir[item][0] = self.cosmic_accretion_element_fractions[i] * cosmic_inflow

		starting_gas = sum(self.sfr) * gas_reservoir_mass_factor - self.infall[0]
		self.gas_reservoir['gas'][0] += starting_gas
		for i,item in enumerate(self.infall_symbols):
			self.gas_reservoir[item][0] += self.infall_fractions[i] * starting_gas	
	
		for i,item in enumerate(self.elements):
			if item not in ['H','He']:
				self.cube['Z'][0] += self.cube[item][0]
				self.gas_reservoir['Z'][0] += self.gas_reservoir[item][0]
		self.cube['Z'][0] = self.cube['Z'][0] / self.cube['gas'][0]	
		self.gas_reservoir['Z'][0] = self.gas_reservoir['Z'][0] / self.gas_reservoir['gas'][0]	
		
		self.cube['alpha'][0] = 0

		self.all_feedback_names = ['mass_in_remnants','sn1a','sn2','pn','bh'] + self.elements
		base = np.zeros((len(self.time),len(self.time))) 
		list_of_arrays = []
		for i in range(len(self.all_feedback_names)):
			list_of_arrays.append(base)
		self.all_feedback = np.core.records.fromarrays(list_of_arrays,names=self.all_feedback_names)

		if self.check_processes:			
			base = np.zeros((len(self.time),len(self.time)))
			list_of_arrays = []
			for i in range(len(self.process_feedback_names)):
				list_of_arrays.append(base) 
			self.sn2_feedback = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)
			self.sn1a_feedback = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)
			self.agb_feedback = np.core.records.fromarrays(list_of_arrays,names=self.process_feedback_names)

		
	def advance_one_step(self,index,ssp_yield,sn2_yield,agb_yield,sn1a_yield):
		'''
		This method advances the chemical evolution one time-step.

		INPUT: 
		
		   index = which time step should be filled up
		
		   ssp_yield = yield of the ssp
		
		   sn2_yield = yield of sn2 only
		
		   agb_yield = yield of agb only
		
		   sn1a_yield = yield of sn1a only
		'''
		### This aligns the SSP yield such that it becomes a simple vector multiplication with a little memory overhead
		# self.all_feedback has the following data structure: [element][time_index_of_the_born_ssp][time_index_of_this_ssp_giving_back_the_elements_mass]

		for i,item in enumerate(self.all_feedback_names):
			self.all_feedback[item][index-1][index-1:] = ssp_yield[item] 
		if self.check_processes:
			for i,item in enumerate(self.process_feedback_names):
				self.sn2_feedback[item][index-1][index-1:] = sn2_yield[item] 
				self.sn1a_feedback[item][index-1][index-1:] = sn1a_yield[item]
				self.agb_feedback[item][index-1][index-1:] = agb_yield[item]
		feedback_mass = []
		for i,item in enumerate(self.all_feedback_names):
			tmp = self.sfr[:index] * self.all_feedback[item][:index,index]
			self.cube[item][index] = self.cube[item][index-1] + (1. - self.outflow_feedback_fraction) * sum(tmp)
			if item not in ['mass_in_remnants','sn1a','sn2','pn','bh']:
				feedback_mass.append(sum(tmp))
				self.gas_reservoir[item][index] = self.gas_reservoir[item][index-1] + (sum(tmp) * self.outflow_feedback_fraction)
		self.cube['stars'][index] = self.cube['stars'][index-1] - sum(feedback_mass) 
		self.cube['feedback'][index] = sum(feedback_mass)		
		if self.check_processes:
			for i,item in enumerate(self.process_feedback_names):
				tmp_sn2 = self.sfr[:index] * self.sn2_feedback[item][:index,index]
				self.sn2_cube[item][index] = sum(tmp_sn2)
				tmp_sn1a = self.sfr[:index] * self.sn1a_feedback[item][:index,index]
				self.sn1a_cube[item][index] = sum(tmp_sn1a)
				tmp_agb = self.sfr[:index] * self.agb_feedback[item][:index,index]
				self.agb_cube[item][index] = sum(tmp_agb)


		### First add the cosmic inflow to the corona
		cosmic_inflow = self.sfr[index] * self.sfr_factor_for_cosmic_accretion
		self.gas_reservoir['gas'][index] =  self.gas_reservoir['gas'][index-1] + cosmic_inflow
		for i,item in enumerate(self.cosmic_accretion_elements):
			self.gas_reservoir[item][index] += self.cosmic_accretion_element_fractions[i] * cosmic_inflow
		### part of the feedback goes into the gas_reservoir
		self.gas_reservoir['gas'][index] += (sum(feedback_mass) * self.outflow_feedback_fraction)
		self.gas_reservoir['feedback'][index] = (sum(feedback_mass) * self.outflow_feedback_fraction)
		#### For infall related to sfr (Schmidt law) calculating the infall.
		if self.starformation_efficiency != 0.:
			gas_needed = np.power(self.sfr[index] / float(self.starformation_efficiency),1./float(self.gas_power))
			gas_there = sum(list(self.cube[self.elements][index]))
			infall_needed = gas_needed - gas_there
			infall_needed *= 1.00000001 # to avoid less gas being requested than needed due to rounding errors (not sure what results from that, too little gas in the corona could be a result. lets see)
			if infall_needed < 0. :
				infall_needed = 0.
			## for few parameter values of gas_power and SFE the infall_needed value could be too small
			if infall_needed + gas_there <= self.sfr[index]:
				print('too few gas requested', 'infall needed= ', infall_needed, 'gas there = ', gas_there, 'total SFR = ', self.sfr, 'gas needed = ', gas_needed, 'corona = ', self.gas_reservoir['gas'][index], 'sfe = ', self.starformation_efficiency , 'dt = ', self.dt) 
				infall_needed = self.sfr[index] - gas_there
				infall_needed *= 1.01 # to avoid the ISM being empty
			if infall_needed > self.gas_reservoir['gas'][index]:
				print('gas reservoir is empty')
				infall_needed = float(self.gas_reservoir['gas'][index])
			self.infall[index] = float(infall_needed)
			self.cube['infall'][index] = float(infall_needed)
		## gas reservoir gas is taken away and infalling on cube_gas
		for i,item in enumerate(self.elements):
			self.cube[item][index] += self.infall[index] * np.divide(self.gas_reservoir[item][index],self.gas_reservoir['gas'][index])
		for i,item in enumerate(self.elements):
			self.gas_reservoir[item][index] -= self.infall[index] * np.divide(self.gas_reservoir[item][index],self.gas_reservoir['gas'][index])
		self.gas_reservoir['gas'][index] -= self.infall[index]

		# sfr will be subtracted in the next step self.sfr[index]
		self.cube['gas'][index] = sum(list(self.cube[self.elements][index]))
		assert self.cube['gas'][index] >= self.sfr[index], ('time index: ', index, 'gas: ', self.cube['gas'][index], 'sfr: ', self.sfr[index], 'total SFR: ', self.sfr, 'gas needed = ', gas_needed, 'corona = ', self.gas_reservoir['gas'][index], 'sfe = ', self.starformation_efficiency , 'dt = ', self.dt)
		for i,item in enumerate(self.elements):
			self.cube[item][index] -= self.sfr[index] * np.divide(self.cube[item][index],self.cube['gas'][index])
		self.cube['gas'][index] -= self.sfr[index]
		self.cube['stars'][index] += self.sfr[index]
		


		# determine metal fraction
		for i,item in enumerate(self.elements):
			if item not in ['H','He']:
				self.cube['Z'][index] += self.cube[item][index]
				self.gas_reservoir['Z'][index] += self.gas_reservoir[item][index]
		self.cube['Z'][index] = self.cube['Z'][index] / float(self.cube['gas'][index])
		self.gas_reservoir['Z'][index] = self.gas_reservoir['Z'][index] / float(self.gas_reservoir['gas'][index])		
		# determine alpha enhancement
		tmp = 0.
		for i,item in enumerate(self.elements):
			if item in ['O','Mg','Si','S','Ca','Ti']:
				self.cube['alpha'][index] += self.cube[item][index]
			#if item in ['Sc','V','Cr','Mn','Fe','Co','Ni','Zn']:
				#tmp += self.cube[item][index]
		self.cube['alpha'][index] = self.cube['alpha'][index] / self.cube['gas'][index]
		#print index, len(ssp_yield),self.cube['Z'][index],self.cube['alpha'][index]
