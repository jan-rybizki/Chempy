import numpy as np 



def imf_mass_fraction_non_nativ(imf_dm,imf_x,mlow,mup):
	'''
	Function to determine the mass fraction of the IMF between two masses

	INPUT:
	
        imf_dm = imf_class.dm
	
        imf_x = imf_class.x
	
        mlow = lower mass
	
        mup = upper mass

	
        OUTPUT:
	
        the mass fraction in this mass range
	'''
	cut = np.where(np.logical_and(imf_x>=mlow,imf_x<mup))
	fraction = sum(imf_dm[cut])
	
	return(fraction)

def lifetime_Argast(m,Z_frac):
	"""
	here we will calculate the MS lifetime of the star after Argast et al., 2000, A&A, 356, 873
	
        INPUT:

        m = mass in Msun
	
        Z_frac = fractions of metals of the stellar composition
	
        Z = metallicity in Zsun
	
        
        OUTPUT:

        returns the lifetime of the star in Gyrs
	"""
	solar_metallicity = 0.015
	Z = Z_frac / solar_metallicity
	lm = np.log10(m)
	a0 =  3.79 + 0.24*Z
	a1 = -3.10 - 0.35*Z
	a2 =  0.74 + 0.11*Z
	tmp = a0 + a1*lm + a2*lm*lm
	return np.divide(np.power(10,tmp),1000)

def lifetime_Raiteri(m,Z):
	"""
	INPUT:

	   m = mass in Msun
	
	   Z = metallicity in Zsun
	
	returns the lifetime of the star in Gyrs
	"""
	lm = np.log10(m)
	lz = np.log10(Z)
	a0 =  10.13 + 0.07547*lz - 0.008084*lz*lz
	a1 = -4.424 - 0.7939*lz - 0.1187*lz*lz
	a2 =  1.262 + 0.3385*lz + 0.05417*lz*lz
	tmp = a0 + a1*lm + a2*lm*lm
	return np.divide(np.power(10,tmp),1e9)

class SSP(object):
	'''
	The simple stellar population class can calculate the enrichment over time for an SSP from a few assumptions and input yield tables.
	'''
	def __init__(self,output, z, imf_x, imf_dm, imf_dn, time_steps, elements_to_trace,  stellar_lifetimes, interpolation_scheme,only_net_yields_in_process_tables):
		'''
		Upon initialisation the table for the SSP evolution and enrichment over time is created. Interpolation from yield tables in mass is made linear.

		INPUT:
		
		   output = bool, should output be plotted
		
		   z = metallicity of the SSP in Z, not normed to solar
		
		   imf_x = class_imf.x
		
		   imf_dm = class_imf.dm
		
		   imf_dn = class_imf.dn
		
		   time_steps = time_steps usually coming from class_sfr.t
		
		   elements_to_trace = which elements should be traced
		
		   stellar_lifetimes = the name of the stellar lifetimes function that should be used ('Raiteri_1996', 'Argast_2000')
		
		   interpolation_scheme = which interpolation in metallicity should be used for the yield table ('linear' or 'logarithmic')
		
		   only_net_yields_in_process_tables = Should the total yields or only the net yields be stored in the nucleosynthetic enrichment tables, bool

		OUTPUT:
		
		the ssp_class.table holds key values of the evolution of the SSP all normalised to a mass of unity (which is the starting mass of the SSP).
		
		mass_in_ms_stars + cumsum(mass_of_ms_stars_dying) = 1 calculated from stellar lifetime with IMF 
		
		The element feedbacks are also given normalised to one. And the number of events as well
		'''
		self.z = z
		self.x = imf_x
		self.dx = self.x[1]-self.x[0]
		self.dm = imf_dm
		self.dn = imf_dn
		self.t = time_steps
		self.dt = time_steps[1] - time_steps[0]
		self.elements = elements_to_trace
		self.interpolation_scheme = interpolation_scheme
		self.stellar_lifetimes = stellar_lifetimes
		self.output = output
		self.net_yields = only_net_yields_in_process_tables
		lifetime_functions = {'Argast_2000': lifetime_Argast, 'Raiteri_1996': lifetime_Raiteri}
		
		######### Calculating and normalising the delay time function of SN1a
		if self.stellar_lifetimes in lifetime_functions:
			self.inverse_imf = np.interp(self.t,lifetime_functions[self.stellar_lifetimes](self.x[::-1],self.z),self.x[::-1])
		else:
			raise Exception("Lifetime function named '%s' not implemented" % self.stellar_lifetimes)
		additional_keys = ['mass_of_ms_stars_dying','mass_in_ms_stars','mass_in_remnants','sn2','sn1a','pn','bh','hydrogen_mass_accreted_onto_white_dwarfs']
		names = additional_keys + self.elements
		base = np.zeros(len(self.t))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		self.table = np.core.records.fromarrays(list_of_arrays,names=names)
		index = len(self.t)-1
		tmp = np.zeros(index)
		for i in range(index):
			mlow = self.inverse_imf[index-i]
			mup = self.inverse_imf[index-i-1]
			tmp[i] = imf_mass_fraction_non_nativ(self.dm,self.x,mlow,mup) 
		self.table['mass_of_ms_stars_dying'][1:] = tmp[::-1]
		self.table['mass_in_ms_stars'] = 1.-np.cumsum(self.table['mass_of_ms_stars_dying'])
		#print self.inverse_imf[1]

	def sn2_feedback(self, sn2_elements, sn2_yields, sn2_metallicities, sn2_mmin, sn2_mmax, fractions_in_gas):	
		'''
		Calculating the CC-SN feedback over time.
		The routine is sensitive to the ordering of the masses in the yield table, it must begin with the smallest increase to the biggest value.

		INPUT:
		
		   sn2_elements = which elements are provided by the yield table, list containing the symbols
		
		   sn2_yields = the yield table provided by Chempys SN2 yield class
		
		   sn2_metallicities = the metallicities of that table
		
		   sn2_mmin = the minimal mass of the CC-SN (default 8) in Msun
		
		   sn2_mmax = the maximum mass of the CC-SN (default 100) in Msun
		
		   fractions_in_gas = the birth material of the SSP (will be mixed into the enrichment as unprocessed material)

		OUTPUT:

		the classes table will be filled up and a sn2_table will be added (so that the individual processes can be tracked)
		'''
		# tracking the elemental feedback of individual processes
		additional_keys = ['kinetic_energy','number_of_events','mass_in_remnants']
		names = additional_keys + self.elements # not sure if all elements should be taken (might be easier to add the 3 tables in order to get total yield)
		base = np.zeros(len(self.t))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		self.sn2_table = np.core.records.fromarrays(list_of_arrays,names=names)
		

		self.sn2_elements = sn2_elements
		self.sn2 = sn2_yields
		self.sn2_mmin = sn2_mmin
		self.sn2_mmax = sn2_mmax
		self.sn2_metallicities = sn2_metallicities

		####### counting the sn2 events
		for i,item in enumerate(self.inverse_imf[:-1]):
			lower_cut = max(self.inverse_imf[i+1],self.sn2_mmin)
			upper_cut = min(self.inverse_imf[i],self.sn2_mmax)
			self.table['sn2'][i+1] += imf_mass_fraction_non_nativ(self.dn,self.x,lower_cut,upper_cut)
			self.sn2_table['number_of_events'][i+1] += imf_mass_fraction_non_nativ(self.dn,self.x,lower_cut,upper_cut)
			if upper_cut<lower_cut and self.inverse_imf[i+1]<self.sn2_mmin:
				break
		### interpolation of 2 yield sets with different metallicities
		self.sn2_metallicities = np.sort(self.sn2_metallicities)
		metallicity_list = []
		if len(self.sn2_metallicities) == 1:
			metallicity_list.append(self.sn2_metallicities[0])
		elif self.z < min(self.sn2_metallicities):
			metallicity_list.append(min(self.sn2_metallicities))
		elif self.z > max(self.sn2_metallicities):
			metallicity_list.append(max(self.sn2_metallicities))
		elif self.z in self.sn2_metallicities:
			metallicity_list.append(self.z)
		else:
			j=1
			while self.sn2_metallicities[j] < self.z:
				j += 1
			metallicity_list.append(self.sn2_metallicities[j-1])
			metallicity_list.append(self.sn2_metallicities[j])
		### the loop will be run through 2 times (if metallicity not outside or exactly at one of the precalculated metallicities) and the values of the two tables will be interpolated according to the prescribed function
		tables_to_interpolate = []
		if self.net_yields:
			net_tables_to_interpolate = []
		for s,metallicity_key in enumerate(metallicity_list):
			tables_to_interpolate.append(np.zeros_like(self.table))
			if self.net_yields:
				net_tables_to_interpolate.append(np.zeros_like(self.table))
			################################## loop which is metallicity independent		
			##### yield table is cut down such that only yields for masses between sn2_mmin and sn2_mmax are left in
			self.sn2[metallicity_key] = self.sn2[metallicity_key][np.where(np.logical_and(self.sn2[metallicity_key]['Mass']>=self.sn2_mmin,self.sn2[metallicity_key]['Mass']<=self.sn2_mmax))]

			self.sn2[metallicity_key] = np.sort(self.sn2[metallicity_key], order = 'Mass')[::-1]

			# tmp_masses holds the masses for which yields are calculated. 'weights' gives the mass fraction of the IMF that is dying and lies in the mass range of a specific yield mass
			tmp_masses = self.sn2[metallicity_key]['Mass']

			# Here the feedback per time step is calculated
			for time_index in range(len(self.table)-1):
				if self.inverse_imf[time_index] < self.sn2_mmin:
					break
				support = []
				support.append(self.sn2_mmax)
				for i in range(len(tmp_masses)-1): 
					support.append(np.mean(tmp_masses[i:i+2]))
				support.append(self.sn2_mmin)
				### support spans the mass ranges where the yields of a specific stellar mass (tmp_masses) will be applied

				weights = np.zeros(shape = len(support)-1)
				# This loop makes that the mass weights is only used for stars that are dying in this time-step
				for i in range(len(support)-1):
					lower_cut = max(support[i+1],self.inverse_imf[time_index+1])
					upper_cut = min(support[i],self.inverse_imf[time_index])
					if upper_cut < lower_cut:
						weights[i] = 0.
						continue
					cut = np.where(np.logical_and(self.x<upper_cut,self.x>=lower_cut))
					weights[i] = sum(self.dm[cut]) 
					# weights hold the weights with which each mass_range of the yield set needs to be multiplied
				###### here the feedback of the elements is calculated
				for item in list(set(self.elements).intersection(self.sn2_elements))+['mass_in_remnants']:
					tables_to_interpolate[s][item][time_index+1] = sum(self.sn2[metallicity_key][item]*weights)				
					if self.net_yields:
						net_tables_to_interpolate[s][item][time_index+1] = sum(self.sn2[metallicity_key][item]*weights)
				######## here the unprocessed wind feedback is calculated
				if 'unprocessed_mass_in_winds' in self.sn2[metallicity_key].dtype.names:
					for element_i,element_n in enumerate(self.elements):
						tables_to_interpolate[s][element_n][time_index+1] += sum(self.sn2[metallicity_key]['unprocessed_mass_in_winds']*weights*fractions_in_gas[element_i])
		metallicity_weight = []
		if len(metallicity_list) == 1:
			metallicity_weight.append(1.)
			# next line is the old way. might be much faster. But now also the elements in the wind are feed back into the ism.
			#for element_name in list(set(self.elements).intersection(self.sn2_elements))+['mass_in_remnants']:
			for element_name in list(self.elements)+['mass_in_remnants']:
				self.table[element_name] += tables_to_interpolate[0][element_name]
				if self.net_yields:
					self.sn2_table[element_name] += net_tables_to_interpolate[0][element_name]
				else:
					self.sn2_table[element_name] += tables_to_interpolate[0][element_name]
		else:
			if self.interpolation_scheme == 'linear':
				distance = metallicity_list[1]-metallicity_list[0]
				metallicity_weight.append((metallicity_list[1] - self.z) / float(distance))
				metallicity_weight.append((self.z - metallicity_list[0]) / float(distance))
			if self.interpolation_scheme == 'logarithmic':
				assert metallicity_list[1] != 0.
				if metallicity_list[0] == 0:
					metallicity_list[0] = 1e-10 ### This was a bug when using a yield table with 0 metallicity the log interpolation could not work.
				distance = np.log10(metallicity_list[1]) - np.log10(metallicity_list[0])
				metallicity_weight.append(( np.log10(metallicity_list[1]) - np.log10(self.z)) / float(distance))
				metallicity_weight.append(( np.log10(self.z) - np.log10(metallicity_list[0])) / float(distance))
			for i,item in enumerate(metallicity_list):
				for element_name in list(self.elements)+['mass_in_remnants']:
					self.table[element_name] += tables_to_interpolate[i][element_name] * metallicity_weight[i]
					if self.net_yields:
						self.sn2_table[element_name] += net_tables_to_interpolate[i][element_name] * metallicity_weight[i]
					else:
						self.sn2_table[element_name] += tables_to_interpolate[i][element_name] * metallicity_weight[i]
		### here the number of stars going sn2 is calculated

	def agb_feedback(self,agb_elements, agb_yields, agb_metallicities, agb_mmin, agb_mmax, fractions_in_gas):
		'''
		AGB enrichment calculation adds the feedback to the total SSP table and also to the self.agb_yield table.

		INPUT:
		
		   agb_elements = which elements are provided by the yield table, list containing the symbols
		
		   agb_yields = the yield table provided by Chempys AGB yield class
		
		   agb_metallicities = the metallicities of that table
		
		   agb_mmin = the minimal mass of the AGB stars (default 0.5) in Msun
		
		   agb_mmax = the maximum mass of the AGB stars (default 8) in Msun
		
		   fractions_in_gas = the birth material of the SSP (will be mixed into the enrichment as unprocessed material)

		OUTPUT:

		the classes table will be filled up and a sn2_table will be added (so that the individual processes can be tracked)
		'''

		# sensitive to the ordering of the masses from the yield. Should be checked in the beginning
		# ATTENTION this loop and interpolation scheme only works with fractional yields which should be used by default
		# metallicity interpolation is implemented
		# tmp_masses is the masses for which a yield is available
		# The breaks and weights do not need to be calculated every time. Only when the stellar lifetimes change significantly, though it does not take up too much time either
		# tracking the elemental feedback of individual processes

		additional_keys = ['kinetic_energy','number_of_events','mass_in_remnants']
		names = additional_keys + self.elements # not sure if all elements should be taken (might be easier to add the 3 tables in order to get total yield)
		base = np.zeros(len(self.t))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		self.agb_table = np.core.records.fromarrays(list_of_arrays,names=names)


		self.agb_mmin = agb_mmin
		self.agb_mmax = agb_mmax
		self.agb_elements = agb_elements
		self.agb = agb_yields
		self.agb_metallicities = np.sort(agb_metallicities)
		
		####### counting the pn events
		count_variable = 0
		for i,item in enumerate(self.inverse_imf):
			if item > agb_mmax:
				continue	
			elif count_variable == 0:
				self.table['pn'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.inverse_imf[i],self.agb_mmax)
				self.agb_table['number_of_events'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.inverse_imf[i],self.agb_mmax)
				count_variable += 1
			else:
				# The last time step is cut off with this method but it was important to add in order to be able to vary agb_mmin and if agb_mmin is below 1Msun effect is negligible
				if item < self.agb_mmin:
					#self.table['pn'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.agb_mmin,self.inverse_imf[i-1])
					#self.agb_table['number_of_events'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.agb_mmin,self.inverse_imf[i-1])
					break
				self.table['pn'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.inverse_imf[i],self.inverse_imf[i-1])
				self.agb_table['number_of_events'][i] += imf_mass_fraction_non_nativ(self.dn,self.x,self.inverse_imf[i],self.inverse_imf[i-1])
		metallicity_list = []
		if len(self.agb_metallicities) == 1:
			metallicity_list.append(self.agb_metallicities[0])
		elif self.z < min(self.agb_metallicities):
			metallicity_list.append(min(self.agb_metallicities))
		elif self.z > max(self.agb_metallicities):
			metallicity_list.append(max(self.agb_metallicities))
		elif self.z in self.agb_metallicities:
			metallicity_list.append(self.z)
		else:
			j=1
			while self.agb_metallicities[j] < self.z:
				j += 1
			metallicity_list.append(self.agb_metallicities[j-1])
			metallicity_list.append(self.agb_metallicities[j])
		### the loop will be run through 2 times (if metallicity not outside or exactly at one of the precalculated metallicities) and the values of the two tables will be interpolated according to the prescribed function
		tables_to_interpolate = []
		if self.net_yields:
			net_tables_to_interpolate = []
		for s,metallicity_key in enumerate(metallicity_list):
			tables_to_interpolate.append(np.zeros_like(self.table))
			if self.net_yields:
				net_tables_to_interpolate.append(np.zeros_like(self.table))
			################################## loop which is metallicity independent
			##### yield table is cut down such that only yields for masses between sn2_mmin and sn2_mmax are left in
			self.agb[metallicity_key] = self.agb[metallicity_key][np.where(np.logical_and(self.agb[metallicity_key]['Mass']>=self.agb_mmin,self.agb[metallicity_key]['Mass']<=self.agb_mmax))]

			tmp_masses = self.agb[metallicity_key]['Mass']
			# support defines ranges in which the yield could be asked. In fact it's defining the borders of each value in tmp_masses for which the yield of tmp_masses will be used
			#print tmp_masses
			support = []
			support.append(self.agb_mmax)
			for i in range(len(tmp_masses)-1): 
				support.append(np.mean(tmp_masses[i:i+2]))
			support.append(self.agb_mmin)
			support = np.array(support)
			#print support
			j=0
			last_item = 0.
			# the inverse_IMF (stars of mass ... are dying in this timestep) is going from top down. Each step is also a time_step
			mass_index_list = []
			mass_weight_list = []
			len_of_mass_weights = []
			##### loop to catch the mass weights and the mass indices for each time step
			for i,item in enumerate(self.inverse_imf):
				# was important to add in order to be able to vary agb_mmin
				if item < self.agb_mmin:
					break
				#gaps defines a temporal support for one time_step
				gaps = []
				#mass index gives the index of the mass from which the yield will be used in the range of 'gaps'
				mass_index = []
				gaps.append(last_item)
				#print item,support
				while support[j] > item: 
					gaps.append(support[j])
					mass_index.append(j-1)
					j += 1

				mass_index.append(j-1)
				last_item = item
				gaps.append(item)



				mass_weight = []
				# now this is the loop where the imf.dm is summed up i.e. the weight is calculated
				for t in range(len(gaps)-1):
					if mass_index[t] < 0:
						mass_weight.append(0.)####just added so that mass_weight and mass_indices have the same structure
						continue
					cut = np.where(np.logical_and(self.x>=gaps[t+1],self.x<gaps[t]))
					weight = sum(self.dm[cut])
					mass_weight.append(weight)

				#print i,j,support[j],item,gaps,mass_index

				len_of_mass_weights.append(len(mass_weight)) 
				mass_weight = np.array(mass_weight)
				mass_index = np.array(mass_index)
				mass_weight_list.append(mass_weight)
				mass_index_list.append(mass_index)

			
			# here the list structure is brought onto arrays to make vector multiplication possible for feedback calculations			
			max_different_masses_per_time_step = max(len_of_mass_weights)
			mass_index_array = np.zeros((len(self.inverse_imf),max_different_masses_per_time_step),dtype=int)
			mass_weight_array = np.zeros((len(self.inverse_imf),max_different_masses_per_time_step))

			for i in range(len(self.inverse_imf)):
				# was important to add in order to be able to vary agb_mmin
				if self.inverse_imf[i] < self.agb_mmin:
					break				
				for t in range(max_different_masses_per_time_step):
					if len(mass_weight_list[i])==t:
						break
					mass_weight_array[i][t] = mass_weight_list[i][t]
					mass_index_array[i][t] = mass_index_list[i][t]
			for element_name in list(set(self.elements).intersection(self.agb_elements))+['mass_in_remnants']:
				for t in range(max_different_masses_per_time_step):
					tables_to_interpolate[s][element_name] += self.agb[metallicity_key][element_name][mass_index_array[:,t]] * mass_weight_array[:,t]
					if self.net_yields:
						net_tables_to_interpolate[s][element_name] += self.agb[metallicity_key][element_name][mass_index_array[:,t]] * mass_weight_array[:,t]
					#print 'mass_index: ', element_name, t, mass_index_array[:,t]
			#print tables_to_interpolate[s][['O','C']]
			if 'unprocessed_mass_in_winds' in self.agb[metallicity_key].dtype.names:
				for element_i,element_n in enumerate(self.elements):
					for t in range(max_different_masses_per_time_step):	
						tables_to_interpolate[s][element_n] += (self.agb[metallicity_key]['unprocessed_mass_in_winds'][mass_index_array[:,t]] * mass_weight_array[:,t]) * fractions_in_gas[element_i]
			#print 'after adding wind' 
			#print tables_to_interpolate[s][['O','C']]



		########## end of loop which is metallicity independent
		#for i,item in enumerate(self.inverse_imf[:]):
		#	
		#	print i, item, mass_weight_array[i],mass_index_array[i]
		##### interpolating metallicity if there were two different used.
		metallicity_weight = []
		if len(metallicity_list) == 1:
			metallicity_weight.append(1.)
			for element_name in list(self.elements)+['mass_in_remnants']:
				self.table[element_name] += tables_to_interpolate[0][element_name]
				if self.net_yields:
					self.agb_table[element_name] += net_tables_to_interpolate[0][element_name]
				else:
					self.agb_table[element_name] += tables_to_interpolate[0][element_name]


		else:
			if self.interpolation_scheme == 'linear':
				distance = metallicity_list[1]-metallicity_list[0]
				metallicity_weight.append((metallicity_list[1] - self.z) / float(distance))
				metallicity_weight.append((self.z - metallicity_list[0]) / float(distance))
			if self.interpolation_scheme == 'logarithmic':
				assert metallicity_list[1] != 0.
				if metallicity_list[0] == 0:
					metallicity_list[0] = 1e-10 ### This was a bug when using a yield table with 0 metallicity the log interpolation could not work.
				distance = np.log10(metallicity_list[1]) - np.log10(metallicity_list[0])
				metallicity_weight.append(( np.log10(metallicity_list[1]) - np.log10(self.z)) / float(distance))
				metallicity_weight.append(( np.log10(self.z) - np.log10(metallicity_list[0])) / float(distance))
			for i,item in enumerate(metallicity_list):
				for element_name in list(self.elements)+['mass_in_remnants']:
					self.table[element_name] += tables_to_interpolate[i][element_name] * metallicity_weight[i]
					if self.net_yields:
						self.agb_table[element_name] += net_tables_to_interpolate[i][element_name] * metallicity_weight[i]
					else:
						self.agb_table[element_name] += tables_to_interpolate[i][element_name] * metallicity_weight[i]

	def sn1a_feedback(self, sn1a_elements, sn1a_metallicities, sn1a_yields, time_delay_functional_form, sn1a_min, sn1a_max, time_delay_parameter,ssp_mass, stochastic_IMF):
		'''
		Calculating the SN1a feedback over time

		INPUT:
		
		   sn1a_elements = Which elements are provided by the yield table
		
		   sn1a_metallicities = metallicities in the yield table 
		
		   sn1a_yields = yield table
		
		   time_delay_functional_form = which functional form of the delay time should be used ('normal','maoz','gamma_function'). Maoz is the default and the others are not tested. Check for functionality
		
		   sn1a_min = the minimum mass from which sn1a can occur (does not matter for maoz)
		
		   sn1a_max = the maximum mass from which SN Ia can occur (does not mater for maoz)
		
		   time_delay_parameter = a tuple containing the parameters for the specific functional form
		
		   ssp_mass = the mass of the SSP
		
		   stochastic_IMF = bool, do we want to use stochastci explosions

		
		OUTPUT:

		the classes table will be filled up and a sn2_table will be added (so that the individual processes can be tracked)

		for MAOZ functional form the following parameters are in time_delay_parameter:
		
		   N_0 = Number of SNIa exploding per Msun over the course of 15Gyr
		
		   tau_8 = The delay time when the first SN Ia explode (usually 40Myr are anticipated because then 8Msun stars start to die but our Prior is more at 160Myr)
		
		   s_eponent = the time decay exponent
		
		   dummy = not in use anymore
		'''
		end_of_time = 15 #Gyrs Over this time-span the SN1a explosions will be distributed, for mass normalisation reasons
		additional_keys = ['kinetic_energy','number_of_events','mass_in_remnants']
		names = additional_keys + self.elements # not sure if all elements should be taken (might be easier to add the 3 tables in order to get total yield)
		base = np.zeros(len(self.t))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		self.sn1a_table = np.core.records.fromarrays(list_of_arrays,names=names)


		if time_delay_functional_form == 'normal':
			pn_number_detonating_until_12Gyr,time_delay_peak,time_delay_time_scale,gauss_beginning = time_delay_parameter
		elif time_delay_functional_form == 'maoz':
			N_0,tau_8,s_exponent,dummy = time_delay_parameter
		elif time_delay_functional_form == 'gamma_function':
			number_factor, a_parameter, loc, scale = time_delay_parameter

		def gamma_function_delay():
			#mass_factor = 0.003,a_parameter = 2, loc = 0, scale = 3
			'''
			the gamma function for a_parameter = 2 and loc = 0 produces a peak at scale so we have a two parameter sfr.
			
			Later we can also release a to have a larger diversity in functional form.
			'''
			from scipy.stats import gamma
			#norm = sum(self.sfr)*number_factor
			#self.infall = np.divide(self.infall*norm,sum(self.infall))
			
			full_time = np.linspace(0, end_of_time, (end_of_time/self.dt)+1)
			feedback_number = gamma.pdf(full_time,a_parameter,loc,scale)

			#for i in range(len(full_time)):
			#	if full_time[i]>=tau_8:
			#		feedback_number[i] = N_0 * np.power(np.divide(full_time[i],tau_8),-1*s_exponent) * np.divide(s_exponent-1,tau_8)		
			feedback_number = np.divide(feedback_number*number_factor,sum(feedback_number))
			
			#number_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dn,self.x,self.sn1a_min,self.sn1a_max)
			#mass_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dm,self.x,self.sn1a_min,self.sn1a_max)
			# for sn1a max and min == 8 and 2 we get:
			self.mean_mass_of_feedback = float(-self.sn1a_yields[0.02]['mass_in_remnants']) ## This mass is turned into the explosion
			#=1.37
			self.mean_mass = 3.21
			#=3.21
			self.mean_remnant_mass = 1.19
			#=1.19
			self.mean_accretion_mass = self.mean_mass_of_feedback - self.mean_remnant_mass ## This comes from Hydrogen feedback
			#=0.18
			feedback_mass = feedback_number * self.mean_mass_of_feedback
			feedback_mass = np.interp(self.t,full_time,feedback_mass)
			feedback_number = np.interp(self.t,full_time,feedback_number)
			self.sn1a_feedback_mass = feedback_mass
			self.sn1a_feedback_number = feedback_number

		def maoz_timedelay():
			'''
			Calculating the delay time distribution of the SNIa explosion. Stochastic sampling is possible if wanted.
			'''
			#number_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dn,self.x,self.sn1a_min,self.sn1a_max) ##Analytic result for Chabrier IMF and sn1a min max 1,8: 0.182189794774
			#mass_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dm,self.x,self.sn1a_min,self.sn1a_max) ##Analytic result for Chabrier IMF and sn1a min max 1,8: 0.398074766434

			full_time = list(self.t)
			while full_time[-1] < end_of_time:
				full_time.append(full_time[-1]+self.dt)
			full_time = np.array(full_time)
			feedback_number = np.zeros_like(full_time)
			for i in range(len(full_time)):
				if full_time[i]>=tau_8:
					feedback_number[i] = np.power(np.divide(full_time[i],tau_8),-1*s_exponent) * np.divide(s_exponent-1,tau_8)# * N_0 * number_of_stars_in_mass_range_for_remnant
			feedback_number = np.divide(feedback_number,sum(feedback_number)) * N_0# * number_of_stars_in_mass_range_for_remnant
			#N_0 now is the number of SN1a per 1Msun
			if stochastic_IMF:
				number_of_potential_sn1a_explosions = int(round(ssp_mass))#int(round(number_of_stars_in_mass_range_for_remnant * ssp_mass) )
				random_sample = np.random.uniform(size = number_of_potential_sn1a_explosions)
				number_of_explosions = len(random_sample[np.where(random_sample<N_0)])
				random_number = np.random.uniform(low = 0.0, high = sum(feedback_number), size = number_of_explosions)
				counting = np.cumsum(feedback_number)
				for i in range(len(counting)-1):
					if i == 0:
						cut = np.where(np.logical_and(random_number>0.,random_number<=counting[i]))
					else:
						cut = np.where(np.logical_and(random_number>counting[i-1],random_number<=counting[i]))
					number_of_stars_in_time_bin = len(random_number[cut])
					feedback_number[i] = number_of_stars_in_time_bin
				feedback_number /= ssp_mass # This way the number of sn1a exploding is NOT fixed by the analytic solution

			# for sn1a max and min == 8 and 2 we get:
			self.mean_mass_of_feedback = float(-self.sn1a_yields[0.02]['mass_in_remnants']) ## This mass is turned into the explosion
			#=1.37
			#if number_of_stars_in_mass_range_for_remnant != 0:
			self.mean_mass = 2.2#mass_of_stars_in_mass_range_for_remnant/number_of_stars_in_mass_range_for_remnant
			#else:
			#	self.mean_mass = 0.
			#=3.21
			self.mean_remnant_mass = self.mean_mass * 0.3 ### approximately the remnant mass of the star exploding in a sn1a. This will be subtracted from the remnant masses
			#=1.19
			self.mean_accretion_mass = self.mean_mass_of_feedback - self.mean_remnant_mass ## This comes from Hydrogen feedback following the single degenerate prescription
			#=0.18
			feedback_number = feedback_number[:len(self.t)]
			feedback_mass = feedback_number * self.mean_mass_of_feedback
			#feedback_number = np.interp(self.t,full_time,feedback_number)
			self.sn1a_feedback_mass = feedback_mass
			self.sn1a_feedback_number = feedback_number

		def normal_timedelay():
			full_time = np.linspace(0, end_of_time, (end_of_time/self.dt)+1)
			feedback_mass = np.zeros_like(full_time)
			for i in range(len(full_time)):
				feedback_mass[i] = np.exp(np.divide(-(full_time[i]-time_delay_peak),time_delay_time_scale))
			j = 0
			while full_time[j]<time_delay_peak:
				feedback_mass[j] = np.exp(-0.5*((full_time[j]-time_delay_peak)/float(gauss_beginning))**2) #1./(std*np.sqrt(2*np.pi))*
				j += 1
	

			number_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dn,self.x,self.sn1a_min,self.sn1a_max)
			mass_of_stars_in_mass_range_for_remnant = imf_mass_fraction_non_nativ(self.dm,self.x,self.sn1a_min,self.sn1a_max)
			
			self.mean_mass_of_feedback = float(-self.sn1a_yields[0.02]['mass_in_remnants']) ## This mass is turned into the explosion
			self.mean_mass = mass_of_stars_in_mass_range_for_remnant/number_of_stars_in_mass_range_for_remnant
			self.mean_remnant_mass = self.mean_mass * 0.37 # interpolation from Karakas yields
			self.mean_accretion_mass = self.mean_mass_of_feedback - self.mean_remnant_mass ## This comes from Hydrogen feedback
			
			mass_fraction_detonating = self.mean_mass_of_feedback * number_of_stars_in_mass_range_for_remnant * pn_number_detonating_until_12Gyr
			number_fraction_detonating = number_of_stars_in_mass_range_for_remnant * pn_number_detonating_until_12Gyr

			feedback_mass = feedback_mass * (mass_fraction_detonating/sum(feedback_mass))
			feedback_number = feedback_mass * (number_fraction_detonating/sum(feedback_mass))
			feedback_mass = np.interp(self.t,full_time,feedback_mass)
			feedback_number = np.interp(self.t,full_time,feedback_number)
			self.sn1a_feedback_mass = feedback_mass
			self.sn1a_feedback_number = feedback_number

		self.sn1a_functional_form = time_delay_functional_form
		timedelays = {'normal': normal_timedelay, 'maoz': maoz_timedelay, 'gamma_function': gamma_function_delay}
		delay_name = self.sn1a_functional_form
		self.sn1a_elements = sn1a_elements
		self.sn1a_metallicities = np.sort(sn1a_metallicities)
		self.sn1a_min = sn1a_min
		self.sn1a_max = sn1a_max
		self.sn1a_yields = sn1a_yields

		######### Calculating and normalising the delay time function of SN1a
		if delay_name in timedelays:
			timedelays[delay_name]() # + argument list of course
		else:
			raise Exception("Time delay '%s' not implemented" % delay_name)

		self.table['mass_in_remnants'] -= self.sn1a_feedback_mass * (self.mean_remnant_mass/self.mean_mass_of_feedback) ### how much mass will be turned from remnants into feedback (single degenerate scenario)
		self.sn1a_table['mass_in_remnants'] -= self.sn1a_feedback_mass * (self.mean_remnant_mass/self.mean_mass_of_feedback)
		self.table['hydrogen_mass_accreted_onto_white_dwarfs'] = self.sn1a_feedback_mass * (1.-(self.mean_remnant_mass/self.mean_mass_of_feedback))
		self.table['H'] -= self.sn1a_feedback_mass * (1.-(self.mean_remnant_mass/self.mean_mass_of_feedback))
		self.table['sn1a'] = self.sn1a_feedback_number
		self.sn1a_table['number_of_events'] = self.sn1a_feedback_number
		
		######### Interpolation of yieldsets for different metallicities
		metallicity_list = []
		if len(self.sn1a_metallicities) == 1:
			metallicity_list.append(self.sn1a_metallicities[0])
		elif self.z < min(self.sn1a_metallicities):
			metallicity_list.append(min(self.sn1a_metallicities))
		elif self.z > max(self.sn1a_metallicities):
			metallicity_list.append(max(self.sn1a_metallicities))
		elif self.z in self.sn1a_metallicities:
			metallicity_list.append(self.z)
		else:
			j=1
			while self.sn1a_metallicities[j] < self.z:
				j += 1
			metallicity_list.append(self.sn1a_metallicities[j-1])
			metallicity_list.append(self.sn1a_metallicities[j])
		### the loop will be run through 2 times (if metallicity not outside or exactly at one of the precalculated metallicities) and the values of the two tables will be interpolated according to the prescribed function
		tables_to_interpolate = []
		for s,metallicity_key in enumerate(metallicity_list):
			tables_to_interpolate.append(np.zeros_like(self.table))
			for element_index, element_name in enumerate(list(set(self.elements).intersection(self.sn1a_elements))):
				tables_to_interpolate[s][element_name] = self.sn1a_yields[metallicity_key][element_name]
			#tables_to_interpolate[s]['mass_in_remnants'] = self.sn1a_yields[metallicity_key]['mass_in_remnants']
			########## end of loop which is metallicity independent
		
		metallicity_weight = []
		if len(metallicity_list) == 1:
			metallicity_weight.append(1.)
			for element_index, element_name in enumerate(list(set(self.elements).intersection(self.sn1a_elements))):
				self.table[element_name] += self.sn1a_feedback_mass * tables_to_interpolate[0][element_name]
				self.sn1a_table[element_name] += self.sn1a_feedback_mass * tables_to_interpolate[0][element_name]
			 #tables_to_interpolate[0]['mass_in_remnants']
			
		else:
			if self.interpolation_scheme == 'linear':
				distance = metallicity_list[1]-metallicity_list[0]
				metallicity_weight.append((metallicity_list[1] - self.z) / float(distance))
				metallicity_weight.append((self.z - metallicity_list[0]) / float(distance))
			if self.interpolation_scheme == 'logarithmic':
				if metallicity_list[0] == 0: # zero metallicity is problematic
					metallicity_list[0] = 0.000000001
				distance = np.log10(metallicity_list[1]) - np.log10(metallicity_list[0])
				metallicity_weight.append(( np.log10(metallicity_list[1]) - np.log10(self.z)) / float(distance))
				metallicity_weight.append(( np.log10(self.z) - np.log10(metallicity_list[0])) / float(distance))
			for i,item in enumerate(metallicity_list):
				tmp=[]
				for element_index, element_name in enumerate(list(set(self.elements).intersection(self.sn1a_elements))):
					self.table[element_name] += (self.sn1a_feedback_mass * tables_to_interpolate[i][element_name]) * metallicity_weight[i]
					self.sn1a_table[element_name] += (self.sn1a_feedback_mass * tables_to_interpolate[i][element_name]) * metallicity_weight[i]
				#self.table['mass_in_remnants'] -= self.sn1a_feedback_mass/float(len(metallicity_list)) #tables_to_interpolate[i]['mass_in_remnants'] * metallicity_weight[i]
				tmp.append(sum((self.sn1a_feedback_mass * tables_to_interpolate[i][element_name])))

## old declarations, the above ones are improved and should be used. BH and PAGB feedback are dummies for feeding back wind without enrichment

	def sn2_feedback_IRA(self, sn2_elements, sn2_yields, sn2_metallicities, sn2_mmin, sn2_mmax,fractions_in_gas):	
		'''
		The mass fraction of the IMF between sn2_mmin and sn2_mmax is fed back instantaneously
		to the ISM according to the relative yields of sn2. The interpolation is linear in mass and metallicity
		Also the mass transformed into remnants is calculated.
		The routine is sensitive to the ordering of the masses in the yield table, it must begin with the smallest increase to the biggest value.
		'''


		# tracking the elemental feedback of individual processes
		additional_keys = ['kinetic_energy','number_of_events','mass_in_remnants']
		names = additional_keys + self.elements # not sure if all elements should be taken (might be easier to add the 3 tables in order to get total yield)
		base = np.zeros(len(self.t))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		self.sn2_table = np.core.records.fromarrays(list_of_arrays,names=names)
		

		#from pycallgraph import PyCallGraph
		#from pycallgraph.output import GraphvizOutput

		#with PyCallGraph(output=GraphvizOutput()):


		self.sn2_elements = sn2_elements
		self.sn2 = sn2_yields
		self.sn2_mmin = sn2_mmin
		self.sn2_mmax = sn2_mmax
		self.sn2_metallicities = sn2_metallicities


		### interpolation of 2 yield sets with different metallicities
		self.sn2_metallicities = np.sort(self.sn2_metallicities)
		metallicity_list = []
		if len(self.sn2_metallicities) == 1:
			metallicity_list.append(self.sn2_metallicities[0])
		elif self.z < min(self.sn2_metallicities):
			metallicity_list.append(min(self.sn2_metallicities))
		elif self.z > max(self.sn2_metallicities):
			metallicity_list.append(max(self.sn2_metallicities))
		elif self.z in self.sn2_metallicities:
			metallicity_list.append(self.z)
		else:
			j=1
			while self.sn2_metallicities[j] < self.z:
				j += 1
			metallicity_list.append(self.sn2_metallicities[j-1])
			metallicity_list.append(self.sn2_metallicities[j])
		### the loop will be run through 2 times (if metallicity not outside or exactly at one of the precalculated metallicities) and the values of the two tables will be interpolated according to the prescribed function
		tables_to_interpolate = []
		for s,metallicity_key in enumerate(metallicity_list):
			tables_to_interpolate.append(np.zeros_like(self.table[1]))
			################################## loop which is metallicity independent		
			##### yield table is cut down such that only yields for masses between sn2_mmin and sn2_mmax are left in
			self.sn2[metallicity_key] = self.sn2[metallicity_key][np.where(np.logical_and(self.sn2[metallicity_key]['Mass']>=self.sn2_mmin,self.sn2[metallicity_key]['Mass']<=self.sn2_mmax))]

			tmp_masses = self.sn2[metallicity_key]['Mass']
			support = []
			support.append(self.sn2_mmin)
			for i in range(len(tmp_masses)-1): 
				support.append(np.mean(tmp_masses[i:i+2]))
			support.append(self.sn2_mmax)
			support = np.array(support)
			### support spans the mass ranges where the yields of a specific stellar mass (tmp_masses) will be applied

			weights = np.zeros_like(tmp_masses)
			for i in range(len(tmp_masses)):
				cut = np.where(np.logical_and(self.x>=support[i],self.x<support[i+1]))
				weights[i] = sum(self.dm[cut]) 
				# weights hold the weights with which each mass_range of the yield set needs to be multiplied
			
			###### here the feedback of the elements is calculated
			for item in list(set(self.elements).intersection(self.sn2_elements))+['mass_in_remnants']:
				tables_to_interpolate[s][item] = sum(self.sn2[metallicity_key][item]*weights)
			######## here the unprocessed wind feedback is calculated
			if 'unprocessed_mass_in_winds' in self.sn2[metallicity_key].dtype.names:
				for i,item in enumerate(self.elements):
					#print i, item, fractions_in_gas[i],weights
					tables_to_interpolate[s][item] += sum(self.sn2[metallicity_key]['unprocessed_mass_in_winds']*weights*fractions_in_gas[i])

		metallicity_weight = []
		if len(metallicity_list) == 1:
			metallicity_weight.append(1.)
			for element_name in list(set(self.elements).intersection(self.sn2_elements))+['mass_in_remnants']:
				self.table[element_name][1] += tables_to_interpolate[0][element_name]
				self.sn2_table[element_name][1] += tables_to_interpolate[0][element_name]
		else:
			if self.interpolation_scheme == 'linear':
				distance = metallicity_list[1]-metallicity_list[0]
				metallicity_weight.append((metallicity_list[1] - self.z) / float(distance))
				metallicity_weight.append((self.z - metallicity_list[0]) / float(distance))
			if self.interpolation_scheme == 'logarithmic':
				distance = np.log10(metallicity_list[1]) - np.log10(metallicity_list[0])
				metallicity_weight.append(( np.log10(metallicity_list[1]) - np.log10(self.z)) / float(distance))
				metallicity_weight.append(( np.log10(self.z) - np.log10(metallicity_list[0])) / float(distance))
			for i,item in enumerate(metallicity_list):
				for element_name in list(set(self.elements).intersection(self.sn2_elements))+['mass_in_remnants']:
					self.table[element_name][1] += tables_to_interpolate[i][element_name] * metallicity_weight[i]
					self.sn2_table[element_name][1] += tables_to_interpolate[i][element_name] * metallicity_weight[i]
		### here the number of stars going sn2 is calculated
		self.table['sn2'][1] = imf_mass_fraction_non_nativ(self.dn,self.x,sn2_mmin,sn2_mmax)
		self.sn2_table['number_of_events'][1] = imf_mass_fraction_non_nativ(self.dn,self.x,sn2_mmin,sn2_mmax)

	def bh_feedback(self,bhmmin,bhmmax,element_list,fractions_in_gas,percentage_of_bh_mass):
		'''
		Old routine, just no enrichment for a specific mass range
		'''
		cut = np.where(np.logical_and(self.x>=bhmmin,self.x<bhmmax))
		weight = sum(self.dm[cut])

		self.table['mass_in_remnants'][1] += percentage_of_bh_mass*weight
		for i,item in enumerate(element_list):
			self.table[item][1] += (1 - percentage_of_bh_mass)*weight*fractions_in_gas[i]
		self.table['bh'][1] = imf_mass_fraction_non_nativ(self.dn,self.x,bhmmin,bhmmax)	

	def post_agb_feedback(self,mmin,mmax,element_list,fractions_in_gas,percentage_to_remnant):
		'''
		just to produce no new elements for stars between agb and sn2, like in kobayashi 2011
		'''

		cut = np.where(np.logical_and(self.x>=mmin,self.x<mmax))
		weight = sum(self.dm[cut])
		if len(self.table) >= 3:
			self.table['mass_in_remnants'][2] += percentage_to_remnant*weight
			for i,item in enumerate(element_list):
				self.table[item][2] += (1 - percentage_to_remnant)*weight*fractions_in_gas[i]
			self.table['pn'][2] += imf_mass_fraction_non_nativ(self.dn,self.x,mmin,mmax)	

