import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline 
import os,os.path
import re
from numpy.lib.recfunctions import append_fields
from . import localpath

class SN1a_feedback(object):
	def __init__(self):    
		"""
		this is the object that holds the feedback table for SN1a
		
		   .masses gives a list of masses
		
		   .metallicities gives a list of possible yield metallicities
		
		   .elements gives the elements considered in the yield table
		
		   .table gives a dictionary where the yield table for a specific metallicity can be queried
		
		   .table[0.02] gives a yield table.
		
		      Keys of this object are ['Mass','mass_in_remnants','elements']
		
		         Mass is in units of Msun
		
		         'mass_in_remnants' in units of Msun but with a '-'
		
		         'elements' yield in Msun normalised to Mass. i.e. integral over all elements is unity 
		"""
	def Seitenzahl(self):
		"""
		Seitenzahl 2013 from Ivo txt
		"""
		y = np.genfromtxt(localpath + 'input/yields/Seitenzahl2013/0.02.txt', names = True, dtype = None)
		self.metallicities = list([0.02])
		self.masses = list([1.4004633930489443])
		names = list(y.dtype.names)
		self.elements = names[2:]
		base = np.zeros(len(self.masses))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
		for name in names:
			if name in ['Mass','mass_in_remnants']:
				yield_tables_final_structure_subtable[name] = y[name]
			else:
				yield_tables_final_structure_subtable[name] = np.divide(y[name],self.masses)
		yield_tables_final_structure = {}
		yield_tables_final_structure[0.02] = yield_tables_final_structure_subtable
		self.table = yield_tables_final_structure
	def Thielemann(self):
		"""
		Thilemann 2003 yields as compiled in Travaglio 2004
		"""	
		y = np.genfromtxt(localpath + 'input/yields/Thielemann2003/0.02.txt', names = True, dtype = None)

		metallicity_list = [0.02]
		self.metallicities = metallicity_list
		self.masses = [1.37409]
		names = y.dtype.names
		base = np.zeros(len(self.masses))
		list_of_arrays = []
		for i in range(len(names)):
			list_of_arrays.append(base)
		yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
		for name in names:
			if name in ['Mass','mass_in_remnants']:
				yield_tables_final_structure_subtable[name] = y[name]
			else:
				yield_tables_final_structure_subtable[name] = np.divide(y[name],self.masses)
		self.elements = list(y.dtype.names[2:])
		yield_tables_final_structure = {}
		yield_tables_final_structure[0.02] = yield_tables_final_structure_subtable
		self.table = yield_tables_final_structure

	def Iwamoto(self):
		'''
		Iwamoto99 yields building up on Nomoto84
		'''
		import numpy.lib.recfunctions as rcfuncs

		tdtype =   [('species1','|S4'),('W7',float),('W70',float),('WDD1',float),('WDD2',float),('WDD3',float),('CDD1',float),('CDD2',float)]
		metallicity_list = [0.02,0.0]
		self.metallicities = metallicity_list
		self.masses = [1.38]
		y = np.genfromtxt(localpath + 'input/yields/Iwamoto/sn1a_yields.txt',dtype = tdtype, names = None)
		## Python3 need transformation between bytes and strings
		element_list2 = []
		for j,jtem in enumerate(y['species1']):
			element_list2.append(jtem.decode('utf8'))
		y = rcfuncs.append_fields(y,'species',element_list2,usemask = False)


		################################
		without_radioactive_isotopes=True
		if without_radioactive_isotopes:### without radioactive isotopes it should be used this way because the radioactive nuclides are already calculated in here
			carbon_list = ['12C','13C']
			nitrogen_list = ['14N','15N']
			oxygen_list = ['16O','17O','18O']
			fluorin_list = ['19F']
			neon_list = ['20Ne','21Ne','22Ne']#,'22Na']
			sodium_list = ['23Na']
			magnesium_list = ['24Mg','25Mg','26Mg']#,'26Al']
			aluminium_list = ['27Al']
			silicon_list = ['28Si','29Si','30Si']
			phosphorus_list = ['31P']
			sulfur_list = ['32S','33S','34S','36S']
			chlorine_list = ['35Cl','37Cl']
			argon_list = ['36Ar','38Ar','40Ar']#, '36Cl']
			potassium_list = ['39K','41K']#, '39Ar', '41Ca']
			calcium_list = ['40Ca','42Ca','43Ca','44Ca','46Ca','48Ca']#, '40K']
			scandium_list = ['45Sc']#,'44Ti']
			titanium_list = ['46Ti','47Ti','48Ti','49Ti','50Ti']#,'48V','49V']
			vanadium_list = ['50V','51V']
			chromium_list = ['50Cr','52Cr','53Cr','54Cr']#,'53Mn']
			manganese_list = ['55Mn']
			iron_list = ['54Fe', '56Fe','57Fe','58Fe']#,'56Co','57Co']
			cobalt_list = ['59Co']#,'60Fe','56Ni','57Ni','59Ni']
			nickel_list = ['58Ni','60Ni','61Ni','62Ni','64Ni']#,'60Co']
			copper_list = ['63Cu','65Cu']#,'63Ni']
			zinc_list = ['64Zn','66Zn','67Zn','68Zn']
		##### with radioactive isotopes (unclear weather they are double, probably not but remnant mass is too big)
		else:
			carbon_list = ['12C','13C']
			nitrogen_list = ['14N','15N']
			oxygen_list = ['16O','17O','18O']
			fluorin_list = ['19F']
			neon_list = ['20Ne','21Ne','22Ne','22Na']
			sodium_list = ['23Na']
			magnesium_list = ['24Mg','25Mg','26Mg','26Al']
			aluminium_list = ['27Al']
			silicon_list = ['28Si','29Si','30Si']
			phosphorus_list = ['31P']
			sulfur_list = ['32S','33S','34S','36S']
			chlorine_list = ['35Cl','37Cl']
			argon_list = ['36Ar','38Ar','40Ar', '36Cl']
			potassium_list = ['39K','41K', '39Ar', '41Ca']
			calcium_list = ['40Ca','42Ca','43Ca','44Ca','46Ca','48Ca', '40K']
			scandium_list = ['45Sc','44Ti']
			titanium_list = ['46Ti','47Ti','48Ti','49Ti','50Ti','48V','49V']
			vanadium_list = ['50V','51V']
			chromium_list = ['50Cr','52Cr','53Cr','54Cr','53Mn']
			manganese_list = ['55Mn']
			iron_list = ['54Fe', '56Fe','57Fe','58Fe','56Co','57Co','56Ni','57Ni']
			cobalt_list = ['59Co','60Fe','59Ni']
			nickel_list = ['58Ni','60Ni','61Ni','62Ni','64Ni','60Co']
			copper_list = ['63Cu','65Cu','63Ni']
			zinc_list = ['64Zn','66Zn','67Zn','68Zn']


		indexing = {}
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Ar'] = argon_list
		indexing['K'] = potassium_list
		indexing['Ca'] = calcium_list
		indexing['Sc'] = scandium_list
		indexing['Ti'] = titanium_list
		indexing['V'] = vanadium_list
		indexing['Cr'] = chromium_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list
		indexing['Cu'] = copper_list
		indexing['Zn'] = zinc_list
		
		self.elements = list(indexing.keys())        
		

		#################################
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(metallicity_list[:]):
			if metallicity == 0.02:
				model = 'W7'
			elif metallicity == 0.0:
				model = 'W70'
			else:
				print('this metallicity is not represented in the Iwamoto yields. They only have solar (0.02) and zero (0.0001)')
			additional_keys = ['Mass', 'mass_in_remnants']
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = self.masses[0]
			
			total_mass = []
			for i,item in enumerate(self.elements):
				for j,jtem in enumerate(indexing[item]):
					cut = np.where(y['species']==jtem)
					yield_tables_final_structure_subtable[item] += y[model][cut]
					total_mass.append(y[model][cut])
			yield_tables_final_structure_subtable['mass_in_remnants'] = -sum(total_mass)
			for i,item in enumerate(self.elements):
				yield_tables_final_structure_subtable[item] = np.divide(yield_tables_final_structure_subtable[item],-yield_tables_final_structure_subtable['mass_in_remnants'])
			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable

		self.table = yield_tables_final_structure

class SN2_feedback(object):
	def __init__(self):
		"""
		This is the object that holds the feedback table for CC-SN.
                Different tables can be loaded by the methods.
		"""
	def Portinari(self):
		'''
		Loading the yield table from Portinari1998.
		'''
		self.metallicities = [0.0004,0.004,0.008,0.02,0.05]
		x = np.genfromtxt(localpath + 'input/yields/Portinari_1998/0.02.txt',names=True)
		self.masses = list(x['Mass'])
		self.elements = list(x.dtype.names[3:])
		
		yield_tables_final_structure = {}
		for metallicity in self.metallicities:
			additional_keys = ['Mass', 'mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = np.array(self.masses)
			
			x = np.genfromtxt(localpath + 'input/yields/Portinari_1998/%s.txt' %(metallicity),names=True)	
			for item in self.elements:
				yield_tables_final_structure_subtable[item] = np.divide(x[item],x['Mass'])
			yield_tables_final_structure_subtable['mass_in_remnants'] = np.divide(x['Mass'] - x['ejected_mass'], x['Mass'])
			
			for i,item in enumerate(self.masses):
				yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][i] = 1. - (yield_tables_final_structure_subtable['mass_in_remnants'][i] + sum(list(yield_tables_final_structure_subtable[self.elements][i])))

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable

		self.table = yield_tables_final_structure



	def francois(self):
		'''
		Loading the yield table of Francois et. al. 2004. Taken from the paper table 1 and 2 and added O H He from WW95 table 5A and 5B
		where all elements are for Z=Zsun and values for Msun > 40 have been stayed the same as for Msun=40.
		Values from 11-25 Msun used case A from WW95 and 30-40 Msun used case B.
		'''
		y = np.genfromtxt(localpath + 'input/yields/Francois04/francois_yields.txt',names=True)
		self.elements = list(y.dtype.names[1:])
		self.masses = y[y.dtype.names[0]]
		self.metallicities = [0.02]
		######### going from absolute ejected masses to relative ejected masses normed with the weight of the initial star
		for i,item in enumerate(y.dtype.names[1:]):
			y[item] = np.divide(y[item],y['Mass'])
		yield_tables = {}
		for i,item in enumerate(self.metallicities):
			yield_tables[item] = y 
		self.table = yield_tables

	def chieffi04(self):

		'''
		Loading the yield table of chieffi04.
		'''
		DATADIR = localpath + 'input/yields/Chieffi04'
		if not os.path.exists(DATADIR):
			os.mkdir(DATADIR)

		MASTERFILE = '{}/chieffi04_yields'.format(DATADIR)

		def _download_chieffi04():
			"""
			Downloads chieffi 04 yields from Vizier.
			"""
			url = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J%2FApJ%2F608%2F405'
			import urllib
			print('Downloading Chieffi 04 yield tables from Vizier (should happen only at the first time)...')
			if os.path.exists(MASTERFILE):
				os.remove(MASTERFILE)
			urllib.urlretrieve(url,MASTERFILE)

			import tarfile
			tar = tarfile.open(MASTERFILE)
			tar.extractall(path=DATADIR)
			tar.close()

		if not os.path.exists(MASTERFILE):
			_download_chieffi04()

		tdtype =   [('metallicity',float),('date_after_explosion',float),('species','|S5'),('13',float),('15',float),('20',float),('25',float),('30',float),('35',float)]
		


		y = np.genfromtxt('%s/yields.dat' %(DATADIR), dtype = tdtype, names = None)
		metallicity_list = np.unique(y['metallicity'])
		self.metallicities = np.sort(metallicity_list)
		number_of_species = int(len(y)/len(self.metallicities))
		tables = []
		for i, item in enumerate(self.metallicities):
			tables.append(y[(i*number_of_species):((i+1)*number_of_species)])
		
		#############################################
		for i in range(len(tables)):
			tables[i] = tables[i][np.where(tables[i]['date_after_explosion']==0)]
		element_list = tables[0]['species'][3:]
		# For python 3 the bytes need to be changed into strings
		element_list2 = []
		for i, item in enumerate(element_list):
			element_list2.append(item.decode('utf8'))
		element_list = np.array(element_list2)
		indexing = [re.split(r'(\d+)', s)[1:] for s in element_list]
		element_position = []
		for i,item in enumerate(element_list):
			element_position.append(indexing[i][1])
		self.elements = list(np.unique(element_position))
		masses = tables[0].dtype.names[3:]
		masses_list = []
		for i,item in enumerate(masses):
			masses_list.append(int(item))
		self.masses = masses_list

		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = tables[metallicity_index]
			additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			
			yield_tables_final_structure_subtable['Mass'] = np.array(self.masses)
			for j,jtem in enumerate(self.masses):
				yield_tables_final_structure_subtable['mass_in_remnants'][j] = yields_for_one_metallicity[str(jtem)][1] / float(jtem) # ,yield_tables_final_structure_subtable['Mass'][i])
				for i,item in enumerate(self.elements):
					################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
					for t,ttem in enumerate(element_position):
						if ttem == item:
							yield_tables_final_structure_subtable[item][j] += yields_for_one_metallicity[str(jtem)][t+3] / float(jtem)
			# remnant + yields of all elements is less than the total mass. In the next loop the wind mass is calculated.
			name_list = list(yield_tables_final_structure_subtable.dtype.names[3:]) + ['mass_in_remnants']
			for i in range(len(yield_tables_final_structure_subtable)):
				tmp = []
				for j,jtem in enumerate(name_list):
					tmp.append(yield_tables_final_structure_subtable[jtem][i])
				tmp = sum(tmp)
				yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][i] = 1 - tmp
			
			yield_tables_final_structure[self.metallicities[metallicity_index]] = yield_tables_final_structure_subtable#[::-1]
		self.table = yield_tables_final_structure


	def chieffi04_net(self):

		'''
		Loading the yield table of chieffi04 corrected for Anders & Grevesse 1989 solar scaled initial yields
		'''
		DATADIR = localpath + 'input/yields/Chieffi04'
		if not os.path.exists(DATADIR):
			os.mkdir(DATADIR)

		MASTERFILE = '{}/chieffi04_yields'.format(DATADIR)

		def _download_chieffi04():
			"""
			Downloads chieffi 04 yields from Vizier.
			"""
			url = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J%2FApJ%2F608%2F405'
			import urllib
			print('Downloading Chieffi 04 yield tables from Vizier (should happen only at the first time)...')
			if os.path.exists(MASTERFILE):
				os.remove(MASTERFILE)
			urllib.urlretrieve(url,MASTERFILE)

			import tarfile
			tar = tarfile.open(MASTERFILE)
			tar.extractall(path=DATADIR)
			tar.close()

		if not os.path.exists(MASTERFILE):
			_download_chieffi04()

		tdtype =   [('metallicity',float),('date_after_explosion',float),('species','|S5'),('13',float),('15',float),('20',float),('25',float),('30',float),('35',float)]
		


		y = np.genfromtxt('%s/yields.dat' %(DATADIR), dtype = tdtype, names = None)
		metallicity_list = np.unique(y['metallicity'])
		self.metallicities = np.sort(metallicity_list)
		number_of_species = int(len(y)/len(self.metallicities))
		tables = []
		for i, item in enumerate(self.metallicities):
			tables.append(y[(i*number_of_species):((i+1)*number_of_species)])
		
		#############################################
		for i in range(len(tables)):
			tables[i] = tables[i][np.where(tables[i]['date_after_explosion']==0)]
		element_list = tables[0]['species'][3:]
		# For python 3 the bytes need to be changed into strings
		element_list2 = []
		for i, item in enumerate(element_list):
			element_list2.append(item.decode('utf8'))
		element_list = np.array(element_list2)
		indexing = [re.split(r'(\d+)', s)[1:] for s in element_list]
		element_position = []
		for i,item in enumerate(element_list):
			element_position.append(indexing[i][1])
		self.elements = list(np.unique(element_position))
		masses = tables[0].dtype.names[3:]
		masses_list = []
		for i,item in enumerate(masses):
			masses_list.append(int(item))
		self.masses = masses_list

		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yield_tables_final_structure[self.metallicities[metallicity_index]] = np.load(DATADIR + '/chieffi_net_met_ind_%d.npy' %(metallicity_index))
		self.table = yield_tables_final_structure

		#############################################
	def Nugrid(self):
		'''
		loading the Nugrid sn2 stellar yields NuGrid stellar data set. I. Stellar yields from H to Bi for stars with metallicities Z = 0.02 and Z = 0.01
		The wind yields need to be added to the *exp* explosion yields.
		No r-process contribution but s and p process from AGB and massive stars
		delayed and rapid SN Explosiom postprocessing is included. Rapid is not consistent with very massive stars so we use the 'delayed' yield set
		mass in remnants not totally consistent with paper table: [ 6.47634087,  2.67590435,  1.98070676] vs. [6.05,2.73,1.61] see table 4
		same with z=0.02 but other elements are implemented in the right way:[ 3.27070753,  8.99349996,  6.12286813,  3.1179861 ,  1.96401573] vs. [3,8.75,5.71,2.7,1.6]
		we have a switch to change between the two different methods (rapid/delay explosion)
		'''
		import numpy.lib.recfunctions as rcfuncs

		tdtype =   [('empty',int),('element1','|S3'),('165',float),('200',float),('300',float),('500',float),('1500',float),('2000',float),('2500',float)]
		tdtype2 =   [('empty',int),('element1','|S3'),('165',float),('200',float),('300',float),('500',float),('1500',float),('2000',float),('2500',float),('3200',float),('6000',float)]
		
		expdtype =   [('empty',int),('element1','|S3'),('15_delay',float),('15_rapid',float),('20_delay',float),('20_rapid',float),('25_delay',float),('25_rapid',float)]
		expdtype2 =   [('empty',int),('element1','|S3'),('15_delay',float),('15_rapid',float),('20_delay',float),('20_rapid',float),('25_delay',float),('32_delay',float),('32_rapid',float),('60_delay',float)]
		
		yield_tables = {}
		self.metallicities = [0.02,0.01]

		which_sn_model_to_use = 'delay' # 'rapid'

		for i,metallicity_index in enumerate([2,1]): 
			if i == 0:
				z = np.genfromtxt(localpath + 'input/yields/NuGrid_AGB_SNII_2013/set1p%d/element_table_set1.%d_yields_winds.txt' %(metallicity_index,metallicity_index),dtype = tdtype2,names = None,skip_header = 3, delimiter = '&', autostrip = True)
				y = np.genfromtxt(localpath + 'input/yields/NuGrid_AGB_SNII_2013/set1p%d/element_table_set1.%d_yields_exp.txt' %(metallicity_index,metallicity_index),dtype = expdtype2,names = None,skip_header = 3, delimiter = '&', autostrip = True)
				y['15_%s' %(which_sn_model_to_use)] += z['1500']
				y['20_%s' %(which_sn_model_to_use)] += z['2000']
				y['25_delay'] += z['2500']
				y['32_%s' %(which_sn_model_to_use)] += z['3200']
				y['60_delay'] += z['6000']
			else:
				z = np.genfromtxt(localpath +'input/yields/NuGrid_AGB_SNII_2013/set1p%d/element_table_set1.%d_yields_winds.txt' %(metallicity_index,metallicity_index),dtype = tdtype,names = None,skip_header = 3, delimiter = '&', autostrip = True)
				y = np.genfromtxt(localpath + 'input/yields/NuGrid_AGB_SNII_2013/set1p%d/element_table_set1.%d_yields_exp.txt' %(metallicity_index,metallicity_index),dtype = expdtype,names = None,skip_header = 3, delimiter = '&', autostrip = True)
				y['15_%s' %(which_sn_model_to_use)] += z['1500']
				y['20_%s' %(which_sn_model_to_use)] += z['2000']
				y['25_%s' %(which_sn_model_to_use)] += z['2500']
				
			# For python 3 the bytes need to be changed into strings
			element_list2 = []
			for j,item in enumerate(y['element1']):
					element_list2.append(item.decode('utf8'))
			y = rcfuncs.append_fields(y,'element',element_list2,usemask = False)
			
			yield_tables[self.metallicities[i]] = y
		
		self.elements = list(yield_tables[0.02]['element']) 
		# For python 3 the bytes need to be changed into strings
		self.masses = np.array((15,20,25,32,60))

		######
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = yield_tables[metallicity]
			final_mass_name_tag = 'mass_in_remnants'
			additional_keys = ['Mass',final_mass_name_tag]
			names = additional_keys + self.elements
			if metallicity == 0.02:
				base = np.zeros(len(self.masses))
			else:
				base = np.zeros(len(self.masses)-2)
			list_of_arrays = []

			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)

			if metallicity == 0.02:
				yield_tables_final_structure_subtable['Mass'] = self.masses
			else:
				yield_tables_final_structure_subtable['Mass'] = self.masses[:-2]

			for i,item in enumerate(self.elements):
				################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
				if metallicity == 0.02:
					line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['element']==item)]
					temp1 = np.zeros(5)
					temp1[0] = line_of_one_element['15_%s' %(which_sn_model_to_use)]
					temp1[1] = line_of_one_element['20_%s' %(which_sn_model_to_use)]
					temp1[2] = line_of_one_element['25_delay']
					temp1[3] = line_of_one_element['32_%s' %(which_sn_model_to_use)]
					temp1[4] = line_of_one_element['60_delay']
					yield_tables_final_structure_subtable[item] = np.divide(temp1,self.masses)

				else:
					line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['element']==item)]
					temp1 = np.zeros(3)
					temp1[0] = line_of_one_element['15_%s' %(which_sn_model_to_use)]
					temp1[1] = line_of_one_element['20_%s' %(which_sn_model_to_use)]
					temp1[2] = line_of_one_element['25_%s' %(which_sn_model_to_use)]
					yield_tables_final_structure_subtable[item] = np.divide(temp1,self.masses[:-2])

			if metallicity == 0.02:
				yield_tables_final_structure_subtable[final_mass_name_tag][0] = (1-sum(yield_tables_final_structure_subtable[self.elements][0]))
				yield_tables_final_structure_subtable[final_mass_name_tag][1] = (1-sum(yield_tables_final_structure_subtable[self.elements][1]))
				yield_tables_final_structure_subtable[final_mass_name_tag][2] = (1-sum(yield_tables_final_structure_subtable[self.elements][2]))
				yield_tables_final_structure_subtable[final_mass_name_tag][3] = (1-sum(yield_tables_final_structure_subtable[self.elements][3]))
				yield_tables_final_structure_subtable[final_mass_name_tag][4] = (1-sum(yield_tables_final_structure_subtable[self.elements][4]))

			else:        
				yield_tables_final_structure_subtable[final_mass_name_tag][0] = (1-sum(yield_tables_final_structure_subtable[self.elements][0]))
				yield_tables_final_structure_subtable[final_mass_name_tag][1] = (1-sum(yield_tables_final_structure_subtable[self.elements][1]))
				yield_tables_final_structure_subtable[final_mass_name_tag][2] = (1-sum(yield_tables_final_structure_subtable[self.elements][2]))

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable#[::-1]
		self.table = yield_tables_final_structure


	def one_parameter(self, elements, element_fractions):
		"""
		This function was introduced in order to find best-fit yield sets where each element has just a single yield (no metallicity or mass dependence).
                One potential problem is that sn2 feedback has a large fraction of Neon ~ 0.01, the next one missing is Argon but that only has 0.05%. This might spoil the metallicity derivation a bit.
		Another problem: He and the remnant mass fraction is not constrained in the APOGEE data. Maybe these can be constrained externally by yield sets or cosmic abundance standard or solar abundances.
		"""
		self.metallicities = [0.01]
		self.masses = np.array([10])
		self.elements = elements 

		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		

		additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
		names = additional_keys + self.elements
		base = np.zeros(len(self.masses))
		list_of_arrays = []
		
		for i in range(len(names)):
			list_of_arrays.append(base)
		yield_table = np.core.records.fromarrays(list_of_arrays,names=names)
		yield_table['Mass'] = self.masses
		yield_table['mass_in_remnants'] = 0.1
		yield_table['unprocessed_mass_in_winds'] = 1 - yield_table['mass_in_remnants']
		for i,item in enumerate(self.elements[1:]):
			yield_table[item] = element_fractions[i+1]
		yield_table['H'] = -sum(element_fractions[1:])


		yield_tables_final_structure[self.metallicities[0]] = yield_table
		self.table = yield_tables_final_structure


	def Nomoto2013(self):
		'''
		Nomoto2013 sn2 yields from 13Msun onwards
		'''
		import numpy.lib.recfunctions as rcfuncs

		dt = np.dtype('a13,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
		yield_tables = {}
		self.metallicities = [0.0500,0.0200,0.0080,0.0040,0.0010]
		self.masses = np.array((13,15,18,20,25,30,40))
		z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=0.0200.dat',dtype=dt,names = True)
		
		yield_tables_dict = {}
		for item in self.metallicities:
			z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=%.4f.dat' %(item),dtype=dt,names = True)
			yield_tables_dict[item]=z

		hydrogen_list = ['H__1','H__2']
		helium_list = ['He_3','He_4']
		lithium_list = ['Li_6','Li_7']
		berillium_list = ['Be_9']
		boron_list = ['B_10','B_11']
		carbon_list = ['C_12','C_13']
		nitrogen_list = ['N_14','N_15']
		oxygen_list = ['O_16','O_17','O_18']
		fluorin_list = ['F_19']
		neon_list = ['Ne20','Ne21','Ne22']
		sodium_list = ['Na23']
		magnesium_list = ['Mg24','Mg25','Mg26']
		aluminium_list = ['Al27']
		silicon_list = ['Si28','Si29','Si30']
		phosphorus_list = ['P_31']
		sulfur_list = ['S_32','S_33','S_34','S_36']
		chlorine_list = ['Cl35','Cl37']
		argon_list = ['Ar36','Ar38','Ar40']
		potassium_list = ['K_39','K_41']
		calcium_list = ['K_40','Ca40','Ca42','Ca43','Ca44','Ca46','Ca48']
		scandium_list = ['Sc45']
		titanium_list = ['Ti46','Ti47','Ti48','Ti49','Ti50']
		vanadium_list = ['V_50','V_51']
		chromium_list = ['Cr50','Cr52','Cr53','Cr54']
		manganese_list = ['Mn55']
		iron_list = ['Fe54', 'Fe56','Fe57','Fe58']
		cobalt_list = ['Co59']
		nickel_list = ['Ni58','Ni60','Ni61','Ni62','Ni64']
		copper_list = ['Cu63','Cu65']
		zinc_list = ['Zn64','Zn66','Zn67','Zn68','Zn70']
		gallium_list = ['Ga69','Ga71']
		germanium_list = ['Ge70','Ge72','Ge73','Ge74']

		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		indexing['Be'] = berillium_list
		indexing['B'] = boron_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Ar'] = argon_list
		indexing['K'] = potassium_list
		indexing['Ca'] = calcium_list
		indexing['Sc'] = scandium_list
		indexing['Ti'] = titanium_list
		indexing['V'] = vanadium_list
		indexing['Cr'] = chromium_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list
		indexing['Cu'] = copper_list
		indexing['Zn'] = zinc_list
		indexing['Ga'] = gallium_list
		indexing['Ge'] = germanium_list

		self.elements = list(indexing.keys())
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = yield_tables_dict[metallicity]
			# For python 3 the bytes need to be changed into strings
			element_list2 = []
			for j,item in enumerate(yields_for_one_metallicity['M']):
					element_list2.append(item.decode('utf8'))
			yields_for_one_metallicity = rcfuncs.append_fields(yields_for_one_metallicity,'element',element_list2,usemask = False)
			additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = self.masses
			#yield_tables_final_structure_subtable['mass_in_remnants'] = yields_for_one_metallicity['M']
			temp1 = np.zeros(len(self.masses))
			temp1[0] = yields_for_one_metallicity[0][21]
			temp1[1] = yields_for_one_metallicity[0][22]
			temp1[2] = yields_for_one_metallicity[0][23]
			temp1[3] = yields_for_one_metallicity[0][24]
			temp1[4] = yields_for_one_metallicity[0][25]
			temp1[5] = yields_for_one_metallicity[0][26]
			temp1[6] = yields_for_one_metallicity[0][27]

			yield_tables_final_structure_subtable['mass_in_remnants'] = np.divide(temp1,self.masses)
			for i,item in enumerate(self.elements):
				yield_tables_final_structure_subtable[item] = 0
				for j,jtem in enumerate(indexing[item]):
						################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
						line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['element']==jtem)][0]
						temp1 = np.zeros(len(self.masses))
						temp1[0] = line_of_one_element[21]
						temp1[1] = line_of_one_element[22]
						temp1[2] = line_of_one_element[23]
						temp1[3] = line_of_one_element[24]
						temp1[4] = line_of_one_element[25]
						temp1[5] = line_of_one_element[26]
						temp1[6] = line_of_one_element[27]
						yield_tables_final_structure_subtable[item] += np.divide(temp1,self.masses)
						
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][0] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][0]-sum(yield_tables_final_structure_subtable[self.elements][0]))#yields_for_one_metallicity[0][21]#
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][1] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][1]-sum(yield_tables_final_structure_subtable[self.elements][1]))#yields_for_one_metallicity[0][22]#
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][2] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][2]-sum(yield_tables_final_structure_subtable[self.elements][2]))#yields_for_one_metallicity[0][23]#divided by mass because 'mass in remnant' is also normalised
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][3] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][3]-sum(yield_tables_final_structure_subtable[self.elements][3]))#yields_for_one_metallicity[0][24]#
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][4] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][4]-sum(yield_tables_final_structure_subtable[self.elements][4]))#yields_for_one_metallicity[0][25]#
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][5] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][5]-sum(yield_tables_final_structure_subtable[self.elements][5]))#yields_for_one_metallicity[0][26]#
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][6] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][6]-sum(yield_tables_final_structure_subtable[self.elements][6]))#yields_for_one_metallicity[0][27]#
	

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable#[::-1]
		self.table = yield_tables_final_structure


	def Nomoto2013_net(self):
		'''
		Nomoto2013 sn2 yields from 13Msun onwards
		'''
		import numpy.lib.recfunctions as rcfuncs

		dt = np.dtype('a13,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
		yield_tables = {}
		self.metallicities = [0.0500,0.0200,0.0080,0.0040,0.0010]
		self.masses = np.array((13,15,18,20,25,30,40))
		z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=0.0200.dat',dtype=dt,names = True)
		
		yield_tables_dict = {}
		for item in self.metallicities:
			z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=%.4f.dat' %(item),dtype=dt,names = True)
			yield_tables_dict[item]=z

		hydrogen_list = ['H__1','H__2']
		helium_list = ['He_3','He_4']
		lithium_list = ['Li_6','Li_7']
		berillium_list = ['Be_9']
		boron_list = ['B_10','B_11']
		carbon_list = ['C_12','C_13']
		nitrogen_list = ['N_14','N_15']
		oxygen_list = ['O_16','O_17','O_18']
		fluorin_list = ['F_19']
		neon_list = ['Ne20','Ne21','Ne22']
		sodium_list = ['Na23']
		magnesium_list = ['Mg24','Mg25','Mg26']
		aluminium_list = ['Al27']
		silicon_list = ['Si28','Si29','Si30']
		phosphorus_list = ['P_31']
		sulfur_list = ['S_32','S_33','S_34','S_36']
		chlorine_list = ['Cl35','Cl37']
		argon_list = ['Ar36','Ar38','Ar40']
		potassium_list = ['K_39','K_41']
		calcium_list = ['K_40','Ca40','Ca42','Ca43','Ca44','Ca46','Ca48']
		scandium_list = ['Sc45']
		titanium_list = ['Ti46','Ti47','Ti48','Ti49','Ti50']
		vanadium_list = ['V_50','V_51']
		chromium_list = ['Cr50','Cr52','Cr53','Cr54']
		manganese_list = ['Mn55']
		iron_list = ['Fe54', 'Fe56','Fe57','Fe58']
		cobalt_list = ['Co59']
		nickel_list = ['Ni58','Ni60','Ni61','Ni62','Ni64']
		copper_list = ['Cu63','Cu65']
		zinc_list = ['Zn64','Zn66','Zn67','Zn68','Zn70']
		gallium_list = ['Ga69','Ga71']
		germanium_list = ['Ge70','Ge72','Ge73','Ge74']

		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		indexing['Be'] = berillium_list
		indexing['B'] = boron_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Ar'] = argon_list
		indexing['K'] = potassium_list
		indexing['Ca'] = calcium_list
		indexing['Sc'] = scandium_list
		indexing['Ti'] = titanium_list
		indexing['V'] = vanadium_list
		indexing['Cr'] = chromium_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list
		indexing['Cu'] = copper_list
		indexing['Zn'] = zinc_list
		indexing['Ga'] = gallium_list
		indexing['Ge'] = germanium_list

		self.elements = list(indexing.keys())
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):

			yield_tables_final_structure[metallicity] = np.load(localpath + 'input/yields/Nomoto2013/nomoto_net_met_ind_%d.npy' %(metallicity_index))
		self.table = yield_tables_final_structure

#######################
class AGB_feedback(object):
	def __init__(self):   
		"""
		This is the object that holds the feedback table for agb stars.
                The different methods load different tables from the literature. They are in the input/yields/ folder.
		"""
	def Ventura(self):
		"""
		Ventura 2013 net yields from Paolo himself
		"""
		self.metallicities = [0.04,0.018,0.008,0.004,0.001,0.0003]
		x = np.genfromtxt(localpath + 'input/yields/Ventura2013/0.018.txt',names=True)
		self.masses = x['Mass']
		self.elements = ['H', 'He', 'Li','C','N','O','F','Ne','Na','Mg','Al','Si']
		###
		yield_tables_final_structure = {}
		for metallicity in self.metallicities:
			x = np.genfromtxt(localpath + 'input/yields/Ventura2013/%s.txt' %(str(metallicity)),names=True)
			additional_keys = ['Mass', 'mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			base = np.zeros(len(x['Mass']))
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = x['Mass']
			yield_tables_final_structure_subtable['mass_in_remnants'] = np.divide(x['mass_in_remnants'],x['Mass'])	
			for item in self.elements:
				if item == 'C':
					yield_tables_final_structure_subtable[item] = x['C12']
					yield_tables_final_structure_subtable[item] += x['C13']
				elif item == 'N':
					yield_tables_final_structure_subtable[item] = x['N14']
				elif item == 'O':
					yield_tables_final_structure_subtable[item] = x['O16']
					yield_tables_final_structure_subtable[item] += x['O17']
					yield_tables_final_structure_subtable[item] += x['O18']
				elif item == 'F':
					yield_tables_final_structure_subtable[item] = x['F19']
				elif item == 'Ne':
					yield_tables_final_structure_subtable[item] = x['NE20']
					yield_tables_final_structure_subtable[item] += x['NE22']
				elif item == 'Na':
					yield_tables_final_structure_subtable[item] = x['NA23']
				elif item == 'Mg':
					yield_tables_final_structure_subtable[item] = x['MG24']
					yield_tables_final_structure_subtable[item] += x['MG25']
					yield_tables_final_structure_subtable[item] += x['MG26']
				elif item == 'Al':
					yield_tables_final_structure_subtable[item] = x['AL26']
					yield_tables_final_structure_subtable[item] += x['AL27']
				elif item == 'Si':
					yield_tables_final_structure_subtable[item] = x['SI28']
				else:
					yield_tables_final_structure_subtable[item] = x[item]
			for item in self.elements:
				yield_tables_final_structure_subtable[item] = np.divide(yield_tables_final_structure_subtable[item],x['Mass'])

			for i,item in enumerate(x['Mass']):
				yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][i] = 1. - (yield_tables_final_structure_subtable['mass_in_remnants'][i] + sum(list(yield_tables_final_structure_subtable[self.elements][i])))

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable

		self.table = yield_tables_final_structure
		###

	def Nomoto2013(self):
		'''
		Nomoto2013 agb yields up to 6.5Msun and are a copy of Karakas2010. Only that the yields here are given as net yields which does not help so much
		'''
		dt = np.dtype('a13,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
		yield_tables = {}
		self.metallicities = [0.0500,0.0200,0.0080,0.0040,0.0010]
		self.masses = np.array((1.,1.2,1.5,1.8,1.9,2.0,2.2,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0))#,6.5,7.0,8.0,10.))
		z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=0.0200.dat',dtype=dt,names = True)
		
		yield_tables_dict = {}
		for item in self.metallicities:
			z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/nomoto_2013_z=%.4f.dat' %(item),dtype=dt,names = True)
			yield_tables_dict[item]=z
		#########################
		hydrogen_list = ['H__1','H__2']
		helium_list = ['He_3','He_4']
		lithium_list = ['Li_6','Li_7']
		berillium_list = ['Be_9']
		boron_list = ['B_10','B_11']
		carbon_list = ['C_12','C_13']
		nitrogen_list = ['N_14','N_15']
		oxygen_list = ['O_16','O_17','O_18']
		fluorin_list = ['F_19']
		neon_list = ['Ne20','Ne21','Ne22']
		sodium_list = ['Na23']
		magnesium_list = ['Mg24','Mg25','Mg26']
		aluminium_list = ['Al27']
		silicon_list = ['Si28','Si29','Si30']
		phosphorus_list = ['P_31']
		sulfur_list = ['S_32','S_33','S_34','S_36']
		chlorine_list = ['Cl35','Cl37']
		argon_list = ['Ar36','Ar38','Ar40']
		potassium_list = ['K_39','K_41']
		calcium_list = ['K_40','Ca40','Ca42','Ca43','Ca44','Ca46','Ca48']
		scandium_list = ['Sc45']
		titanium_list = ['Ti46','Ti47','Ti48','Ti49','Ti50']
		vanadium_list = ['V_50','V_51']
		chromium_list = ['Cr50','Cr52','Cr53','Cr54']
		manganese_list = ['Mn55']
		iron_list = ['Fe54', 'Fe56','Fe57','Fe58']
		cobalt_list = ['Co59']
		nickel_list = ['Ni58','Ni60','Ni61','Ni62','Ni64']
		copper_list = ['Cu63','Cu65']
		zinc_list = ['Zn64','Zn66','Zn67','Zn68','Zn70']
		gallium_list = ['Ga69','Ga71']
		germanium_list = ['Ge70','Ge72','Ge73','Ge74']

		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		indexing['Be'] = berillium_list
		indexing['B'] = boron_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Ar'] = argon_list
		indexing['K'] = potassium_list
		indexing['Ca'] = calcium_list
		indexing['Sc'] = scandium_list
		indexing['Ti'] = titanium_list
		indexing['V'] = vanadium_list
		indexing['Cr'] = chromium_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list
		indexing['Cu'] = copper_list
		indexing['Zn'] = zinc_list
		indexing['Ga'] = gallium_list
		indexing['Ge'] = germanium_list

		self.elements = indexing.keys()
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = yield_tables_dict[metallicity]
			final_mass_name_tag = 'mass_in_remnants'
			additional_keys = ['Mass',final_mass_name_tag]
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = self.masses
			
			for i,item in enumerate(self.elements):
				yield_tables_final_structure_subtable[item] = 0
				for j,jtem in enumerate(indexing[item]):
						################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
						line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['M']==jtem)][0]
						temp1 = np.zeros(len(self.masses))
						for s in range(len(self.masses)): 
							temp1[s] = line_of_one_element[s+2]
						yield_tables_final_structure_subtable[item] += np.divide(temp1,self.masses)

			for t in range(len(self.masses)):
				yield_tables_final_structure_subtable[final_mass_name_tag][t] = (1-sum(yield_tables_final_structure_subtable[self.elements][t]))#yields_for_one_metallicity[0][21]#

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable#[::-1]
		self.table = yield_tables_final_structure

	def Nugrid(self):
		'''
		loading the Nugrid intermediate mass stellar yields NuGrid stellar data set. I. Stellar yields from H to Bi for stars with metallicities Z = 0.02 and Z = 0.01
		'''
		import numpy.lib.recfunctions as rcfuncs
		tdtype =   [('empty',int),('element1','|S3'),('165',float),('200',float),('300',float),('500',float),('1500',float),('2000',float),('2500',float)]
		yield_tables = {}
		self.metallicities = [0.02,0.01]

		for i,metallicity_index in enumerate([2,1]): 
			y = np.genfromtxt(localpath + 'input/yields/NuGrid_AGB_SNII_2013/set1p%d/element_table_set1.%d_yields_winds.txt' %(metallicity_index,metallicity_index),dtype = tdtype,names = None,skip_header = 3, delimiter = '&', autostrip = True)

			## Python3 need transformation between bytes and strings
			element_list2 = []
			for j,jtem in enumerate(y['element1']):
					element_list2.append(jtem.decode('utf8'))
			y = rcfuncs.append_fields(y,'element',element_list2,usemask = False)
			
			yield_tables[self.metallicities[i]] = y
		
		
		self.elements = list(yield_tables[0.02]['element']) 
		self.masses = np.array((1.65,2.0,3.0,5.0))

		######
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = yield_tables[metallicity]
			final_mass_name_tag = 'mass_in_remnants'
			additional_keys = ['Mass',final_mass_name_tag]
			names = additional_keys + self.elements
			
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)


			yield_tables_final_structure_subtable['Mass'] = self.masses
			for i,item in enumerate(self.elements):
				################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
				line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['element']==item)]
				temp1 = np.zeros(4)
				temp1[0] = line_of_one_element['165']
				temp1[1] = line_of_one_element['200']
				temp1[2] = line_of_one_element['300']
				temp1[3] = line_of_one_element['500']
				yield_tables_final_structure_subtable[item] = np.divide(temp1,self.masses)

			yield_tables_final_structure_subtable[final_mass_name_tag][0] = (1-sum(yield_tables_final_structure_subtable[self.elements][0]))
			yield_tables_final_structure_subtable[final_mass_name_tag][1] = (1-sum(yield_tables_final_structure_subtable[self.elements][1]))
			yield_tables_final_structure_subtable[final_mass_name_tag][2] = (1-sum(yield_tables_final_structure_subtable[self.elements][2]))
			yield_tables_final_structure_subtable[final_mass_name_tag][3] = (1-sum(yield_tables_final_structure_subtable[self.elements][3]))
			

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable[::-1]
		self.table = yield_tables_final_structure
		######

	def Karakas(self):
		'''
		loading the yield table of Karakas 2010.
		'''
		import numpy.lib.recfunctions as rcfuncs

		DATADIR = localpath + 'input/yields/Karakas2010'
		if not os.path.exists(DATADIR):
			os.mkdir(DATADIR)

		MASTERFILE = '{}/karakas_yields'.format(DATADIR)

		def _download_karakas():
			"""
			Downloads Karakas yields from Vizier.
			"""
			#url = 'http://zenodo.org/record/12800/files/dartmouth.h5'
			url = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J%2FMNRAS%2F403%2F1413'
			import urllib
			print('Downloading Karakas 2010 yield tables from Vizier (should happen only at the first time)...')
			if os.path.exists(MASTERFILE):
				os.remove(MASTERFILE)
			urllib.urlretrieve(url,MASTERFILE)

			import tarfile
			tar = tarfile.open(MASTERFILE)
			tar.extractall(path=DATADIR)
			tar.close()

		if not os.path.exists(MASTERFILE):
			_download_karakas()



		tdtype =   [('imass',float),('metallicity',float),('fmass',float),('species1','|S4'),('A',int),('net_yield',float),('ejected_mass',float),('initial_wind',float),('average_wind',float),('initial_mass_fraction',float),('production_factor',float)]
		metallicity_list = [0.02, 0.008, 0.004 ,0.0001]
		self.metallicities = metallicity_list
		


		tables = []
		for i,item in enumerate(metallicity_list):
			y = np.genfromtxt('%s/tablea%d.dat' %(DATADIR,i+2), dtype = tdtype, names = None)
			## Python3 need transformation between bytes and strings
			element_list2 = []
			for j,jtem in enumerate(y['species1']):
					element_list2.append(jtem.decode('utf8'))
			y = rcfuncs.append_fields(y,'species',element_list2,usemask = False)


			tables.append(y)
		

		### easy to extend to other species just make a new list of isotopes (see karakas tables)
		### and then also extend the indexing variable. 
		### The choice for specific elements can be done later when just using specific species
		hydrogen_list = ['n','p','d']
		helium_list = ['he3','he4']
		lithium_list = ['li7','be7','b8']
		
		carbon_list = ['c12','c13','n13']
		nitrogen_list = ['n14','n15','c14','o14','o15']
		oxygen_list = [ 'o16','o17','o18','f17','f18']
		fluorin_list = ['ne19','f19','o19']
		neon_list = ['ne20','ne21','ne22','f20','na21','na22']
		sodium_list = ['na23','ne23','mg23']
		magnesium_list = ['mg24','mg25','mg26','al-6','na24','al25']
		aluminium_list = ['mg27','al*6','al27','si27']
		silicon_list = ['al28','si28','si29','si30','p29','p30']
		phosphorus_list = ['si31','si32','si33','p31']
		sulfur_list = ['s32','s33','s34','p32','p33','p34']
		chlorine_list = ['s35']
		iron_list = ['fe54', 'fe56','fe57','fe58']
		manganese_list = ['fe55']
		cobalt_list = ['ni59','fe59','co59']
		nickel_list = ['ni58','ni60','ni61','ni62','co60','co61','fe60','fe61']


		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list

		#indexing['S_el'] = ni_to_bi
		
		self.elements = list(indexing.keys())
		#### little fix for karakas tablea5.dat: 6.0 M_sun is written two times. We chose the first one
		#tables[3]['imass'][-77:] = 6.5 # this is the fix if the second 6msun line was interpreted as 6.5 msun
		tables[3] = tables[3][:-77]

		#### making the general feedback table with yields for the individual elements
		### loop for the different metallicities
		yield_tables = {}
		for metallicity_index,metallicity in enumerate(metallicity_list[:]):
			### loop for the different elements
			yields_002 = {}
			for i,item1 in enumerate(indexing):
				unique_masses = len(np.unique(tables[metallicity_index]['imass']))
				element = np.zeros((unique_masses,), dtype=[('imass',float),('species','|S4'),('fmass',float),('net_yield',float),('ejected_mass',float),('initial_mass_fraction',float),('initial_wind',float),('average_wind',float),('production_factor',float)])
				for j,item in enumerate(indexing[item1]):
					cut = np.where(tables[metallicity_index]['species']==item)
					temp = tables[metallicity_index][cut]
					if j == 0:
						element['imass'] = temp['imass']
						element['fmass'] = temp['fmass']
						element['species'] = temp['species'] ### just for test purposes
					element['net_yield'] += temp['net_yield']
					element['ejected_mass'] += temp['ejected_mass']
					element['initial_mass_fraction'] += temp['initial_mass_fraction']
					element['initial_wind'] += temp['initial_wind']
					element['average_wind'] += temp['average_wind']
					element['production_factor'] += temp['production_factor']
				yields_002[item1] = element
			yield_tables[metallicity] = yields_002
		
		self.masses = np.unique(tables[0]['imass']) ## table a3 and a4 and maybe a5 are missing 6.5 Msun its probably easier to skip the 6.5 Msun entries altogether for interpolation reasons

		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(metallicity_list[:]):
			yields_for_one_metallicity = yield_tables[metallicity]
			final_mass_name_tag = 'mass_in_remnants'
			additional_keys = ['Mass',final_mass_name_tag]
			names = additional_keys + self.elements
			if metallicity == 0.02: #or metallicity == 0.0001:
				base = np.zeros(len(self.masses))
			else:
				base = np.zeros(len(self.masses)-1)
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)

			yield_tables_final_structure_subtable['Mass'] = yields_for_one_metallicity[self.elements[0]]['imass']
			yield_tables_final_structure_subtable[final_mass_name_tag] = np.divide(yields_for_one_metallicity[self.elements[0]]['fmass'],yield_tables_final_structure_subtable['Mass'])#yields_for_one_metallicity[self.elements[0]]['fmass']
			for i,item in enumerate(self.elements):
				################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
				yield_tables_final_structure_subtable[item] = np.divide(yields_for_one_metallicity[item]['ejected_mass'],yield_tables_final_structure_subtable['Mass'])

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable[::-1]
		self.table = yield_tables_final_structure

	def Karakas16_net(self):
		"""
		load the Karakas 2016 yields send by Amanda and Fishlock 2014 for Z = 0.001. With slight inconsistencies in the mass normalisation and not sure which Asplund2009 solar abundances she uses
		"""
		import numpy.lib.recfunctions as rcfuncs
		import sys


		list_of_metallicities = [0.001,0.007, 0.014, 0.03 ]
		self.metallicities = list_of_metallicities
		data_path = localpath + 'input/yields/Karakas2016/'
		yield_tables = {}
		for metallicity in list_of_metallicities:
			metallicity_name = str(metallicity)[2:]
			if metallicity == 0.001:
				dt = np.dtype([('element1', '|S4'), ('atomic_number', np.int),('yield', np.float),('mass_lost', np.float),('mass_0', np.float),('xi', np.float),('x0', np.float),('log_xi_x0', np.float)])
			else:
				dt = np.dtype([('element1', '|S4'), ('atomic_number', np.int),('log_e', np.float),('xh', np.float),('xfe', np.float),('xi', np.float),('massi', np.float)])
			### yield
			y = np.genfromtxt('%syield_z%s.dat' %(data_path,metallicity_name), dtype=dt)
			
			## Python3 need transformation between bytes and strings
			if sys.version[0] == '3':
				element_list2 = []
				for j,jtem in enumerate(y['element1']):
						element_list2.append(jtem.decode('utf8'))
				y = rcfuncs.append_fields(y,'element',element_list2,usemask = False)
			elif sys.version[0] == '2':
				y = rcfuncs.append_fields(y,'element',y['element1'],usemask = False)
			else:
				print('not a valid python version')


			dt = np.dtype([('element1', '|S4'), ('atomic_number', np.int),('log_e', np.float),('xh', np.float),('xfe', np.float),('xo', np.float),('xi', np.float)])
			### surface
			s = np.genfromtxt('%ssurf_z%s.dat' %(data_path,metallicity_name), dtype=dt)
			## Python3 need transformation between bytes and strings
			if sys.version[0] == '3':
				element_list2 = []
				for j,jtem in enumerate(s['element1']):
						element_list2.append(jtem.decode('utf8'))
				s = rcfuncs.append_fields(s,'element',element_list2,usemask = False)
			elif sys.version[0] == '2':
				s = rcfuncs.append_fields(s,'element',s['element1'],usemask = False)
			else:
				print('not a valid python version')

			t =  np.where(s['element']== 'p')
			len_elements = t[0][2]-1
			elements = list(s['element'][:len_elements])
			for i,item in enumerate(elements):
				if len(elements[i]) == 2:
					elements[i] = str.upper(elements[i][0]) + elements[i][1]
				else:
					elements[i] = str.upper(elements[i][0])
			elements[0] = 'H'
			additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + elements
			base = np.zeros(1)
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			initial_abundances = np.core.records.fromarrays(list_of_arrays,names=names)
			initial_abundances['Mass'] = 1.
			for i,item in enumerate(elements):
				initial_abundances[item] = s['xi'][i]
			### renormalising because the mass fractions add to more than 1
			metals_fraction = sum(list(initial_abundances[0])[5:])
			sum_all = sum(list(initial_abundances[0])[3:])
			for i,item in enumerate(elements):
				initial_abundances[item] /= sum_all
			#### just copied out of the files. Also several masses and other overshootfactors had to be excluded. 
			if metallicity == 0.001:
				list_of_masses = [1.,1.25,1.5,2.0,2.25,2.5,2.75,3.,3.25,3.5,4.,4.5,5.,5.5,6.,7.]
				list_of_remnant = [0.678,0.669,0.657,0.668,0.839,0.948,1.057,1.189,1.403,1.176,1.726,1.659,1.740,1.962,1.725,2.062]
			if metallicity == 0.014:
				list_of_masses = [1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.5,6.,7.,8.]
				list_of_remnant = [0.585,0.605,0.616,0.638,0.66,0.675,0.679,0.684,0.694,0.708,0.73,0.766,0.813,0.853,0.862,0.87,0.879,0.9,0.921,0.976,1.062]
			if metallicity == 0.03:
				list_of_masses = [1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.5,6.,7.,8.]
				list_of_remnant = [0.573,0.590,0.607,0.625,0.643,0.661,0.650,0.670,0.691,0.713,0.727,0.744,0.744,0.806,0.848,0.858,0.867,0.886,0.907,0.963,1.053]
			if metallicity == 0.007:
				list_of_masses = [1.,1.25,1.5,1.75,1.9,2.1,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,5.5,6.,7.,7.5]
				list_of_remnant = [0.606,0.629,0.646,0.641,0.657,0.659,0.663,0.668,0.679,0.698,0.728,0.766,0.802,0.849,0.859,0.873,0.883,0.895,0.921,0.956,1.040,1.116]
			if metallicity == 0.001:
				t = np.where(y['element']=='H')
				len_elements = t[0][1]
				elements = list(y['element'][:len_elements])
			else:
				t =  np.where(y['element']== 'p')
				len_elements = t[0][2]
				elements = list(y['element'][:len_elements])
				for i,item in enumerate(elements):
					if len(elements[i]) == 2:
						elements[i] = str.upper(elements[i][0]) + elements[i][1]
					else:
						elements[i] = str.upper(elements[i][0])
				elements[0] = 'H'
			additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + elements
			base = np.zeros(len(list_of_masses))
			list_of_arrays = []
			for i in range(len(names)):
				list_of_arrays.append(base)
			table_for_one_metallicity = np.core.records.fromarrays(list_of_arrays,names=names)
			table_for_one_metallicity['Mass'] = np.array(list_of_masses)
			table_for_one_metallicity['mass_in_remnants'] = np.array(list_of_remnant)
			for i,item in enumerate(elements):
				for j,jtem in enumerate(list_of_masses):
					table_for_one_metallicity[item][j] = y['xi'][i+j*len_elements]
			for i,item in enumerate(table_for_one_metallicity["Mass"]):
				table_for_one_metallicity['mass_in_remnants'][i] /= item
				table_for_one_metallicity['unprocessed_mass_in_winds'][i] = 1.- table_for_one_metallicity['mass_in_remnants'][i]
				temp = sum(list(table_for_one_metallicity[i])[3:])
				for j,jtem in enumerate(elements):
					table_for_one_metallicity[jtem][i] /= temp
			for i,item in enumerate(elements):
				table_for_one_metallicity[item] -= initial_abundances[item][0]
			yield_tables[metallicity] = table_for_one_metallicity[::-1]
		self.masses = table_for_one_metallicity['Mass'][::-1]
		self.elements = elements
		self.table = yield_tables

	def Karakas_net_yield(self):
		'''
		loading the yield table of Karakas 2010.
		'''
		import numpy.lib.recfunctions as rcfuncs


		DATADIR = localpath + 'input/yields/Karakas2010'
		if not os.path.exists(DATADIR):
			os.mkdir(DATADIR)

		MASTERFILE = '{}/karakas_yields'.format(DATADIR)

		def _download_karakas():
			"""
			Downloads Karakas yields from Vizier.
			"""
			#url = 'http://zenodo.org/record/12800/files/dartmouth.h5'
			url = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J%2FMNRAS%2F403%2F1413'
			import urllib
			print('Downloading Karakas 2010 yield tables from Vizier (should happen only at the first time)...')
			if os.path.exists(MASTERFILE):
				os.remove(MASTERFILE)
			urllib.urlretrieve(url,MASTERFILE)

			import tarfile
			tar = tarfile.open(MASTERFILE)
			tar.extractall(path=DATADIR)
			tar.close()

		if not os.path.exists(MASTERFILE):
			_download_karakas()



		tdtype =   [('imass',float),('metallicity',float),('fmass',float),('species1','|S4'),('A',int),('net_yield',float),('ejected_mass',float),('initial_wind',float),('average_wind',float),('initial_mass_fraction',float),('production_factor',float)]
		metallicity_list = [0.02, 0.008, 0.004 ,0.0001]
		self.metallicities = metallicity_list
		


		tables = []
		for i,item in enumerate(metallicity_list):
			y = np.genfromtxt('%s/tablea%d.dat' %(DATADIR,i+2), dtype = tdtype, names = None)
			## Python3 need transformation between bytes and strings
			element_list2 = []
			for j,jtem in enumerate(y['species1']):
					element_list2.append(jtem.decode('utf8'))
			y = rcfuncs.append_fields(y,'species',element_list2,usemask = False)


			tables.append(y)
		

		### easy to extend to other species just make a new list of isotopes (see karakas tables)
		### and then also extend the indexing variable. 
		### The choice for specific elements can be done later when just using specific species
		hydrogen_list = ['n','p','d']
		helium_list = ['he3','he4']
		lithium_list = ['li7','be7','b8']
		
		carbon_list = ['c12','c13','n13']
		nitrogen_list = ['n14','n15','c14','o14','o15']
		oxygen_list = [ 'o16','o17','o18','f17','f18']
		fluorin_list = ['ne19','f19','o19']
		neon_list = ['ne20','ne21','ne22','f20','na21','na22']
		sodium_list = ['na23','ne23','mg23']
		magnesium_list = ['mg24','mg25','mg26','al-6','na24','al25']
		aluminium_list = ['mg27','al*6','al27','si27']
		silicon_list = ['al28','si28','si29','si30','p29','p30']
		phosphorus_list = ['si31','si32','si33','p31']
		sulfur_list = ['s32','s33','s34','p32','p33','p34']
		chlorine_list = ['s35']
		iron_list = ['fe54', 'fe56','fe57','fe58']
		manganese_list = ['fe55']
		cobalt_list = ['ni59','fe59','co59']
		nickel_list = ['ni58','ni60','ni61','ni62','co60','co61','fe60','fe61']


		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list

		#indexing['S_el'] = ni_to_bi
		
		self.elements = list(indexing.keys())
		#### little fix for karakas tablea5.dat: 6.0 M_sun is written two times. We chose the first one
		#tables[3]['imass'][-77:] = 6.5 # this is the fix if the second 6msun line was interpreted as 6.5 msun
		tables[3] = tables[3][:-77]

		#### making the general feedback table with yields for the individual elements
		### loop for the different metallicities
		yield_tables = {}
		for metallicity_index,metallicity in enumerate(metallicity_list[:]):
			### loop for the different elements
			yields_002 = {}
			for i,item1 in enumerate(indexing):
				unique_masses = len(np.unique(tables[metallicity_index]['imass']))
				element = np.zeros((unique_masses,), dtype=[('imass',float),('species','|S4'),('fmass',float),('net_yield',float),('ejected_mass',float),('initial_mass_fraction',float),('initial_wind',float),('average_wind',float),('production_factor',float)])
				for j,item in enumerate(indexing[item1]):
					cut = np.where(tables[metallicity_index]['species']==item)
					temp = tables[metallicity_index][cut]
					if j == 0:
						element['imass'] = temp['imass']
						element['fmass'] = temp['fmass']
						element['species'] = temp['species'] ### just for test purposes
					element['net_yield'] += temp['net_yield']
					element['ejected_mass'] += temp['ejected_mass']
					element['initial_mass_fraction'] += temp['initial_mass_fraction']
					element['initial_wind'] += temp['initial_wind']
					element['average_wind'] += temp['average_wind']
					element['production_factor'] += temp['production_factor']

				yields_002[item1] = element
			yield_tables[metallicity] = yields_002
		
		self.masses = np.unique(tables[0]['imass']) ## table a3 and a4 and maybe a5 are missing 6.5 Msun its probably easier to skip the 6.5 Msun entries altogether for interpolation reasons

		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(metallicity_list[:]):
			yields_for_one_metallicity = yield_tables[metallicity]
			final_mass_name_tag = 'mass_in_remnants'
			additional_keys = ['Mass',final_mass_name_tag,'unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			if metallicity == 0.02: #or metallicity == 0.0001:
				base = np.zeros(len(self.masses))
			else:
				base = np.zeros(len(self.masses)-1)
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)

			yield_tables_final_structure_subtable['Mass'] = yields_for_one_metallicity[self.elements[0]]['imass']
			yield_tables_final_structure_subtable[final_mass_name_tag] = np.divide(yields_for_one_metallicity[self.elements[0]]['fmass'],yield_tables_final_structure_subtable['Mass'])#np.divide(yields_for_one_metallicity[self.elements[0]]['fmass'],yield_tables_final_structure_subtable['Mass'])
			temp = np.zeros_like(yield_tables_final_structure_subtable['Mass'])
			for i,item in enumerate(self.elements):
				################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
				yield_tables_final_structure_subtable[item] = np.divide(yields_for_one_metallicity[item]['net_yield'],yield_tables_final_structure_subtable['Mass'])
				temp += yield_tables_final_structure_subtable[item]
			yield_tables_final_structure_subtable['unprocessed_mass_in_winds'] = 1. - (yield_tables_final_structure_subtable[final_mass_name_tag] + temp )

			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable[::-1]
		self.table = yield_tables_final_structure

	def one_parameter(self, elements, element_fractions):
		"""
		Another problem: He and the remnant mass fraction is not constrained in the APOGEE data. Maybe these can be constrained externally by yield sets or cosmic abundance standard or solar abundances.
		"""
		self.metallicities = [0.01]
		self.masses = np.array([3])
		self.elements = elements 

		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		

		additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
		names = additional_keys + self.elements
		base = np.zeros(len(self.masses))
		list_of_arrays = []
		
		for i in range(len(names)):
			list_of_arrays.append(base)
		yield_table = np.core.records.fromarrays(list_of_arrays,names=names)
		yield_table['Mass'] = self.masses
		yield_table['mass_in_remnants'] = 0.27
		yield_table['unprocessed_mass_in_winds'] = 1 - yield_table['mass_in_remnants']
		for i,item in enumerate(self.elements[1:]):
			yield_table[item] = element_fractions[i+1]
		yield_table['H'] = -sum(element_fractions[1:])


		yield_tables_final_structure[self.metallicities[0]] = yield_table
		self.table = yield_tables_final_structure


class Hypernova_feedback(object):
	def __init__(self):
		"""
		this is the object that holds the feedback table for Hypernova
		"""
	def Nomoto2013(self):
		'''
		Nomoto2013 sn2 yields from 13Msun onwards
		'''
		import numpy.lib.recfunctions as rcfuncs

		dt = np.dtype('a13,f8,f8,f8,f8')
		yield_tables = {}
		self.metallicities = [0.0500,0.0200,0.0080,0.0040,0.0010]
		self.masses = np.array((20,25,30,40))
		z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/hn_z=0.0200.dat',dtype=dt,names = True)
		
		yield_tables_dict = {}
		for item in self.metallicities:
			z = np.genfromtxt(localpath + 'input/yields/Nomoto2013/hn_z=%.4f.dat' %(item),dtype=dt,names = True)
			yield_tables_dict[item]=z
		#########################
		hydrogen_list = ['H__1','H__2']
		helium_list = ['He_3','He_4']
		lithium_list = ['Li_6','Li_7']
		berillium_list = ['Be_9']
		boron_list = ['B_10','B_11']
		carbon_list = ['C_12','C_13']
		nitrogen_list = ['N_14','N_15']
		oxygen_list = ['O_16','O_17','O_18']
		fluorin_list = ['F_19']
		neon_list = ['Ne20','Ne21','Ne22']
		sodium_list = ['Na23']
		magnesium_list = ['Mg24','Mg25','Mg26']
		aluminium_list = ['Al27']
		silicon_list = ['Si28','Si29','Si30']
		phosphorus_list = ['P_31']
		sulfur_list = ['S_32','S_33','S_34','S_36']
		chlorine_list = ['Cl35','Cl37']
		argon_list = ['Ar36','Ar38','Ar40']
		potassium_list = ['K_39','K_41']
		calcium_list = ['K_40','Ca40','Ca42','Ca43','Ca44','Ca46','Ca48']
		scandium_list = ['Sc45']
		titanium_list = ['Ti46','Ti47','Ti48','Ti49','Ti50']
		vanadium_list = ['V_50','V_51']
		chromium_list = ['Cr50','Cr52','Cr53','Cr54']
		manganese_list = ['Mn55']
		iron_list = ['Fe54', 'Fe56','Fe57','Fe58']
		cobalt_list = ['Co59']
		nickel_list = ['Ni58','Ni60','Ni61','Ni62','Ni64']


		copper_list = ['Cu63','Cu65']
		zinc_list = ['Zn64','Zn66','Zn67','Zn68','Zn70']
		gallium_list = ['Ga69','Ga71']
		germanium_list = ['Ge70','Ge72','Ge73','Ge74']

		indexing = {}
		indexing['H'] = hydrogen_list
		indexing['He'] = helium_list
		indexing['Li'] = lithium_list
		indexing['Be'] = berillium_list
		indexing['B'] = boron_list
		
		indexing['C'] = carbon_list
		indexing['N'] = nitrogen_list
		indexing['O'] = oxygen_list
		indexing['F'] = fluorin_list
		indexing['Ne'] = neon_list
		indexing['Na'] = sodium_list
		indexing['Mg'] = magnesium_list
		indexing['Al'] = aluminium_list
		indexing['Si'] = silicon_list
		indexing['P'] = phosphorus_list
		indexing['S'] = sulfur_list
		indexing['Cl'] = chlorine_list
		indexing['Ar'] = argon_list
		indexing['K'] = potassium_list
		indexing['Ca'] = calcium_list
		indexing['Sc'] = scandium_list
		indexing['Ti'] = titanium_list
		indexing['V'] = vanadium_list
		indexing['Cr'] = chromium_list
		indexing['Mn'] = manganese_list
		indexing['Fe'] = iron_list
		indexing['Co'] = cobalt_list
		indexing['Ni'] = nickel_list
		indexing['Cu'] = copper_list
		indexing['Zn'] = zinc_list
		indexing['Ga'] = gallium_list
		indexing['Ge'] = germanium_list

		self.elements = list(indexing.keys())
		### restructuring the tables such that it looks like the sn2 dictionary: basic_agb[metallicicty][element]
		yield_tables_final_structure = {}
		for metallicity_index,metallicity in enumerate(self.metallicities):
			yields_for_one_metallicity = yield_tables_dict[metallicity]
			## Python3 need transformation between bytes and strings
			element_list2 = []
			for j,item in enumerate(yields_for_one_metallicity['M']):
					element_list2.append(item.decode('utf8'))
			yields_for_one_metallicity = rcfuncs.append_fields(yields_for_one_metallicity,'element',element_list2,usemask = False)

			additional_keys = ['Mass','mass_in_remnants','unprocessed_mass_in_winds']
			names = additional_keys + self.elements
			base = np.zeros(len(self.masses))
			list_of_arrays = []
			
			for i in range(len(names)):
				list_of_arrays.append(base)
			yield_tables_final_structure_subtable = np.core.records.fromarrays(list_of_arrays,names=names)
			yield_tables_final_structure_subtable['Mass'] = self.masses
			temp1 = np.zeros(len(self.masses))
			for i in range(len(self.masses)):
				temp1[i] = yields_for_one_metallicity[0][i+1]
			yield_tables_final_structure_subtable['mass_in_remnants'] = np.divide(temp1,self.masses)
			
			for i,item in enumerate(self.elements):
				yield_tables_final_structure_subtable[item] = 0
				for j,jtem in enumerate(indexing[item]):
						################### here we can change the yield that we need for processing. normalising 'ejected_mass' with the initial mass to get relative masses
						line_of_one_element = yields_for_one_metallicity[np.where(yields_for_one_metallicity['element']==jtem)][0]
						temp1 = np.zeros(len(self.masses))
						for i in range(len(self.masses)):
							temp1[i] = line_of_one_element[i+1]
						yield_tables_final_structure_subtable[item] += np.divide(temp1,self.masses)

			for i in range(len(self.masses)):
				yield_tables_final_structure_subtable['unprocessed_mass_in_winds'][i] = (1-yield_tables_final_structure_subtable['mass_in_remnants'][i]-sum(yield_tables_final_structure_subtable[self.elements][i]))#yields_for_one_metallicity[0][21]#			
			yield_tables_final_structure[metallicity] = yield_tables_final_structure_subtable#[::-1]
		self.table = yield_tables_final_structure
