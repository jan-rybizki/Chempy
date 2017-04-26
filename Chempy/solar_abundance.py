import numpy as np 
from .making_abundances import abundance_to_mass_fraction
import numpy.lib.recfunctions as rcfuncs
from . import localpath

class solar_abundances(object):
	'''
	The solar abundances class. Holds information on the element names, symbols, mass numbers and photospheric abundances.
	'''
	def __init__(self):    
		'''
		Upon initialization the element names and numbers and masses are loaded and the solar table loaded. The values are filled after using one of the methods which will load solar abundance literature values into the table.
		'''
		self.table = np.load(localpath + 'input/elemental_table.npy')
		## Python3 need transformation between bytes and strings
		element_list = []
		for j,jtem in enumerate(self.table['Symbol']):
			element_list.append(jtem.decode('utf8'))
		self.all_elements = element_list
		self.table = rcfuncs.drop_fields(self.table,'Symbol',usemask = False)
		self.table = rcfuncs.append_fields(self.table,'Symbol',element_list,usemask = False)
		self.all_element_numbers = list(self.table['Number'])
		self.all_element_masses = list(self.table['Mass'])
		self.dimensions = ['Symbol','Number','Mass']
		self.extra = ['photospheric','error']
		self.names = self.dimensions + self.extra
		self.base = np.zeros(len(self.all_elements))

	def Lodders09(self):#								    O 8.73              Mg 7.54 Fe 7.46
		'''
		Photospheric abundances and errors are loaded from Lodders+ 2009. Also the elment fractions are calculated together with X, Y and Z the Hydrogen, Helium and metallicity fraction.
		'''
		abundances = [12.00,10.93,3.28,1.32,2.81,8.39,7.86,8.73,4.44,8.05,6.29,7.54,6.46,\
		7.53,5.45,7.16,5.25,6.50,5.11,6.31,3.07,4.93,3.99,5.65,5.50,7.46,\
		4.90,6.22,4.27,4.65,3.10,3.59,2.32,3.36,2.56,3.28,2.38,2.90,2.20,\
		2.57,1.42,1.94,1.78,1.10,1.67,1.22,1.73,0.78,2.09,1.03,2.20,1.57,\
		2.27,1.10,2.18,1.19,1.60,0.77,1.47,0.96,0.53,1.09,0.34,1.14,0.49,\
		0.95,0.14,0.94,0.11,0.73,-0.14,0.67,0.28,1.37,1.36,1.64,0.82,1.19,\
		0.79,2.06,0.67,0.08,-0.52]

		errors = [0., 2, 5, 3, 4, 4, 12, 7, 6, 10, 4, 6,  7,\
		6, 5, 2, 6, 10, 4, 2, 2, 3, 3, 2, 1, 8,\
		8, 4, 4, 4, 2, 6, 4, 3, 6, 8, 3, 3, 4,\
		4, 7, 6, 3, 13, 4, 2, 3, 3, 6, 6, 3, 8, 8,\
		2, 7, 2, 6, 5, 2, 2, 4, 6, 3, 6, 3, 6,\
		3, 2, 2, 2, 4, 4, 4, 3, 6, 3, 4, 8,\
		3, 3, 4, 3, 3]

		errors = list(np.array(errors)*0.01)
		numbers = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12,  13,\
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,\
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,\
		40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,\
		54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67,\
		68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,\
		81, 82, 83, 90, 92]

		for i,item in enumerate(numbers):
			self.table['photospheric'][np.where(self.table['Number']==item)] = abundances[i]
			self.table['error'][np.where(self.table['Number']==item)] = errors[i]

		self.fractions = abundance_to_mass_fraction(np.hstack(self.all_elements),np.hstack(self.all_element_masses),self.table['photospheric'],self.table['photospheric'],np.hstack(self.all_elements))
		self.x = self.fractions[0]
		self.y = self.fractions[1]
		self.z = sum(self.fractions[2:])
		self.errors = errors

	def Asplund09(self):#									O 8.69             Mg 7.60 Fe 7.50
		'''
		Photospheric abundances and errors are loaded from Asplund+ 2009. Also the elment fractions are calculated together with X, Y and Z the Hydrogen, Helium and metallicity fraction.
		'''
		abundances = [12.00,10.93,3.26,1.30,2.79,8.43,7.83,8.69,4.42,7.93,6.24,7.60,6.45,\
		7.51,5.41,7.12,5.23,6.40,5.03,6.34,3.15,4.95,3.93,5.64,5.43,7.50,\
		4.99,6.22,4.19,4.56,3.04,3.65,2.30,3.34,2.54,3.25,2.52,2.87,2.21,\
		2.57,1.46,1.88,1.75,0.91,1.57,0.94,1.71,0.80,2.04,1.01,2.18,1.55,\
		2.24,1.08,2.18,1.10,1.58,0.72,1.42,0.96,0.52,1.07,0.30,1.10,0.48,\
		0.92,0.10,0.84,0.10,0.85,-0.12,0.85,0.26,1.40,1.38,1.62,0.92,1.17,\
		0.90,1.75,0.65,0.02,-0.54]

		
		errors = [0, 1, 5, 3, 4, 5, 5, 5 , 6, 10, 4, 4, 3,\
		3, 3, 3, 6, 13, 9, 4, 4, 5, 8, 4, 4, 4,\
		7, 4, 4, 5, 9, 10, 4, 3, 6, 6, 10, 7, 5,\
		4, 4, 8, 8, 10, 10, 10, 3, 20, 10, 6, 3, 8,\
		6, 2, 9, 4, 4, 4, 4, 4, 4, 4, 10, 4, 11,\
		5, 4, 11, 9, 4, 4, 12, 4, 8, 7, 3, 10, 8,\
		20, 10, 4, 10, 3]

		errors = list(np.array(errors)*0.01)
		
		numbers = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12,  13,\
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,\
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,\
		40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,\
		54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67,\
		68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,\
		81, 82, 83, 90, 92]
		
		for i,item in enumerate(numbers):
			self.table['photospheric'][np.where(self.table['Number']==item)] = abundances[i]
			self.table['error'][np.where(self.table['Number']==item)] = errors[i]		
		self.fractions = abundance_to_mass_fraction(np.hstack(self.all_elements),np.hstack(self.all_element_masses),self.table['photospheric'],self.table['photospheric'],np.hstack(self.all_elements))

		self.x = self.fractions[0]
		self.y = self.fractions[1]
		self.z = sum(self.fractions[2:])


	def AG89(self):#									O 8.69             Mg 7.60 Fe 7.50
		'''
		Photospheric abundances and errors are loaded from Anders & Grevesse 1989. Also the elment fractions are calculated together with X, Y and Z the Hydrogen, Helium and metallicity fraction.
		'''
		abundances = [12.00,10.99,3.31,1.42,2.88,8.56,8.05,8.93,4.56,8.09,6.33,7.58,6.47,\
		7.55,5.45,7.21,5.27,6.56,5.12,6.36,3.10,4.99,4.00,5.67,5.39,7.67,\
		4.92,6.25,4.21,4.60,3.13,3.41,2.37,3.35,2.63,3.23,2.60,2.90,2.24,\
		2.60,1.42,1.92,1.84,1.12,1.69,0.94,1.86,0.82,2.14,1.04,2.24,1.51,\
		2.23,1.12,2.13,1.22,1.55,0.71,1.50,1.00,0.51,1.12,0.33,1.15,0.50,\
		0.93,0.13,0.95,0.12,0.73,0.13,0.68,0.27,1.38,1.37,1.68,0.83,1.09,\
		0.82,1.85,0.71,0.08,-0.49]

		
		errors = [0, 3, 4, 4, 4, 4, 4, 4, 3, 10, 3, 5, 7,\
		5, 4, 6, 6, 10, 13, 2, 4, 2, 2, 3, 3, 3,\
		4, 4, 4, 8, 3, 14, 5, 3, 8, 7, 3, 6, 3,\
		3, 6, 5, 7, 12, 4, 1, 15, 3, 4, 7, 4, 8,\
		8, 2, 5, 9, 20, 8, 6, 8, 8, 4, 1, 1, 1,\
		6, 1, 1, 1, 1, 1, 2, 4, 3, 3, 3, 6, 5,\
		4, 3, 3, 2, 4]

		errors = list(np.array(errors)*0.01)
		
		numbers = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12,  13,\
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,\
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,\
		40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,\
		54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67,\
		68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,\
		81, 82, 83, 90, 92]
		
		for i,item in enumerate(numbers):
			self.table['photospheric'][np.where(self.table['Number']==item)] = abundances[i]
			self.table['error'][np.where(self.table['Number']==item)] = errors[i]		
		self.fractions = abundance_to_mass_fraction(np.hstack(self.all_elements),np.hstack(self.all_element_masses),self.table['photospheric'],self.table['photospheric'],np.hstack(self.all_elements))

		self.x = self.fractions[0]
		self.y = self.fractions[1]
		self.z = sum(self.fractions[2:])

	def Asplund05_pure_solar(self):
		'''
		Photospheric abundances and errors are loaded from Asplund+ 2005. Also the elment fractions are calculated together with X, Y and Z the Hydrogen, Helium and metallicity fraction.
		It is not sure for which elements from Asplund+ 2005 the apogee consortium has used the photospheric or the meteoritic abundances.
		Therefore I try here to use only the photospheric except for elements without photospheric values.
		'''
						#									O 8.69             Mg 7.60 Fe 7.50
		abundances = [12.00,10.93,1.05,1.38,2.70,8.39,7.78,8.66,4.56,7.84,6.17,7.53,6.37,\
		7.51,5.36,7.14,5.50,6.18,5.08,6.31,3.05,4.90,4.00,5.64,5.39,7.45,\
		4.92,6.23,4.21,4.60,2.88,3.58,2.29,3.33,2.56,3.28,2.60,2.92,2.21,\
		2.59,1.42,1.92,1.84,1.12,1.69,0.94,1.77,1.60,2.00,1.00,2.19,1.51,\
		2.27,1.07,2.17,1.13,1.58,0.71,1.45,1.01,0.52,1.12,0.28,1.14,0.51,\
		0.93,0.00,1.08,0.06,0.88,-0.17,1.11,0.23,1.45,1.38,1.64,1.01,1.13,\
		0.90,2.00,0.65,0.06,-0.52]
		# these errors just copied from Aslpund 2009
		errors = [0, 1, 5, 3, 4, 5, 5, 5 , 6, 10, 4, 4, 3,\
		3, 3, 3, 6, 13, 9, 4, 4, 5, 8, 4, 4, 4,\
		7, 4, 4, 5, 9, 10, 4, 3, 6, 6, 10, 7, 5,\
		4, 4, 8, 8, 10, 10, 10, 3, 20, 10, 6, 3, 8,\
		6, 2, 9, 4, 4, 4, 4, 4, 4, 4, 10, 4, 11,\
		5, 4, 11, 9, 4, 4, 12, 4, 8, 7, 3, 10, 8,\
		20, 10, 4, 10, 3]

		errors = list(np.array(errors)*0.01)

		
		numbers = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12,  13,\
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,\
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,\
		40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,\
		54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67,\
		68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,\
		81, 82, 83, 90, 92]
		
		for i,item in enumerate(numbers):
			self.table['photospheric'][np.where(self.table['Number']==item)] = abundances[i]
			self.table['error'][np.where(self.table['Number']==item)] = errors[i]
		self.fractions = abundance_to_mass_fraction(np.hstack(self.all_elements),np.hstack(self.all_element_masses),self.table['photospheric'],self.table['photospheric'],np.hstack(self.all_elements))

		self.x = self.fractions[0]
		self.y = self.fractions[1]
		self.z = sum(self.fractions[2:])

	def Asplund05_apogee_correction(self):
		'''
		Photospheric abundances and errors are loaded from Asplund+ 2005 but corrected for the APOGEE scale. Also the elment fractions are calculated together with X, Y and Z the Hydrogen, Helium and metallicity fraction.
		It is not sure for which elements from Asplund+ 2005 the apogee consortium has used the photospheric or the meteoritic abundances.
		Therefore I try here to use only the photospheric except for elements without photospheric values.
		'''
		
		"""
		After an email from Carlos Allende Prieto only the synthetic spectra are made with Asplund 2005. If trying to normalise to zero we need to see the results of a solar twin with the apogee pipeline.
		Carlos send me the results for Vesta
		C  				N  				O  			Na 			Mg 			Al 
		0.0349240    0.0772880   0.00191380    0.0307620    0.0467750 -0.0318830
   		Si 					S  			K  			Ca 			Ti 			V  
   		0.00338050     0.235970   -0.0714790  -0.00784230 -0.0997200     0.105860
    	Mn 				Fe 				Ni
    	0.0379600    0.0133580   0.00974230

    	For C, N, O, Mg, Si, S, Ca and Ti, the abundances are [X/METALS] (where
		METALS is [Fe/H] measured from all metal lines, and as listed above [M/H]=0.026).
		For the rest of the elements (Na, Al, K, V, Mg, Fe and Ni), the entries are [X/H]. 

		For example, the Vesta carbon abundance we derive is [C/H]=[C/METALS]-[METALS/H]=0.035 - 0.026 = 0.009, and since the Asplund et al. value for the Sun is log(epsilon)+12=8.39, ours is 8.39+0.009=8.40.
		
		if we do this for all the elements we arrive at new normalisations:
					metals              Vesta offset  Asplund05 corrected Asplund09
		C  = 0.035 - 0.026 = 0.009  --> C  = 0.009  + 8.39 = 		8.40  ~ 8.43
		N  = 0.077 - 0.026 = 0.051  --> N  = 0.051  + 7.78 = 		7.83  ~ 7.83
		O  = 0.002 - 0.026 = -0.024 --> O  = -0.024 + 8.66 = 		8.64  ~ 8.69
		Mg = 0.047 - 0.026 = 0.021  --> Mg = 0.021  + 7.53 = 		7.55  ~ 7.60
		Si = 0.003 - 0.026 = -0.023 --> Si = -0.023 + 7.51 = 		7.49  ~ 7.51
		S  = 0.236 - 0.026 = 0.21   --> S  = 0.21   + 7.14 = 		7.35  ~ 7.12  !
		Ca = -0.008- 0.026 = -0.034 --> Ca = -0.034 + 6.31 = 		6.28  ~ 6.34
		Ti = -0.100- 0.026 = -0.126 --> Ti = -0.126 + 4.90 = 		4.77  ~ 4.95  !

		Na 				   = 0.031  --> Na = 0.031  + 6.17 = 		6.20  ~ 6.24
		Al 				   = -0.032 --> Al = -0.032 + 6.37 = 		6.34  ~ 6.45  !
		K  				   = -0.071 --> K  = -0.071 + 5.08 = 		5.01  ~ 5.03
		V  				   = 0.106  --> V  = 0.106  + 4.00 = 		4.11  ~ 3.93  !
		Mn 				   = 0.038  --> Mn = 0.038  + 5.39 = 		5.43  ~ 5.43
		Fe 				   = 0.013  --> Fe = 0.013  + 7.45 = 		7.46  ~ 7.50
		Ni 				   = 0.010  --> Ni = 0.010  + 6.23 = 		6.24  ~ 6.22


		However, one needs to keep in mind that the vast majority of APOGEE stars are in the range 3500<Teff<4500 K (and most are giants), so they are cooler than the Sun, and systematic errors between then Sun and them are most likely. 
		"""
						#									O 8.69             Mg 7.60 Fe 7.50
		abundances = [12.00,10.93,1.05,1.38,2.70,8.40,7.83,8.64,4.56,7.84,6.20,7.55,6.34,\
		7.49,5.36,7.35,5.50,6.18,5.01,6.28,3.05,4.77,4.11,5.64,5.43,7.46,\
		4.92,6.24,4.21,4.60,2.88,3.58,2.29,3.33,2.56,3.28,2.60,2.92,2.21,\
		2.59,1.42,1.92,1.84,1.12,1.69,0.94,1.77,1.60,2.00,1.00,2.19,1.51,\
		2.27,1.07,2.17,1.13,1.58,0.71,1.45,1.01,0.52,1.12,0.28,1.14,0.51,\
		0.93,0.00,1.08,0.06,0.88,-0.17,1.11,0.23,1.45,1.38,1.64,1.01,1.13,\
		0.90,2.00,0.65,0.06,-0.52]

		# these errors just copied from Aslpund 2009
		errors = [0, 1, 5, 3, 4, 5, 5, 5 , 6, 10, 4, 4, 3,\
		3, 3, 3, 6, 13, 9, 4, 4, 5, 8, 4, 4, 4,\
		7, 4, 4, 5, 9, 10, 4, 3, 6, 6, 10, 7, 5,\
		4, 4, 8, 8, 10, 10, 10, 3, 20, 10, 6, 3, 8,\
		6, 2, 9, 4, 4, 4, 4, 4, 4, 4, 10, 4, 11,\
		5, 4, 11, 9, 4, 4, 12, 4, 8, 7, 3, 10, 8,\
		20, 10, 4, 10, 3]

		errors = list(np.array(errors)*0.01)
		
		numbers = [1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12,  13,\
		14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,\
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,\
		40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,\
		54, 55, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67,\
		68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,\
		81, 82, 83, 90, 92]
		
		for i,item in enumerate(numbers):
			self.table['photospheric'][np.where(self.table['Number']==item)] = abundances[i]
			self.table['error'][np.where(self.table['Number']==item)] = errors[i]
		self.fractions = abundance_to_mass_fraction(np.hstack(self.all_elements),np.hstack(self.all_element_masses),self.table['photospheric'],self.table['photospheric'],np.hstack(self.all_elements))

		self.x = self.fractions[0]
		self.y = self.fractions[1]
		self.z = sum(self.fractions[2:])
