import numpy as np 

def slope_imf(x,p1,p2,p3,kn1,kn2):
	'''
	Is calculating a three slope IMF
	
	INPUT:

	   x = An array of masses for which the IMF should be calculated
	
	   p1..p3 = the slopes of the power law
	
	   kn1, kn2 = Where the breaks of the power law are

	OUTPUT:
	
	   An array of frequencies matching the mass base array x
	'''
	if(x > kn2):
		t = (pow(kn2,p2)/pow(kn2,p3))*pow(x,p3+1)
	elif (x < kn1):
		t = (pow(kn1,p2)/pow(kn1,p1))*pow(x,p1+1)
	else:
		t = pow(x,p2+1)
	return t


def lifetime(m,Z):
	"""
	here we will calculate the MS lifetime of the star after Argast et al., 2000, A&A, 356, 873
	INPUT:
	
	   m = mass in Msun
	
	   Z = metallicity in Zsun
	

	OUTPUT:
	
	   returns the lifetime of the star in Gyrs
	"""
	lm = np.log10(m)
	a0 =  3.79 + 0.24*Z
	a1 = -3.10 - 0.35*Z
	a2 =  0.74 + 0.11*Z
	tmp = a0 + a1*lm + a2*lm*lm
	return np.divide(np.power(10,tmp),1000)


class IMF(object):
	'''
	This class represents the IMF normed to 1 in units of M_sun. 
	
	Input for initialisation:

	   mmin = minimal mass of the IMF
	
	   mmax = maximal mass of the IMF
	
	   intervals = how many steps inbetween mmin and mmax should be given
	
	Then one of the IMF functions can be used
	
	   self.x = mass base
	
	   self.dn = the number of stars at x 
	
	   self.dm = the masses for each mass interval x
	'''
	def __init__(self, mmin = 0.08 , mmax = 100., intervals = 5000):
		self.mmin = mmin
		self.mmax = mmax
		self.intervals = intervals
		self.x = np.linspace(mmin,mmax,intervals)
		self.dx = self.x[1]-self.x[0]

	def normed_3slope(self,paramet = (-1.3,-2.2,-2.7,0.5,1.0)):
		'''
		Three slope IMF, Kroupa 1993 as a default
		'''
		s1,s2,s3,k1,k2 = paramet
		u = np.zeros_like(self.x)
		v = np.zeros_like(self.x)
		for i in range(len(self.x)):
			u[i] = slope_imf(self.x[i],s1,s2,s3,k1,k2)
		v = np.divide(u,self.x)
		self.dm = np.divide(u,sum(u))
		self.dn = np.divide(self.dm,self.x)
		return(self.dm,self.dn)

	def Chabrier_1(self, paramet = (0.69, 0.079, -2.3)):
		'''
		Chabrier IMF from Chabrier 2003 equation 17 field IMF with variable high mass slope and automatic normalisation
		'''
		sigma, m_c, expo = paramet
		dn = np.zeros_like(self.x)
		for i in range(len(self.x)):
			if self.x[i] <= 1:
				index_with_mass_1 = i 
				dn[i] = (1. / float(self.x[i])) * np.exp(-1*(((np.log10(self.x[i] / m_c))**2)/(2*sigma**2)))
			else:
				dn[i] = (pow(self.x[i],expo))
		# Need to 'attach' the upper to the lower branch
		derivative_at_1 = dn[index_with_mass_1] - dn[index_with_mass_1 - 1]
		target_y_for_m_plus_1 = dn[index_with_mass_1] + derivative_at_1
		rescale = target_y_for_m_plus_1 / dn[index_with_mass_1 + 1] 
		dn[np.where(self.x>1.)] *= rescale
		# Normalizing to 1 in mass space
		self.dn = np.divide(dn,sum(dn))
		dm = dn*self.x
		self.dm = np.divide(dm,sum(dm))
		self.dn = np.divide(self.dm,self.x)
		return(self.dm,self.dn)       
	
	def Chabrier_2(self,paramet = (22.8978, 716.4, 0.25,-2.3)):
		'''
		Chabrier IMF from Chabrier 2001, IMF 3 = equation 8 parameters from table 1
		'''

		A,B,sigma,expo = paramet
		expo -= 1. ## in order to get an alpha index normalisation
		dn = np.zeros_like(self.x)
		for i in range(len(self.x)):
			dn[i] = A*(np.exp(-pow((B/self.x[i]),sigma)))*pow(self.x[i],expo)
		self.dn = np.divide(dn,sum(dn))
		dm = dn*self.x
		self.dm = np.divide(dm,sum(dm))
		self.dn = np.divide(self.dm,self.x)
		return(self.dm,self.dn)

	def salpeter(self, alpha = (2.35)):
		'''
		Salpeter IMF

		Input the slope of the IMF
		'''
		self.alpha = alpha
		temp = np.power(self.x,-self.alpha)
		norm = sum(temp)
		self.dn = np.divide(temp,norm)
		u = self.dn*self.x
		self.dm = np.divide(u,sum(u))
		self.dn = np.divide(self.dm,self.x)
		return (self.dm,self.dn)

	def BrokenPowerLaw(self, paramet):
		breaks,slopes = paramet
		if len(breaks) != len(slopes)-1:
			print("error in the precription of the power law. It needs one more slope than break value")
		else:
			dn = np.zeros_like(self.x)
			self.breaks = breaks
			self.slopes = slopes
			self.mass_range = np.hstack((self.mmin,breaks,self.mmax))
			for i,item in enumerate(self.slopes):
				cut = np.where(np.logical_and(self.x>=self.mass_range[i],self.x<self.mass_range[i+1]))
				dn[cut] = np.power(self.x[cut],item)
				if i != 0:
					renorm = np.divide(last_dn,dn[cut][0])
					dn[cut] = dn[cut]*renorm
				last_dn = dn[cut][-1]
				last_x = self.x[cut][-1]
			self.dn = np.divide(dn,sum(dn))
			u = self.dn*self.x
			self.dm = np.divide(u,sum(u))
			self.dn = np.divide(self.dm,self.x)

	def imf_mass_fraction(self,mlow,mup):
		'''
		Calculates the mass fraction of the IMF sitting between mlow and mup
		'''
		norm = sum(self.dm)
		cut = np.where(np.logical_and(self.x>=mlow,self.x<mup))
		fraction = np.divide(sum(self.dm[cut]),norm)
		return(fraction)
	def imf_number_fraction(self,mlow,mup):
		'''
		Calculating the number fraction of stars of the IMF sitting between mlow and mup
		'''
		norm = sum(self.dn)
		cut = np.where(np.logical_and(self.x>=mlow,self.x<mup))
		fraction = np.divide(sum(self.dn[cut]),norm)
		return(fraction)
	def imf_number_stars(self,mlow,mup):
		cut = np.where(np.logical_and(self.x>=mlow,self.x<mup))
		number = sum(self.dn[cut])
		return(number)
	def stochastic_sampling(self, mass):
		'''
		The analytic IMF will be resampled according to the mass of the SSP.
		The IMF will still be normalised to 1

		Stochastic sampling is realised by fixing the number of expected stars and then drawing from the probability distribution of the number density
		Statistical properties are tested for this sampling and are safe: number of stars and masses converge.
		'''
		number_of_stars = int(round(sum(self.dn) * mass))
		self.dm_copy = np.copy(self.dm)
		self.dn_copy = np.copy(self.dn)

		#self.dn_copy = np.divide(self.dn_copy,sum(self.dn_copy))
		random_number = np.random.uniform(low = 0.0, high = sum(self.dn_copy), size = number_of_stars)
		self.dm = np.zeros_like(self.dm)
		self.dn = np.zeros_like(self.dn)
		'''
		### This could be favourable if the number of stars drawn is low compared to the imf resolution
		for i in range(number_of_stars):
			### the next line randomly draws a mass according to the number distribution of stars
			cut = np.where(np.abs(np.cumsum(self.dn_copy)-random_number[i])== np.min(np.abs(np.cumsum(self.dn_copy)-random_number[i])))
			x1 = self.x[cut][0]
			#self.dn[cut] += 0.5
			self.dn[cut[0]] += 1
			self.dm[cut[0]] += x1 + self.dx/2.
			t.append(x1 + self.dx/2.)
		'''
		counting = np.cumsum(self.dn_copy)
		for i in range(len(counting)-1):
			if i == 0:
				cut = np.where(np.logical_and(random_number>0.,random_number<=counting[i]))
			else:
				cut = np.where(np.logical_and(random_number>counting[i-1],random_number<=counting[i]))
			number_of_stars_in_mass_bin = len(random_number[cut])
			self.dm[i] = number_of_stars_in_mass_bin * self.x[i]
		self.dm = np.divide(self.dm,sum(self.dm))
		self.dn = np.divide(self.dm,self.x)
