import numpy as np
class SFR(object):
	'''
	The SFR class holds the star formation history and the time-steps of Chempy
	'''
	def __init__(self,start,end,time_steps):
		'''
		Upon initialization the time steps need to be provided

		INPUT:
		
		   start = beginning of the simulation
		
		   end = end of the simulation
		
		   time_steps = number of time_steps

		OUTPUT:
		
		   dt, timespan and t will be exposed by the class.
		'''
		self.start = start
		self.end = end
		self.time_steps = time_steps
		self.t = np.linspace(self.start, self.end, self.time_steps)
		self.dt = self.t[1]-self.t[0]
		self.timespan = self.end - self.start

	def model_A(self,S0 = 45.07488,t0 = 5.6,t1 = 8.2):
		'''
		This was the method to load the Just Jahreiss 2010 Model A from txt
		'''
		#Model A SFR from Just & Jahreiss 2010
		
		## this function can be used to read in the model A at different radii
		def read_model(time,r):
			radius = sum(r) * 0.5
			list_of_radii = np.array([4,5,6,7,8,9,10,11,12])
			radius = list_of_radii[np.where(np.abs(list_of_radii-radius)==np.min(np.abs(list_of_radii-radius)))][0]
			age_column = 'sfr_' + str(radius)
			x = np.genfromtxt('input/model/SFR_1.txt', names = True)
			time_model = x['timeGyr']
			age_distribution_model = x[age_column]

			return np.interp(time,time_model,age_distribution_model)
		self.sfr = (self.t + t0)/(self.t**2 + t1**2)**2
		self.sfr = np.divide(self.sfr,sum(self.sfr)/(np.divide(1.,self.dt)*S0))
  
	def doubly_peaked(self,S0 = 45.07488, peak_ratio = 1., decay = 2., t0 = 2., peak1t0 = 0.5, peak1sigma = 0.5):
		'''
		a doubly peaked SFR with quite a few parameters
		'''
		from scipy import signal
		from scipy.stats import norm, expon
		peak1 = norm.pdf(self.t,loc = peak1t0, scale = peak1sigma)
		peak1 = peak_ratio * np.divide(peak1,sum(peak1))
		peak2 = expon.pdf(self.t,loc = t0,scale = decay)
		peak2 = np.divide(peak2,sum(peak2))
		sig = peak1 + peak2
		self.sfr = sig
		self.sfr = np.divide(self.sfr,sum(self.sfr)/(np.divide(1.,self.dt)*S0))

	def gamma_function(self,S0 = 45.07488,a_parameter = 2, loc = 0, scale = 3):
		'''
		the gamma function for a_parameter = 2 and loc = 0 produces a peak at scale so we have a two parameter sfr.
		Later we can also release a to have a larger diversity in functional form.
		'''
		from scipy.stats import gamma
		self.sfr = gamma.pdf(self.t,a_parameter,loc,scale)
		self.sfr[np.where(self.sfr == 0.)] = np.min(self.sfr[np.where(self.sfr != 0.)])*0.01 ## So that no 0 sfr is there because sfr-related infall prescription fails in that case
		self.sfr = np.divide(self.sfr,sum(self.sfr)/(np.divide(1.,self.dt)*S0))
	def prescribed(self, mass_factor,name_of_file):
		'''
		a method to read in prescribed SFR from textfile
		x time is given in log years. our time is in linear Gyrs
		'''
		x = np.genfromtxt(name_of_file, names = True)
		x['time_l'] = np.power(10,x['time_l'])
		x['time_u'] = np.power(10,x['time_u'])
		total_sfr = []
		for i in range(len(x)):
			total_sfr.append((x['time_u'][i] - x['time_l'][i]) * x['SFR'][i])
		total_sfr = sum(total_sfr)

		time_temp = np.linspace(x['time_l'][0],x['time_u'][-1],10000)
		sfr = np.zeros_like(time_temp)
		for i in range(len(x)):
			sfr[np.where(np.logical_and(time_temp>=x['time_l'][i],time_temp<x['time_u'][i]))] = x['SFR'][i]
		self.sfr = np.interp(self.t,time_temp/1e9,sfr)
		self.sfr = np.divide(self.sfr * total_sfr, sum(self.sfr))[::-1]
		self.sfr /= 10000

