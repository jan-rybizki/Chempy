import matplotlib.pyplot as plt
import numpy as np

def gaussian(x,x0,xsig):
	'''
	function to calculate the gaussian probability (its normed to Pmax and given in log)
	
	INPUT:
	
	   x = where is the data point or parameter value
	
	   x0 = mu
	
	   xsig = sigma
	'''
	factor = 1. / (np.sqrt(xsig * xsig * 2. * np.pi))
	exponent = -np.divide((x - x0) * (x - x0),2 * xsig * xsig)
	return factor * np.exp(exponent)

def likelihood_evaluation(model_error, star_error_list, abundance_list, star_abundance_list):
	'''
	This function evaluates the Gaussian for the prediction/observation comparison and returns the resulting log likelihood. The model error and the observed error are added quadratically
	
	INPUT:

	   model_error = the error coming from the models side

	   star_error_list = the error coming from the observations

	   abundance_list = the predictions

	   star_abundance_list = the observations

	OUTPUT:

	   likelihood = the summed log likelihood
	'''

	error = np.sqrt(model_error * model_error + star_error_list * star_error_list)
	list_of_likelihoods = gaussian(star_abundance_list,abundance_list,error)
	log_likelihood_list = np.log(list_of_likelihoods)
	likelihood = sum(log_likelihood_list)
	return likelihood

def sample_stars(weight,selection,element1,element2,error1,error2,nsample):
	'''
	This function samples stars along a chemical evolution track properly taking into account the SFR and the selection function of the stellar population (e.g. red-clump stars). It can be used to produce mock observations which can be compared to survey data.

	INPUT

	   weight: The SFR of the model
	
	   selection: The age-distribution of a stellar population (e.g. red-clump stars). The time intervals need to be the same as for 'weight'.
	
	   element1 = the values of one element for the ISM of the model (same time-intervals as SFR)
	
	   element2 = the values of the other element for the ISM
	
	   error1 = the measurement error of that element
	
	   error2 = measurement error
	
	   nsample = number of stars that should be realized
	'''
	weight = np.cumsum(weight*selection)
	weight /= weight[-1]
	sample = np.random.random(nsample)
	sample = np.sort(sample)
	stars = np.zeros_like(weight)
	for i,item in enumerate(weight):
		if i == 0:
			count = len(sample[np.where(np.logical_and(sample>0.,sample<=item))])
			stars[i] = count
		else:
			count = len(sample[np.where(np.logical_and(sample>weight[i-1],sample<=item))])
			stars[i] = count
	sun_feh = []
	sun_mgfe = []
	for i in range(len(weight)):
		if stars[i] != 0:
			for j in range(int(stars[i])):
				sun_feh.append(element1[i])
				sun_mgfe.append(element2[i])
	sun_feh = np.array(sun_feh)
	sun_mgfe = np.array(sun_mgfe)
	perturbation = np.random.normal(0,error1,len(sun_feh))
	sun_feh += perturbation
	perturbation = np.random.normal(0,error2,len(sun_feh))
	sun_mgfe += perturbation
	return sun_feh,sun_mgfe

def gaussian_1d_log(x,x0,xsig):
	return -np.divide((x-x0)*(x-x0),2*xsig*xsig)

def yield_plot(name_string, yield_class, solar_class, element):
	'''
	This function plots [X/Fe] for the complete mass and metallicity range of a yield class.

	INPUT:
	
	   name_string = a string which is included in the saved file name
	
	   yield_class = a Chempy yield class (e.g 'Nomoto2013' from SN2_Feedback, see tutorial)
	
	   solar_class = a Chempy solar class (needed for normalisation)
	
	   element = the element for which to plot the yield table

	OUTPUT
	
	   the figure will be saved into the current directory
	'''
	elements = np.hstack(solar_class.all_elements)
	solar_fe_fraction = float(solar_class.fractions[np.where(elements == 'Fe')])
	solar_element_fraction = float(solar_class.fractions[np.where(elements == element)])
	plt.clf()
	fig = plt.figure(figsize=(13,8), dpi=100)
	ax = fig.add_subplot(111)

	ax.set_title('Yields of %s' %(name_string))
	ax.set_xlabel(r'metallicity in $\log_{10}\left(\mathrm{Z}/\mathrm{Z}_\odot\right)$')
	ax.set_ylabel('[%s/Fe]' %(element))

	for item in yield_class.metallicities:
		for j,jtem in enumerate(list(yield_class.table[item]['Mass'])):
			ejecta_fe = yield_class.table[item]['Fe'][j]
			ejecta_element = yield_class.table[item][element][j]
			#print "log Z, mass, [X/Fe]"
			#print np.log10(float(item)/solar_class.z), jtem,  np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			if item == 0:
				metallicity = np.log10(float(1e-7)/solar_class.z)
			else:
				metallicity = np.log10(float(item)/solar_class.z)
			alpha_enhancement = np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			mass = jtem
			ax.scatter(metallicity, alpha_enhancement, s=20, c=None, marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, edgecolors=None)
			ax.text(metallicity, alpha_enhancement, mass)

	plt.savefig('yields_%s.png' %(name_string),bbox_inches='tight')

def yield_comparison_plot(yield_name1, yield_name2, yield_class, yield_class2, solar_class, element):
	'''
	a function to plot a comparison between two yield sets. It is similar to 'yield_plot' only that a second yield set can be plotted.

	INPUT
	
	   yield_name1 = Name of the first yield class
	
	   yield_name2 = Name of the second yield class
	
	   yield_class = a Chempy yield class (e.g 'Nomoto2013' from SN2_Feedback, see tutorial)
	
	   yield_class2 = a Chempy yield class (e.g 'Nomoto2013' from SN2_Feedback, see tutorial)
	
	   solar_class = a Chempy solar class (needed for normalisation)
	
	   element = the element for which to plot the yield table

	OUTPUT
	
	   the figure will be saved into the current directory
	'''
	elements = np.hstack(solar_class.all_elements)
	solar_fe_fraction = float(solar_class.fractions[np.where(elements == 'Fe')])
	solar_element_fraction = float(solar_class.fractions[np.where(elements == element)])
	plt.clf()
	fig = plt.figure(figsize=(13,8), dpi=100)
	ax = fig.add_subplot(111)

	ax.set_title('Yields of %s in blue vs %s in red' %(yield_name1,yield_name2))
	ax.set_xlabel(r'metallicity in $\log_{10}\left(\mathrm{Z}/\mathrm{Z}_\odot\right)$')
	ax.set_ylabel('[%s/Fe]' %(element))

	for item in yield_class.metallicities:
		for j,jtem in enumerate(list(yield_class.table[item]['Mass'])):
			ejecta_fe = yield_class.table[item]['Fe'][j]
			ejecta_element = yield_class.table[item][element][j]
			#print "log Z, mass, [X/Fe]"
			#print np.log10(float(item)/solar_class.z), jtem,  np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			if item == 0:
				metallicity = np.log10(float(1e-7)/solar_class.z)
			else:
				metallicity = np.log10(float(item)/solar_class.z)
			alpha_enhancement = np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			mass = jtem
			ax.scatter(metallicity, alpha_enhancement, s=20*mass, c='b', marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=0.5, linewidths=None, verts=None, edgecolors=None)
			ax.annotate(xy = (metallicity, alpha_enhancement), s = mass, color = 'b')

	for item in yield_class2.metallicities:
		for j,jtem in enumerate(list(yield_class2.table[item]['Mass'])):
			ejecta_fe = yield_class2.table[item]['Fe'][j]
			ejecta_element = yield_class2.table[item][element][j]
			#print "log Z, mass, [X/Fe]"
			#print np.log10(float(item)/solar_class.z), jtem,  np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			if item == 0:
				metallicity = np.log10(float(1e-7)/solar_class.z)
			else:
				metallicity = np.log10(float(item)/solar_class.z)
			alpha_enhancement = np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			mass = jtem
			ax.scatter(metallicity, alpha_enhancement, s=20*mass, c='r', marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=0.5, linewidths=None, verts=None, edgecolors=None)
			ax.annotate(xy = (metallicity, alpha_enhancement), s = mass, color = 'r')

	plt.savefig('yields_comparison_%s_vs_%s_for_%s.png' %(yield_name1,yield_name2,element),bbox_inches='tight')

def fractional_yield_comparison_plot(yield_name1, yield_name2, yield_class, yield_class2, solar_class, element):
	'''
	a function to plot a comparison between the fractional yield of two yield sets. The fractional yield is the mass fraction of the star that is expelled as the specific element. Depending on the yield class it will be net or total yield.

	INPUT
	
	   yield_name1 = Name of the first yield class
	
	   yield_name2 = Name of the second yield class
	
	   yield_class = a Chempy yield class (e.g 'Nomoto2013' from SN2_Feedback, see tutorial)
	
	   yield_class2 = a Chempy yield class (e.g 'Nomoto2013' from SN2_Feedback, see tutorial)
	
	   solar_class = a Chempy solar class (needed for normalisation)
	
	   element = the element for which to plot the yield table

	OUTPUT
	
	   the figure will be saved into the current directory
	'''

	plt.clf()
	fig = plt.figure(figsize=(13,8), dpi=100)
	ax = fig.add_subplot(111)

	ax.set_title('Yields of %s in blue vs %s in red' %(yield_name1,yield_name2))
	ax.set_xlabel(r'metallicity in $\log_{10}\left(\mathrm{Z}/\mathrm{Z}_\odot\right)$')
	ax.set_ylabel(r'$\log_{10}$(fractional %s yield)' %(element))

	for item in yield_class.metallicities:
		for j,jtem in enumerate(list(yield_class.table[item]['Mass'])):
			ejecta_element = yield_class.table[item][element][j]
			#print "log Z, mass, [X/Fe]"
			#print np.log10(float(item)/solar_class.z), jtem,  np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			if item == 0:
				metallicity = np.log10(float(1e-7)/solar_class.z)
			else:
				metallicity = np.log10(float(item)/solar_class.z)
			fractional_feedback = np.log10(ejecta_element)
			mass = jtem
			ax.scatter(metallicity, fractional_feedback, s=20*mass, c='b', marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=0.5, linewidths=None, verts=None, edgecolors=None)
			ax.annotate(xy = (metallicity, fractional_feedback), s = mass, color = 'b')

	for item in yield_class2.metallicities:
		for j,jtem in enumerate(list(yield_class2.table[item]['Mass'])):
			ejecta_element = yield_class2.table[item][element][j]
			#print "log Z, mass, [X/Fe]"
			#print np.log10(float(item)/solar_class.z), jtem,  np.log10(ejecta_element/solar_element_fraction) - np.log10(ejecta_fe/solar_fe_fraction)
			if item == 0:
				metallicity = np.log10(float(1e-7)/solar_class.z)
			else:
				metallicity = np.log10(float(item)/solar_class.z)
			fractional_feedback = np.log10(ejecta_element)
			mass = jtem
			ax.scatter(metallicity, fractional_feedback, s=20*mass, c='r', marker=u'o', cmap=None, norm=None, vmin=None, vmax=None, alpha=0.5, linewidths=None, verts=None, edgecolors=None)
			ax.annotate(xy = (metallicity, fractional_feedback), s = mass, color = 'r')

	plt.savefig('fractional_yields_comparison_%s_vs_%s_for_%s.png' %(yield_name1,yield_name2,element),bbox_inches='tight')


def elements_plot(name_string,agb, sn2, sn1a,elements_to_trace, all_elements,max_entry):
	'''
	This function plots the available elements for specific yield sets.

	INPUT:
	
	   name_string = a string that will be added to the file name
	   
	   agb = agb yield class
	
	   sn2 = sn2 yield class
	
	   sn1a = sn1a yield class
	
	   elements_to_trace = which elements do we want to follow (a list)
	
	   all_elements = Symbols of all elements (available in the solar abundances class)
	
	   max_entry = until which element number the figure should be plotted

	OUTPUT:
	
	   a figure in the current directory
	'''
	plt.clf()
	fig = plt.figure(figsize=(max_entry,10.27), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_title('Elements in yields and processes')
	ax.set_xlabel('element number')
	ax.set_ylabel('groups of elements')


	apogee = ['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Ni']
	s_process = ['Sr','Y','Zr','Nb','Mo','Sn','Ba','La','Ce','Nd','W','Pb']
	r_process = ['Ge','Ru','Rh','Pd','Ag','Pr','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt','Au','Bi','Th','U']
	alpha = ['C','N','O','Ne','Mg','Si','S','Ar','Ca','Ti']
	iron = ['V','Cr','Mn','Fe','Co','Ni','Zn']
	BB = ['H','He','Li','Be']

	#ax.axis([-3, 118, -3.5, 8.5])
	ax.set_xlim((-3,max_entry))
	ax.set_ylim((-3.5,8.5))
	ax.text(-2, 3, 'elements', fontsize=15, clip_on=True)
	ax.text(-2, 4, 'We trace', fontsize=15)
	ax.text(-2, 8, 'SN Ia', fontsize=15)
	ax.text(-2, 6, 'CC-SN', fontsize=15)
	ax.text(-2, 7, 'AGB', fontsize=15)
	ax.text(-2, 5, 'mutual', fontsize=15)
	ax.text(-2, -3, r'$\alpha$-element', fontsize=15)
	ax.text(-2, 1, 'iron-peak', fontsize=15)
	ax.text(-2, 0, 's-process', fontsize=15)
	ax.text(-2, -1, 'r-process', fontsize=15)
	ax.text(-2, -2, 'Big Bang', fontsize=15)
	ax.text(-2, 2, 'Apogee', fontsize=15)
	for i, item in enumerate(all_elements['Symbol']):
		ax.text(all_elements['Number'][i], 3, item, fontsize=15,bbox={'facecolor':'red', 'alpha':0.5, 'pad':10}, clip_on=True)
		if item in elements_to_trace:
			ax.text(all_elements['Number'][i], 4, "X", fontsize=15, clip_on=True)		
		if item in sn1a:
			ax.text(all_elements['Number'][i], 8, "X", fontsize=15, clip_on=True)	
		if item in sn2:
			ax.text(all_elements['Number'][i], 6, "X", fontsize=15, clip_on=True)	
		if item in agb:
			ax.text(all_elements['Number'][i], 7, "X", fontsize=15, clip_on=True)	
		if item in agb and item in sn2 and item in sn1a:
			ax.text(all_elements['Number'][i], 5, "X", fontsize=15, clip_on=True)	
		if item in BB:
			ax.text(all_elements['Number'][i], -2, "X", fontsize=15, clip_on=True)
		if item in iron:
			ax.text(all_elements['Number'][i], 1, "X", fontsize=15, clip_on=True)	
		if item in alpha:
			ax.text(all_elements['Number'][i], -3, "X", fontsize=15, clip_on=True)	
		if item in s_process:
			ax.text(all_elements['Number'][i], 0, "X", fontsize=15, clip_on=True)	
		if item in r_process:
			ax.text(all_elements['Number'][i], -1, "X", fontsize=15, clip_on=True)	
		if item in apogee:
			ax.text(all_elements['Number'][i], 2, "X", fontsize=15, clip_on=True)		

	plt.savefig('elements_%s.png' %(name_string),bbox_inches='tight')
	return [0.]

def cosmic_abundance_standard(summary_pdf,name_string,abundances,cube,elements_to_trace,solar,number_of_models_overplotted,produce_mock_data,use_mock_data,error_inflation):
	'''
	This is a likelihood function for the B-stars (a proxy for the present-day ISM) with data from nieva przybilla 2012 cosmic abundance standard paper.

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	
	   abundances = abundances of the IMS (can be calculated from cube)
	
	   cube = the ISM mass fractions per element (A Chempy class containing the model evolution)
	
	   elements_to_trace = Which elements should be used for the analysis
	
	   solar = solar abundance class
	
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	   produce_mock_data = should the predictions be saved with an error added to them (default=False)
	
	   use_mock_data = instead of the real data use the formerly produced mock data (default = False)
	
	   error_inflation = a factor with which the mock data should be perturbed with the observational data (default = 1.0)

	OUTPUT:
	   
	   probabilities = each elements likelihood in a list 
	
	   model_abundances = each elements model prediction as a list
	
	   elements_list = the names of the elements in a list
	'''
	log_abundances_cosmic = np.array([10.99,8.76,7.56,7.5,7.52,8.09,8.33,7.79])#
	error_abundances = np.array([0.01,0.05,0.05,0.05,0.03,0.05,0.04,0.04])# original
	elements= np.array(['He','O','Mg','Si','Fe','Ne','C','N'])
	probabilities = []
	probabilities_max = []
	elements_list = []
	error_list = []
	cas_abundances = []
	model_abundances = []

	for i,item in enumerate(elements):
		if item in elements_to_trace:
			elements_list.append(item)
			cas_abundances.append(log_abundances_cosmic[np.where(elements==item)]-solar['photospheric'][np.where(solar['Symbol']==item)])
			model_abundances.append(abundances[item][-1])
			error_list.append(error_abundances[np.where(elements==item)])


	if produce_mock_data:
		mock_abundance_list = list(np.random.normal(loc = list(np.hstack(model_abundances)),scale = list(error_inflation*np.hstack(error_list))))
		np.save('mock_data_temp/cas_abundances' ,mock_abundance_list)

	if use_mock_data:
		error_list = list(np.hstack(error_list)*error_inflation)
		cas_abundances = np.load('mock_data_temp/cas_abundances.npy')


	for i,item in enumerate(elements_list):
		probabilities.append(float(gaussian_1d_log(model_abundances[i],cas_abundances[i],error_list[i])))
	probability = np.sum(probabilities)


	if number_of_models_overplotted >1:
		if os.path.isfile('output/comparison/cas.npy'):
			old = np.load('output/comparison/cas.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array(model_abundances))
		np.save('output/comparison/cas',old)
		
		if os.path.isfile('output/comparison/cas_likelihood.npy'):
			old_likelihood = np.load('output/comparison/cas_likelihood.npy')
			old_likelihood = list(old_likelihood)
		else:
			old_likelihood = []
		old_likelihood.append(np.array(probabilities))
		np.save('output/comparison/cas_likelihood',old_likelihood)

	if summary_pdf:
		plt.clf()
		fig = plt.figure(figsize=(11.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.errorbar(np.arange(len(elements_list)),cas_abundances,xerr=None,yerr=error_list,mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label='CAS')
		plt.plot(np.arange(len(elements_list)),model_abundances,label='model after %.1f Gyr' %(cube['time'][-1]),linestyle='-')

		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(np.arange(len(elements_list)),np.array(item),linestyle='-', color = 'g', alpha = 0.2)


		for i in range(len(elements_list)):
			plt.annotate(xy=(i,-0.1),s= '%.2f' %(probabilities[i]),size=12)
		plt.grid("on")
		elements_list1 = ['[%s/H]' %(item) for item in elements_list]
		plt.ylim((-0.5,0.5))
		plt.xticks(np.arange(len(elements_list1)), elements_list1)
		plt.ylabel("abundance relative to solar in dex")
		plt.xlabel("Element")
		plt.title('ln(probability) of fulfilling cas= %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('cas_%s.png' %(name_string))
	return probabilities, model_abundances, elements_list

def sol_norm(summary_pdf,name_string,abundances,cube,elements_to_trace, element_names, sol_table, number_of_models_overplotted,produce_mock_data,use_mock_data,error_inflation):
	'''
	This is a likelihood function for solar abundances compared to the model ISM abundances from 4.5Gyr ago.

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	
	   abundances = abundances of the IMS (can be calculated from cube)
	   
	   cube = the ISM mass fractions per element (A Chempy class containing the model evolution)
	
	   elements_to_trace = which elements are tracked by Chempy
	   
	   element_names = which elements should be used for the likelihood

	   sol_table = solar abundance class
	
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	   produce_mock_data = should the predictions be saved with an error added to them (default=False)
	
	   use_mock_data = instead of the real data use the formerly produced mock data (default = False)
	
	   error_inflation = a factor with which the mock data should be perturbed with the observational data (default = 1.0)

	OUTPUT:
	
	   probabilities = each elements likelihood in a list 
	
	   model_abundances = each elements model prediction as a list
	
	   elements_list = the names of the elements in a list
	'''


	'''
	solar abundances of the sun from the photospheric abundances of the basic_solar.table
	compared to the model abundances after 7.5Gyr. If the sun travelled in the galaxy radial abundance gradiants need to be added.
	'''
	elements_to_trace = element_names

	if 'C+N' in element_names:
		new_array = np.log10( np.power(10,abundances['C']) + np.power(10,abundances['N']))
		abundances = append_fields(abundances,'C+N',new_array)


	time_sun = cube['time'][-1] - 4.5
	cut = [np.where(np.abs(cube['time']-time_sun)==np.min(np.abs(cube['time']-time_sun)))]
	if len(cut[0][0]) != 1:
		cut = cut[0][0][0]
	time_model = cube['time'][cut]
	probabilities = []
	abundance_list = []
	error_list = []
	sun_list = []
	for i,item in enumerate(elements_to_trace):
		abundance_list.append(float(abundances[item][cut]))
		error = sol_table['error'][np.where(sol_table['Symbol']==item)] + 0.01 # added because of uncertainty in diffusion correction (Asplund2009)
		error_list.append(error)
		if item != "C+N":
			if item == 'He':
				sun_list.append(0.05)
			else:
				sun_list.append(0.04) ## add 0.04dex to get protosolar abundances (Asplund 2009)
		else:
			sun_list.append(np.log10(2.))

	if produce_mock_data:
		mock_abundance_list = list(np.random.normal(loc = list(np.hstack(abundance_list)),scale = list(error_inflation*np.hstack(error_list))))
		np.save('mock_data_temp/solar_abundances' ,mock_abundance_list)

	if use_mock_data:
		error_list = list(np.hstack(error_list)*error_inflation)
		sun_list = np.load('mock_data_temp/solar_abundances.npy')

	for i,item in enumerate(elements_to_trace):
		probabilities.append(float(gaussian_1d_log(abundance_list[i],sun_list[i],error_list[i])))

	probability = np.sum(probabilities)

	if number_of_models_overplotted > 1:
		if os.path.isfile('output/comparison/sun.npy'):
			old = np.load('output/comparison/sun.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array(abundance_list))
		np.save('output/comparison/sun',old)
		if os.path.isfile('output/comparison/sun_likelihood.npy'):
			old_likelihood = np.load('output/comparison/sun_likelihood.npy')
			old_likelihood = list(old_likelihood)
		else:
			old_likelihood = []
		old_likelihood.append(np.array(probabilities))
		np.save('output/comparison/sun_likelihood',old_likelihood)

	if summary_pdf:
		text_size = 12
		plt.rc('font', family='serif',size = text_size)
		plt.rc('xtick', labelsize=text_size)
		plt.rc('ytick', labelsize=text_size)
		plt.rc('axes', labelsize=text_size, lw=2.0)
		plt.rc('lines', linewidth = 2)
		plt.rcParams['ytick.major.pad']='8'
		plt.clf()
		fig = plt.figure(figsize=(30.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.errorbar(np.arange(len(elements_to_trace)),sun_list,xerr=None,yerr=error_list,linestyle = '',mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label='solar')
		
		plt.plot(np.arange(len(elements_to_trace)),np.array(abundance_list),label='model after %.2f Gyr' %(time_model),linestyle='-')
		
		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(np.arange(len(elements_to_trace)),np.array(item),linestyle='-', color = 'g', alpha = 0.2)


		for i in range(len(elements_to_trace)):
			plt.annotate(xy=(i,-0.4),s= '%.2f' %(probabilities[i]))
		plt.grid("on")
		plt.ylim((-0.5,0.5))
		elements_to_trace = ['[%s/H]' %(item) for item in elements_to_trace]
		plt.xticks(np.arange(len(elements_to_trace)), elements_to_trace)
		plt.ylabel("abundance relative to solar in dex")
		plt.xlabel("Element")
		plt.title('joint probability of agreeing with the sun (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('sol_norm_%s.png' %(name_string))
	return probabilities, abundance_list, elements_to_trace



def arcturus(summary_pdf,name_string,abundances,cube,elements_to_trace, element_names, sol_table, number_of_models_overplotted, arcturus_age,produce_mock_data,use_mock_data,error_inflation):
	'''
	This is a likelihood function for the arcturus abundances compared to the model ISM abundances from some time ago (the age of arcturus which needs to be specified). The abundances are taken from Ramirez+ 2011 and all elements except for Fe are given in [X/Fe]

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	   
	   abundances = abundances of the IMS (can be calculated from cube)
	
	   cube = the ISM mass fractions per element (A Chempy class containing the model evolution)
	
	   elements_to_trace = Which elements should be used for the analysis
	
	   sol_table = solar abundance class
	
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	   produce_mock_data = should the predictions be saved with an error added to them (default = False)
	
	   use_mock_data = instead of the real data use the formerly produced mock data (default = False)
	   
	   error_inflation = a factor with which the mock data should be perturbed with the observational data (default = 1.0)

	OUTPUT:
	
	   probabilities = each elements likelihood in a list 
	
	   model_abundances = each elements model prediction as a list
	   
	   elements_list = the names of the elements in a list
	'''
	#### Arcturus abundances taken from Worley+ 2009 MNRAS all elements except Fe are given [X/Fe]
	#element_list_arcturus = ['Fe', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Sc', 'Ti', 'Zn', 'Y', 'Zr', 'Ba', 'La', 'Nd', 'Eu']
	#element_abundances_arcturus = [-0.60,0.57,0.15,0.34,0.25,0.24,0.19,0.24,0.34,-0.04,0.10,0.06,-0.19,0.04,0.10,0.36]
	#element_errors_arcturus = [0.11,0.02,0.04,0.15,0.07,0.14,0.06,0.01,0.11,0.09,0.17,0.08,0.08,0.08,0.07,0.04]
	#### Arcturus abundances and Age taken from Ramirez+ 2011 all elements except Fe given in [X/Fe]
	element_list_arcturus = np.array(['Fe', 'O', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Co', 'Ni', 'Zn'])#, 'C'
	element_abundances_arcturus = np.array([-0.52, 0.50, 0.11, 0.37, 0.34, 0.33, 0.20, 0.11, 0.19, 0.24, 0.20, -0.05, -0.21, 0.09, 0.06, 0.22])#, 0.43
	element_errors_arcturus = np.array([0.04, 0.03, 0.03, 0.03, 0.03, 0.04, 0.07, 0.04, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.03, 0.06])#,0.07

	assert len(element_list_arcturus) == len(element_abundances_arcturus) == len(element_errors_arcturus)
	elements_to_trace = element_names

	time_arcturus = cube['time'][-1] - arcturus_age # 7.1 +1.5 -1.2
	cut = [np.where(np.abs(cube['time']-time_arcturus)==np.min(np.abs(cube['time']-time_arcturus)))]
	if len(cut[0][0]) != 1:
		cut = cut[0][0][0]
	time_model = cube['time'][cut]

	abundance_list = []
	error_list = []
	arcturus_list = []
	elements_in_common = []
	for i,item in enumerate(elements_to_trace):
		if item in element_list_arcturus:
			elements_in_common.append(item)
			if item == 'Fe':
				abundance_list.append(float(abundances[item][cut]))
			else:
				abundance_list.append(float(abundances[item][cut]-abundances['Fe'][cut]))
			error_list.append(element_errors_arcturus[np.where(element_list_arcturus == item)])
			arcturus_list.append(element_abundances_arcturus[np.where(element_list_arcturus == item)])

	if produce_mock_data:
		mock_abundance_list = list(np.random.normal(loc = list(np.hstack(abundance_list)),scale = list(error_inflation*np.hstack(error_list))))
		np.save('mock_data_temp/arcturus_abundances' ,mock_abundance_list)

	if use_mock_data:
		error_list = list(np.hstack(error_list)*error_inflation)
		arcturus_list = np.load('mock_data_temp/arcturus_abundances.npy')


	probabilities = []
	for i,item in enumerate(elements_in_common):
		probabilities.append(float(gaussian_1d_log(abundance_list[i],arcturus_list[i],error_list[i])))

	probability = np.sum(probabilities)

	if number_of_models_overplotted > 1:
		if os.path.isfile('output/comparison/arc.npy'):
			old = np.load('output/comparison/arc.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array(abundance_list))
		np.save('output/comparison/arc',old)

		if os.path.isfile('output/comparison/arc_likelihood.npy'):
			old_likelihood = np.load('output/comparison/arc_likelihood.npy')
			old_likelihood = list(old_likelihood)
		else:
			old_likelihood = []
		old_likelihood.append(np.array(probabilities))
		np.save('output/comparison/arc_likelihood',old_likelihood)

	if summary_pdf:
		text_size = 8
		plt.rc('font', family='serif',size = text_size)
		plt.rc('xtick', labelsize=text_size)
		plt.rc('ytick', labelsize=text_size)
		plt.rc('axes', labelsize=text_size, lw=1.0)
		plt.rc('lines', linewidth = 1)
		plt.rcParams['ytick.major.pad']='8'
		plt.clf()
		fig = plt.figure(figsize=(10.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.errorbar(np.arange(len(elements_in_common)),arcturus_list,xerr=None,yerr=error_list,mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label='arcturus')
		plt.plot(np.arange(len(elements_in_common)),np.array(abundance_list),label='model after %.2f Gyr' %(time_model),linestyle='-')

		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(np.arange(len(elements_in_common)),np.array(item),linestyle='-', color = 'g', alpha = 0.2)
		for i in range(len(elements_in_common)):
			plt.annotate(xy=(i,-0.4),s= '%.2f' %(probabilities[i]))
		plt.grid("on")
		plt.ylim((-0.5,0.5))
		element_labels = ['[%s/Fe]' %(item) for item in elements_in_common]
		element_labels = np.array(element_labels)
		element_labels[np.where(element_labels == '[Fe/Fe]')] = '[Fe/H]'
		plt.xticks(np.arange(len(elements_in_common)), element_labels)
		plt.ylabel("abundance relative to solar in dex")
		plt.xlabel("Element")
		plt.title('joint probability of agreeing with the sun (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('arcturus_%s.png' %(name_string))

	return probabilities, abundance_list, elements_in_common

def gas_reservoir_metallicity(summary_pdf,name_string,abundances,cube,elements_to_trace,gas_reservoir,number_of_models_overplotted,produce_mock_data,use_mock_data, error_inflation,solZ):
	'''
	This is a likelihood function for the present-day metallicity of the corona gas. Data is taken from the Smith cloud

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	
	   abundances = abundances of the IMS (can be calculated from cube)
	
	   cube = the ISM mass fractions per element (A Chempy class containing the model evolution)
	
	   elements_to_trace = Which elements should be used for the analysis
	
	   gas_reservoir = the gas_reservoir class containing the gas_reservoir evolution
	
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	   produce_mock_data = should the predictions be saved with an error added to them (default=False)
	
	   use_mock_data = instead of the real data use the formerly produced mock data (default = False)
	   
	   error_inflation = a factor with which the mock data should be perturbed with the observational data (default = 1.0)
	
	   solZ = solar metallicity

	OUTPUT:
	
	   probabilities = each elements likelihood in a list 
	
	   model_abundances = each elements model prediction as a list
	
	   elements_list = the names of the elements in a list
	'''
	### data from Smith cloud 
	log10_metallicity_at_end = -0.28
	log10_std = 0.14
	if produce_mock_data:
		mock_data = np.random.normal(loc = np.log10(gas_reservoir['Z'][-1]/solZ), scale = log10_std*error_inflation)
		np.save('mock_data_temp/metallicity_at_end' , mock_data)

	if use_mock_data:
		log10_std *= error_inflation
		log10_metallicity_at_end = np.load('mock_data_temp/metallicity_at_end.npy')

	probability = gaussian_1d_log(np.log10(gas_reservoir['Z'][-1]/solZ),log10_metallicity_at_end,log10_std)


	if number_of_models_overplotted > 1:
		if os.path.isfile('output/comparison/gas_reservoir.npy'):
			old = np.load('output/comparison/gas_reservoir.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array(np.log10(gas_reservoir['Z']/solZ)))
		np.save('output/comparison/gas_reservoir',old)

		if os.path.isfile('output/comparison/gas_reservoir_likelihood.npy'):
			old_likelihood = np.load('output/comparison/gas_reservoir_likelihood.npy')
			old_likelihood = list(old_likelihood)
		else:
			old_likelihood = []
		old_likelihood.append(np.array(probability))
		np.save('output/comparison/gas_reservoir_likelihood',old_likelihood)
		
		if os.path.isfile('output/comparison/gas_metallicity.npy'):
			old1 = np.load('output/comparison/gas_metallicity.npy')
			old1 = list(old1)
		else:
			old1 = []
		old1.append(np.array(np.log10(cube['Z']/solZ)))
		np.save('output/comparison/gas_metallicity',old1)	

	if summary_pdf:
		plt.clf()
		text_size = 8
		plt.rc('font', family='serif',size = text_size)
		plt.rc('xtick', labelsize=text_size)
		plt.rc('ytick', labelsize=text_size)
		plt.rc('axes', labelsize=text_size, lw=1.0)
		plt.rc('lines', linewidth = 1)
		plt.rcParams['ytick.major.pad']='8'
		fig = plt.figure(figsize=(11.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.plot(gas_reservoir['time'],np.log10(gas_reservoir['Z']/solZ),label = "Z_gas_reservoir %.4f" %(np.log10(gas_reservoir['Z'][-1]/solZ)))

		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(gas_reservoir['time'],np.array(item),linestyle='-', color = 'b', alpha = 0.2)
			for item in old1:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'r', alpha = 0.2)

		plt.errorbar(gas_reservoir['time'][-1],log10_metallicity_at_end,xerr=None,yerr=log10_std,mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label='Smith cloud')
		plt.plot(cube['time'],np.log10(cube['Z']/solZ),label = "Z_gas %.4f" %(np.log10(cube['Z'][-1]/solZ)))
		plt.grid("on")
		plt.ylabel("log 10 metallicity in Z (mass fraction of gas in metals)")
		plt.xlabel("time in Gyr")
		plt.title('ln(probability) gas content (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('gas_reservoir_%s.png' %(name_string))

	return [probability],[gas_reservoir['Z'][-1]/solZ],['Corona metallicity']



def ratio_function(summary_pdf,name_string,abundances,cube,elements_to_trace,gas_reservoir, number_of_models_overplotted,produce_mock_data,use_mock_data, error_inflation):
	'''
	This is a likelihood function for the present-day SN-ratio. Data is approximated from Manucci+2005.

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	
	   abundances = abundances of the IMS (can be calculated from cube)
	
	   cube = the ISM mass fractions per element (A Chempy class containing the model evolution)
	
	   elements_to_trace = Which elements should be used for the analysis
	
	   gas_reservoir = the gas_reservoir class containing the gas_reservoir evolution
	
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	   produce_mock_data = should the predictions be saved with an error added to them (default=False)
	
	   use_mock_data = instead of the real data use the formerly produced mock data (default = False)
	
	   error_inflation = a factor with which the mock data should be perturbed with the observational data (default = 1.0)

	OUTPUT:
	
	   probabilities = each elements likelihood in a list 
	
	   model_abundances = each elements model prediction as a list
	
	   elements_list = the names of the elements in a list
	'''
	log10_ratio_at_end = 0.7
	log10_std = 0.37

	if produce_mock_data:
		mock_data = np.random.normal(loc = np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])[-1]),scale = error_inflation*log10_std)
		np.save('mock_data_temp/sn_ratio' ,mock_data)

	if use_mock_data:
		log10_std *= error_inflation
		log10_ratio_at_end = np.load('mock_data_temp/sn_ratio.npy')

	probability = gaussian_1d_log( np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])[-1]),log10_ratio_at_end,log10_std)
	if number_of_models_overplotted > 1:
		if os.path.isfile('output/comparison/ratio.npy'):
			old = np.load('output/comparison/ratio.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array( np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:]))))
		np.save('output/comparison/ratio',old)		

		if os.path.isfile('output/comparison/ratio_likelihood.npy'):
			old_likelihood = np.load('output/comparison/ratio_likelihood.npy')
			old_likelihood = list(old_likelihood)
		else:
			old_likelihood = []
		old_likelihood.append(np.array(probability ))
		np.save('output/comparison/ratio_likelihood',old_likelihood)		

		if os.path.isfile('output/comparison/number_sn2.npy'):
			old1 = np.load('output/comparison/number_sn2.npy')
			old1 = list(old1)
		else:
			old1 = []
		old1.append(cube['sn2'])
		np.save('output/comparison/number_sn2',old1)	

		if os.path.isfile('output/comparison/number_sn1a.npy'):
			old2 = np.load('output/comparison/number_sn1a.npy')
			old2 = list(old2)
		else:
			old2 = []
		old2.append(cube['sn1a'])
		np.save('output/comparison/number_sn1a',old2)	


	if summary_pdf:
		plt.clf()
		fig = plt.figure(figsize=(11.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.plot( cube['time'][2:],np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])), label = 'sn2/sn1a of the model')
		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(cube['time'][2:],np.array(item),linestyle='-', color = 'b', alpha = 0.2)


		plt.annotate(xy=(cube['time'][-1],np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])[-1])),s= '%.2f' %(np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])[-1])))
		plt.errorbar(cube['time'][-1],log10_ratio_at_end,xerr=None,yerr=log10_std,mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label='Prantzos+11')
		plt.grid("on")
		plt.ylabel("sn2/sn1a")
		plt.xlabel("time in Gyr")
		plt.title('ln(probability) of sn ratio (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('sn_ratio_%s.png' %(name_string))


		plt.clf()
	return [probability],[np.log10(np.divide(np.diff(cube['sn2'])[1:],np.diff(cube['sn1a'])[1:])[-1])],['SN-ratio']

def star_function(summary_pdf,name_string,abundances,cube,elements_to_trace,gas_reservoir,number_of_models_overplotted):
	'''
	!!! Should only be used for plotting purposes

	A likelihood function for the stellar surface mass density. It is not updated. Should only be used for plotting purposes as it shows the infall and SFR of a Chempy model
	'''
	# data from... Just Jahreiss 2010
	stars_at_end = 28.
	std = 2.
	dt = cube['time'][1] - cube['time'][0]
	probability = gaussian_1d_log(cube['stars'][-1],stars_at_end,std)
	if number_of_models_overplotted > 1:
		if os.path.isfile('output/comparison/sfr.npy'):
			old = np.load('output/comparison/sfr.npy')
			old = list(old)
		else:
			old = []
		old.append(np.array(cube['sfr']*(1./dt)))
		np.save('output/comparison/sfr',old)
		if os.path.isfile('output/comparison/infall.npy'):
			old1 = np.load('output/comparison/infall.npy')
			old1 = list(old1)
		else:
			old1 = []
		old1.append(np.array(cube['infall']*(1./dt)))
		np.save('output/comparison/infall',old1)
		
		if os.path.isfile('output/comparison/gas_mass.npy'):
			old2 = np.load('output/comparison/gas_mass.npy')
			old2 = list(old2)
		else:
			old2 = []
		old2.append(np.array(cube['gas']))
		np.save('output/comparison/gas_mass',old2)

		if os.path.isfile('output/comparison/star_mass.npy'):
			old3 = np.load('output/comparison/star_mass.npy')
			old3 = list(old3)
		else:
			old3 = []
		old3.append(np.array(cube['stars']))
		np.save('output/comparison/star_mass',old3)

		if os.path.isfile('output/comparison/remnant_mass.npy'):
			old4 = np.load('output/comparison/remnant_mass.npy')
			old4 = list(old4)
		else:
			old4 = []
		old4.append(np.array(cube['mass_in_remnants']))
		np.save('output/comparison/remnant_mass',old4)

		if os.path.isfile('output/comparison/corona_mass.npy'):
			old5 = np.load('output/comparison/corona_mass.npy')
			old5 = list(old5)
		else:
			old5 = []
		old5.append(np.array(gas_reservoir['gas']))
		np.save('output/comparison/corona_mass',old5)

		if os.path.isfile('output/comparison/feedback_mass.npy'):
			old6 = np.load('output/comparison/feedback_mass.npy')
			old6 = list(old6)
		else:
			old6 = []
		old6.append(np.array(cube['feedback']*(1./dt)))
		np.save('output/comparison/feedback_mass',old6)

	if summary_pdf:
		plt.clf()
		fig = plt.figure(figsize=(11.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.plot(cube['time'],cube['gas'],label = "gas %.2f" %(cube['gas'][-1]), color = 'b')
		plt.plot(cube['time'],cube['stars'],label = "stars (only thin disc including remnants) %.2f" %(cube['stars'][-1]), color = 'r')
		plt.plot(cube['time'],cube['mass_in_remnants'],label = "remnants %.2f" %(cube['mass_in_remnants'][-1]), color = 'k')
		plt.plot(cube['time'],gas_reservoir['gas'], label = 'corona %.2f'%(gas_reservoir['gas'][-1]), color = 'y')
		if number_of_models_overplotted > 1:
			for item in old2:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'b', alpha = 0.2)
			for item in old3:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'r', alpha = 0.2)
			for item in old4:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'g', alpha = 0.2)
			for item in old5:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'y', alpha = 0.2)
		plt.grid("on")
		plt.yscale('log')
		plt.ylabel(r"M$_\odot$")
		plt.xlabel("time in Gyr")
		plt.title('ln(probability) star content (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('stars_%s.png' %(name_string))

		plt.clf()
		fig = plt.figure(figsize=(11.69,8.27), dpi=100)
		ax = fig.add_subplot(111)
		plt.plot(cube['time'],cube['infall']*(1./dt),linestyle='-', color = 'r',label = "infall %.2f" %(sum(cube['infall'])))
		plt.plot(cube['time'],cube['sfr']*(1./dt),linestyle='-', color = 'b',label = "sfr %.2f" %(sum(cube['sfr'])))
		
		if number_of_models_overplotted > 1:
			for item in old:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'b', alpha = 0.2)
			for item in old1:
				plt.plot(cube['time'],np.array(item),linestyle='-', color = 'r', alpha = 0.2)


		plt.grid("on")
		plt.ylabel(r"M$_\odot$Gyr$^{-1}$")
		plt.xlabel("time in Gyr")
		plt.title('ln(probability) star content (normed to pmax) = %.2f' %(probability))	
		plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
		plt.savefig('infall_%s.png' %(name_string))
	return [0.]


def plot_processes(summary_pdf,name_string,sn2_cube,sn1a_cube,agb_cube,elements,cube1,number_of_models_overplotted):
	'''
	This is a plotting routine showing the different nucleosynthetic contributions to the individual elements.

	INPUT:

	   summary_pdf = boolean, should a pdf be created?
	
	   name_string = string to be added in the saved file name
	
	   sn2_cube = the sn2_feeback class
	
	   sn1a_cube = the sn1a_feeback class
	
	   agb_cube = the agb feedback class
	   
	   elements = which elements should be plotted
	
	   cube1 = the ISM mass fractions per element (A Chempy class containing the model evolution)
	   
	   number_of_models_overplotted = default is 1, if more the results will be saved and in the last iteration all former models will be plotted at once
	
	OUTPUT:
	
	   A plotfile in the current directory
	'''

	probability = 0

	sn2 = []
	agb = []
	sn1a= []
	for i,item in enumerate(elements):
		if item == 'C+N':
			sn2_temp = np.sum(sn2_cube['C'])
			sn2_temp += np.sum(sn2_cube['N'])
			sn2.append(sn2_temp)
			sn1a_temp = np.sum(sn1a_cube['C'])
			sn1a_temp += np.sum(sn1a_cube['N'])
			sn1a.append(sn1a_temp)	
			agb_temp = np.sum(agb_cube['C'])
			agb_temp += np.sum(agb_cube['N'])
			agb.append(agb_temp)				
		else:
			sn2.append(np.sum(sn2_cube[item]))
			sn1a.append(np.sum(sn1a_cube[item]))
			agb.append(np.sum(agb_cube[item]))
	sn2 = np.array(sn2)
	agb = np.array(agb)
	sn1a = np.array(sn1a)
	
	total_feedback = sn2 + sn1a + agb
	
	
	all_4 = np.vstack((sn2,sn1a,agb,total_feedback))
	if number_of_models_overplotted > 1:
		np.save('output/comparison/elements', elements)
		if os.path.isfile('output/comparison/temp_default.npy'):
			old = np.load('output/comparison/temp_default.npy')
			all_4 = np.dstack((old,all_4))
			np.save('output/comparison/temp_default', all_4)
		else:
			np.save('output/comparison/temp_default', all_4)

	if number_of_models_overplotted > 1:
		medians = np.median(all_4, axis = 2)
		stds = np.std(all_4,axis = 2)
	else:
		medians = all_4
		stds = np.zeros_like(all_4)

	if summary_pdf:
		if number_of_models_overplotted == 1:
			plt.clf()
			fig ,ax1 = plt.subplots( figsize = (20,10))
			l1 = ax1.bar(np.arange(len(elements)),np.divide(sn2,total_feedback), color = 'b' ,label='sn2', width = 0.2)
			l2 = ax1.bar(np.arange(len(elements))+0.25,np.divide(sn1a,total_feedback), color = 'g' ,label='sn1a', width = 0.2)
			l3 = ax1.bar(np.arange(len(elements))+0.5,np.divide(agb,total_feedback), color = 'r' ,label='agb', width = 0.2)
			ax1.set_ylim((0,1))
			plt.xticks(np.arange(len(elements))+0.4, elements)
			ax1.vlines(np.arange(len(elements)),0,1)
			ax1.set_ylabel("fractional feedback")
			ax1.set_xlabel("element")
			
			ax2 = ax1.twinx()
			l4 = ax2.bar(np.arange(len(elements)),total_feedback, color = 'k', alpha = 0.2 ,label='total', width = 1)
			ax2.set_yscale('log')
			ax2.set_ylabel('total mass fed back')
			lines = [l1,l2,l3,l4]
			labels = ['sn2', 'sn1a', 'agb', 'total']
			plt.legend(lines, labels,loc='upper right',numpoints=1).get_frame().set_alpha(0.5)
			plt.savefig('processes_%s.png' %(name_string))
		else:
			plt.clf()
			fig = plt.figure(figsize=(11.69,8.27), dpi=100)
			ax1 = fig.add_subplot(111)
			plt.bar(np.arange(len(elements))-0.05,medians[0], yerr=stds[0],error_kw=dict(elinewidth=2,ecolor='k'), color = 'b' ,label='sn2', width = 0.2)
			plt.bar(np.arange(len(elements))+0.2,medians[1], yerr=stds[1],error_kw=dict(elinewidth=2,ecolor='k'), color = 'g' ,label='sn1a', width = 0.2)
			plt.bar(np.arange(len(elements))+0.45,medians[2], yerr=stds[2],error_kw=dict(elinewidth=2,ecolor='k'), color = 'r' ,label='agb', width = 0.2)
			plt.bar(np.arange(len(elements))+0.45,-medians[2], yerr=stds[2],error_kw=dict(elinewidth=2,ecolor='k'), color = 'r' ,alpha = 0.5, width = 0.2)
			plt.bar(np.arange(len(elements))-0.15,medians[3], yerr=stds[3],error_kw=dict(elinewidth=2,ecolor='y'), color = 'y' ,alpha = 0.1,label = 'total', width = 1)
			plt.ylim((1e-5,1e7))
			plt.xticks(np.arange(len(elements))+0.4, elements)
			plt.vlines(np.arange(len(elements))-0.15,0,1e7)
			plt.yscale('log')
			plt.ylabel("total feedback (arbitrary units)")
			plt.xlabel("element")
			plt.legend(loc='upper right',numpoints=1).get_frame().set_alpha(0.5)
			plt.savefig('processes_%s.png' %(name_string))

	return [probability]

def save_abundances(summary_pdf,name_string,abundances):
	'''
	a function that saves the abundances in the current directory

	INPUT:
	
	   summary_pdf = boolean should the abundances be saved?
	
	   name_string = name for the saved file
	
	   abundances = the abundance instance derived from the cube_class

	OUTPUT:
	
	   Saves the abundances as a npy file        
	'''
	if summary_pdf:
		np.save('abundances_%s' %(name_string),abundances)
	return [0.]

def produce_wildcard_stellar_abundances(stellar_identifier, age_of_star, sigma_age, element_symbols, element_abundances, element_errors):
    '''
    This produces a structured array that can be used by the Chempy wildcard likelihood function.
    
    INPUT:
    
       stellar_identifier = name of the files
    
       age_of_star = age of star in Gyr
    
       sigma_age = the gaussian error of the age (so far not implemented in the likelihood function)
    
       element_symbols = a list of the element symbols
    
       element_abundances = the corresponding abundances in [X/Fe] except for Fe where it is [Fe/H]
    
       element_errors = the corresponding gaussian errors of the abundances
    
    OUTPUT:
    
       it will produce a .npy file in the current directory with stellar_identifier as its name.
    '''    
    names = element_symbols + ['age']
    base = np.zeros(2)
    base_list = []
    for item in names:
        base_list.append(base)
    wildcard = np.core.records.fromarrays(base_list,names=names)
    for i,item in enumerate(element_symbols):
        wildcard[item][0] = element_abundances[i]
        wildcard[item][1] = element_errors[i]
    wildcard['age'][0] = age_of_star
    wildcard['age'][1] = sigma_age
    np.save(stellar_identifier,wildcard)

def plot_abundance_wildcard(stellar_identifier, wildcard,abundance_list, element_list, probabilities, time_model):
    '''
    this function plots the abundances of a stellar wildcard and the abundances of a chempy model together with the resulting likelihood values

    INPUT:
    
       stellar_identifier = str, name of the Star
    
       wildcard = the wildcard recarray
    
       abundance_list = the abundance list from the model
    
       element_list = the corresponding element symbols
    
       probabilities = the corresponding single likelihoods as a list
    
       time_model = the time of the chempy model which is compared to the stellar abundance (age of the star form the wildcard)

    OUTPUT:
    
       This saves a png plot to the current directory
    '''
    star_abundance_list = []
    star_error_list = []
    for item in element_list:
        star_abundance_list.append(float(wildcard[item][0]))
        star_error_list.append(float(wildcard[item][1]))
        
    text_size = 14
    plt.rc('font', family='serif',size = text_size)
    plt.rc('xtick', labelsize=text_size)
    plt.rc('ytick', labelsize=text_size)
    plt.rc('axes', labelsize=text_size, lw=1.0)
    plt.rc('lines', linewidth = 1)
    plt.rcParams['ytick.major.pad']='8'
    plt.clf()
    fig = plt.figure(figsize=(10.69,8.27), dpi=100)
    ax = fig.add_subplot(111)
    plt.errorbar(np.arange(len(element_list)),star_abundance_list,xerr=None,yerr=star_error_list,mew=3,marker='x',capthick =3,capsize = 20, ms = 10,elinewidth=3,label = stellar_identifier)
    plt.plot(np.arange(len(element_list)),np.array(abundance_list),label='prediction after %.2f Gyr' %(time_model),linestyle='-')

    for i in range(len(element_list)):
        plt.annotate(xy=(i,-0.4),s= '%.2f' %(probabilities[i]))
    element_labels = ['[%s/Fe]' %(item) for item in element_list]
    element_labels = np.array(element_labels)
    element_labels[np.where(element_labels == '[Fe/Fe]')] = '[Fe/H]'
    plt.xticks(np.arange(len(element_list)), element_labels)
    plt.ylabel("abundance relative to solar in dex")
    plt.xlabel("Element")
    plt.title('joint probability of agreeing with %s (normed to pmax) = %.2f' %(stellar_identifier,sum(probabilities)))	
    plt.legend(loc='best',numpoints=1).get_frame().set_alpha(0.5)
    plt.savefig('%s.png' %(stellar_identifier))
	

def wildcard_likelihood_function(summary_pdf, stellar_identifier, abundances):    
    '''
    This function produces Chempy conform likelihood output for a abundance wildcard that was produced before with 'produce_wildcard_stellar_abundances'.
    
    INPUT:
    
       summary_pdf = bool, should there be an output
    
       stellar_identifier = str, name of the star
    
       abundances = the abundances instance from a chempy chemical evolution

    OUTPUT:
    
       probabilities = a list of the likelihoods for each element
    
       abundance_list = the abundances of the model for each element
    
       element_list = the symbols of the corresponding elements

    These list will be used to produce the likelihood and the blobs. See cem_function.py
    '''

    wildcard = np.load(stellar_identifier + '.npy')
    
    star_time = abundances['time'][-1] - wildcard['age'][0]
    cut = [np.where(np.abs(abundances['time'] - star_time) == np.min(np.abs(abundances['time'] - star_time)))]
    if len(cut[0][0]) != 1:
        cut = cut[0][0][0]
    time_model = abundances['time'][cut]
    
    abundance_list = []
    element_list = []
    probabilities = []
    for item in list(wildcard.dtype.names):
        if item in list(abundances.dtype.names):
            if item != 'Fe':
                abundance = abundances[item][cut]-abundances['Fe'][cut]
            else:
                abundance = abundances['Fe'][cut]
            element_list.append(item)
            abundance_list.append(float(abundance))
            probabilities.append(float(gaussian_1d_log(abundance,wildcard[item][0],wildcard[item][1])))
    if summary_pdf:
        plot_abundance_wildcard(stellar_identifier,wildcard,abundance_list, element_list, probabilities, time_model)
    return probabilities, abundance_list, element_list


def likelihood_function(stellar_identifier, list_of_abundances, elements_to_trace):    
	'''
	This function calculates analytically an optimal model error and the resulting likelihood from the comparison of predictions and observations

	INPUT:

	   list_of_abundances = a list of the abundances coming from Chempy

	   elements_to_trace = a list of the elemental symbols

	OUTPUT:

	   likelihood = the added log likelihood

	   element_list = the elements that were in common between the predictions and the observations, has the same sequence as the following arrays
	   
	   model_error = the analytic optimal model error
	   
	   star_error_list = the observed error
	   
	   abundance_list = the predictions
	   
	   star_abundance_list = the observations
	'''
	try:
		wildcard = np.load(stellar_identifier + '.npy')
	except Exception as ex:
		from . import localpath
		wildcard = np.load(localpath + 'input/stars/' + stellar_identifier + '.npy')
	# Brings the model_abundances, the stellar abundances and the associated error into the same sequence
	abundance_list = []
	element_list = []
	star_abundance_list = []
	star_error_list = []
	for i,item in enumerate(elements_to_trace):
		if item in list(wildcard.dtype.names):
			element_list.append(item)
			abundance_list.append(float(list_of_abundances[i]))
			star_abundance_list.append(float(wildcard[item][0]))
			star_error_list.append(float(wildcard[item][1]))
	abundance_list = np.hstack(abundance_list)
	star_abundance_list = np.hstack(star_abundance_list)
	star_error_list = np.hstack(star_error_list)
	assert len(abundance_list) == len(star_abundance_list) == len(star_error_list), 'no equal length, something went wrong'
	###################

	# Introduces the optimal model error which can be derived analytically
	model_error = []
	for i, item in enumerate(element_list):
		if (abundance_list[i] - star_abundance_list[i]) * (abundance_list[i] - star_abundance_list[i]) <= star_error_list[i] * star_error_list[i]:
			model_error.append(0.)
		else:
			model_error.append(np.sqrt((abundance_list[i] - star_abundance_list[i]) * (abundance_list[i] - star_abundance_list[i]) - star_error_list[i] * star_error_list[i]))
	model_error = np.hstack(model_error)

	## Uncomment next line if you do not want to use model errors
	#model_error = np.zeros_like(model_error)

	# Now the likelihood is evaluated
	likelihood = likelihood_evaluation(model_error, star_error_list, abundance_list, star_abundance_list)
	return likelihood, element_list, model_error, star_error_list, abundance_list, star_abundance_list