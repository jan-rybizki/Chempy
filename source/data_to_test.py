import matplotlib.pyplot as plt
import numpy as np


def yield_plot(name_string, yield_class, solar_class, element):
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

	plt.savefig('output/yields_%s.png' %(name_string),bbox_inches='tight')

def yield_comparison_plot(yield_name1, yield_name2, yield_class, yield_class2, solar_class, element):
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

	plt.savefig('output/yields_comparison_%s_vs_%s_for_%s.png' %(yield_name1,yield_name2,element),bbox_inches='tight')

def fractional_yield_comparison_plot(yield_name1, yield_name2, yield_class, yield_class2, solar_class, element):
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

	plt.savefig('output/fractional_yields_comparison_%s_vs_%s_for_%s.png' %(yield_name1,yield_name2,element),bbox_inches='tight')


def elements_plot(name_string,agb, sn2, sn1a,elements_to_trace, all_elements,max_entry):
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
	for i, item in enumerate(all_elements["Symbol"]):
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

	plt.savefig('output/elements_%s.png' %(name_string),bbox_inches='tight')
	return [0.]