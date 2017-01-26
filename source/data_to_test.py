import matplotlib.pyplot as plt
import numpy as np

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