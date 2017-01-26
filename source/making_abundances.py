import numpy as np



def mass_fraction_to_abundances(cube, solar_abundances):
	element_names = list(set(solar_abundances['Symbol']).intersection(cube.dtype.names))
	element_number = []
	element_masses = []
	for item in element_names:
		element_number.append(int(solar_abundances['Number'][np.where(solar_abundances['Symbol']==item)]))
		element_masses.append(solar_abundances['Mass'][np.where(solar_abundances['Symbol']==item)])
	
	sorted_index = np.argsort(np.array(element_number))
	element_number = [element_number[i] for i in sorted_index]
	element_masses = [element_masses[i] for i in sorted_index]
	element_names = [element_names[i] for i in sorted_index]

	base = np.zeros(len(cube))
	list_of_arrays = []
	for i in range(len(element_names)):
		list_of_arrays.append(base)

	cube_abundances = np.core.records.fromarrays(list_of_arrays,names=element_names)

	for i,item in enumerate(element_names):
		cube_abundances[item] = np.divide(cube[item],float(element_masses[i]))

	normalisation = np.copy(cube_abundances['H'])

	for i,item in enumerate(element_names):
		cube_abundances[item] = np.divide(cube_abundances[item],normalisation)
	
	for i,item in enumerate(element_names):
		cube_abundances[item] = np.log10(cube_abundances[item]) + 12.

	for i,item in enumerate(element_names):
		cube_abundances[item] -= solar_abundances['photospheric'][np.where(solar_abundances['Symbol']==item)] 

	return (cube_abundances,element_names,element_number)

def abundance_to_mass_fraction(all_elements,all_masses,all_abundances,abundances,symbols):
	fractions = []
	for i,item in enumerate(symbols):
		fractions.append(abundances[i])
		fractions[i] -= 12
		fractions[i] = np.power(10,fractions[i])
		fractions[i] *= all_masses[np.where(all_elements == item)]
	tmp = sum(fractions)
	for i,item in enumerate(symbols):
		fractions[i] /= tmp
	return np.hstack(fractions)

def abundance_to_mass_fraction_normed_to_solar(all_elements,all_masses,all_abundances,abundances,symbols):
	fractions = []
	for i,item in enumerate(symbols):
		fractions.append(abundances[i] + all_abundances[np.where(all_elements == item)])
		fractions[i] -= 12
		fractions[i] = np.power(10,fractions[i])
		fractions[i] *= all_masses[np.where(all_elements == item)]
	tmp = sum(fractions)
	for i,item in enumerate(symbols):
		fractions[i] /= tmp
	return np.hstack(fractions)