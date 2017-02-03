import sys,time,os,inspect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import shutil
import subprocess
from os.path import isfile, join


with_blobs = True
producing_positions_for_plot = True
how_many_plot_samples = 100
how_many_MCMC_samples = 5000
#how_many_MCMC_samples = 1000
plot_mcmc_space = True



directory = 'mcmc/0/'

set_scale = False
use_scale = True

def restructure_chain(directory, set_scale, use_scale):

	positions = np.load('%sflatchain.npy' %(directory))
	posterior = np.load('%sflatlnprobability.npy' %(directory))
	mean_posterior = np.load('%sflatmeanposterior.npy' %(directory))
	std_posterior = np.load('%sflatstdposterior.npy' %(directory))
	if with_blobs:
		blobs = np.load('%sflatblobs.npy' %(directory))
		blobs = np.swapaxes(blobs,0,1)
	if len(blobs.shape)!=3:
		print 'blob shape = ', blobs.shape, 'probably some runs did not return results and were stored anyway.'
		with_blobs = False

	nwalkers = positions.shape[0]
	dimensions = positions.shape[2]
	iterations = positions.shape[1]
	print 'we have %d iterations' %(iterations)
	if with_blobs:
		assert nwalkers == blobs.shape[0]
		assert iterations == blobs.shape[1]
		blob_dimensions = blobs.shape[2]
	if posterior.shape[0] != positions.shape[0] or posterior.shape[1] != positions.shape[1]:
		raise Exception('chain and probability do not have the same number of walkers and or iterations')


	plt.figure(figsize=(10.69,8.27), dpi=100)
	plt.plot(mean_posterior+10.-np.max(mean_posterior), label= 'mean (maximum shifted to 10)')
	plt.plot(std_posterior, label= 'std')
	plt.yscale('log')
	plt.ylim((1.,None))
	plt.title('Statistical moments of the %d walkers at each step' %(nwalkers))
	plt.legend(loc = 'best')
	plt.grid('on')
	plt.savefig('%sposterior_evolution.png' %(directory))
	plt.clf()
	print np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[-1] 
	if True:#np.mean(posterior, axis = 0)[0] < np.mean(posterior, axis = 0)[-1]:
		print 'chain is inverted' ##Sometimes the chain is stored differently depending on the system
		posterior = posterior[:,::-1]
		positions = positions[:,::-1]
		if with_blobs:
			blobs = blobs[:,::-1]

	print np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[-1] 
	
	keeping = int(how_many_MCMC_samples / nwalkers)
	positions = positions[:,:keeping, :]
	posterior = posterior[:,:keeping]
	if with_blobs:
		blobs = blobs[:,:keeping,:]
	print 'now the burn-in tail is cut out:'
	print np.mean(posterior, axis = 0)[0], np.mean(posterior, axis = 0)[-1] 
	### shaping back
	positions = positions.reshape((-1, dimensions), order = 'F')
	posterior = posterior.reshape(-1, order = 'F')
	if with_blobs:
		blobs = blobs.reshape((-1, blob_dimensions), order = 'F')


	print 'after throwing out burn in %d iterations' %(len(posterior))
	assert np.any(np.isinf(posterior)) == False

	total_iterations = len(posterior)
	#
	cut = np.where(posterior==np.max(posterior))
	positions_max = positions[cut]
	posterior_max = posterior[cut]
	if with_blobs:
		blobs_max = blobs[cut]
	posterior_max,indices_max = np.unique(posterior_max, return_index=True)
	positions_max = positions_max[indices_max]


	data_cut = np.max(posterior) - 15 ## throwing out very odd values
	cut = np.where(posterior>data_cut)
	positions_before = len(positions)
	positions = positions[cut]
	positions_after = len(positions)
	throw_out = positions_before - positions_after
	posterior = posterior[cut]
	if with_blobs:
		blobs = blobs[cut]
	print 'we have %d iterations good enough posterior' %(len(posterior))
	if throw_out > 0:
		print '%d runs of the stabilised MCMC had a posterior that was worse -15 ln' %(throw_out)
	### Drawing 100 random posterior positions
	if producing_positions_for_plot:
		random_indices = np.random.choice(a = np.arange(len(positions)),size = how_many_plot_samples,replace = False)
		plot_positions = positions[random_indices]
		plot_posterior = posterior[random_indices]
		if with_blobs:
			plot_blobs = blobs[random_indices]
		np.save('positions_for_plotting',plot_positions)
		np.save('%spositions_for_plotting' %(directory),plot_positions)


	np.save('%sposteriorPDF' %(directory),positions)
	np.save('%sposteriorvalues' %(directory),posterior)
	np.save('%sblobs_distribution' %(directory),blobs)
	# IF ARRAYS SHOULD BE SORTED, looks better on scatter plot:
	if True:
		sorted_indices = np.argsort(posterior)
		posterior = posterior[sorted_indices]
		positions = positions[sorted_indices]
		if with_blobs:
			blobs = blobs[sorted_indices]

	vmax = np.max(posterior)
	vmin = np.min(posterior)
	np.save('best_parameter_values',positions_max)
	np.save('%sbest_parameter_values' %(directory),positions_max)
	print vmax,vmin
	print 'vmax was obtained at: ', positions_max
	print 'for plotting we use the best: ',len(posterior), ' values'
	print 'unique posterior values: ', len(np.unique(posterior))

	x = np.array([np.mean(positions,axis = 0),np.std(positions,axis = 0)])
	if with_blobs:
		y = np.array([np.mean(blobs,axis = 0),np.std(blobs,axis = 0)])
	np.savetxt('%sparameter_moments.csv' %(directory), list(np.hstack(x.T).T), fmt='%.4f', delimiter=',')
	if with_blobs:
		np.savetxt('%slikelihood_moments.csv' %(directory), list(np.hstack(y.T).T), fmt='%.4f', delimiter=',')
		np.save('%sblobs_max' %(directory), blobs_max)
		plt.figure(figsize=(30.69,8.27), dpi=100)
		#plt.plot(blobs_max, label= 'blobs of best run')
		plt.plot(y[0], label= 'mean blobs')
		plt.plot(y[1], label= 'std blobs')
		#plt.yscale('log')
		#plt.ylim((1.,None))
		plt.title('Statistical moments of the blobs for the %d plot runs' %(how_many_plot_samples))
		plt.legend(loc = 'best')
		plt.grid('on')
		plt.savefig('%sblobs_distribution.png' %(directory))
		plt.clf()


	parameter_names = [r'$\alpha_\mathrm{IMF}$',r'$\log_{10}\left(\mathrm{N}_\mathrm{Ia}\right)$',r'$\log_{10}\left(\tau_\mathrm{Ia}\right)$',r'$\log_{10}\left(\mathrm{SFE}\right)$',r'$\mathrm{SFR}_\mathrm{peak}$',r'$\mathrm{x}_\mathrm{out}$',r'$\log_{10}\left(\mathrm{f}_\mathrm{corona}\right)$']
	np.save("%sparameter_names" %(directory), parameter_names)
	
	if len(parameter_names) != dimensions:
		raise Exception('parameter_names not equally numbered as parameter in chain')
	for j in range(len(parameter_names)):
		print parameter_names[j], positions[:,j].mean(), '+-', positions[:,j].std()


	fig, axes = plt.subplots(nrows=nparameter+1, ncols=1,figsize=(14.69,30.27), dpi=100,sharex=True)
	for i in range(nparameter):
		axes[i].plot(positions[:,i])
		axes[i].set_ylabel(parameter_names[i])
	axes[i+1].plot(posterior)
	axes[i+1].set_ylabel('posterior')	

	fig.savefig("%schain.png" %(directory))
	plt.clf()
	plt.figure(figsize=(14.69,30.27), dpi=100)
	plt.hist(posterior,bins=100)
	plt.savefig('%sposterior_histo.png' %(directory))
	plt.clf()
	plt.close()

def plot_mcmc_chain(directory, set_scale):

	plt.clf()
	m=plt.cm.jet
	text_size = 16
	cor_text = 22
	plt.rc('font', family='serif',size = text_size)
	plt.rc('xtick', labelsize=text_size)
	plt.rc('ytick', labelsize=text_size)
	plt.rc('axes', labelsize=text_size, lw=1.0)
	plt.rc('lines', linewidth = 1)
	plt.rcParams['ytick.major.pad']='8'
	plt.rcParams['text.latex.preamble']=[r"\usepackage{libertine}"]
	params = {'text.usetex' : True,
	          'font.size' : 10,
	          'font.family' : 'libertine',
	          'text.latex.unicode': True,
	          }
	plt.rcParams.update(params)
	positions = np.load('%sposteriorPDF.npy' %(directory))
	positions_max = np.load('%sbest_parameter_values.npy' %(directory))
	parameter_names = np.load("%sparameter_names" %(directory))

	nparameter = len(positions[0])
	cor_matrix = np.zeros(shape = (nparameter,nparameter))
	for i in range(nparameter):
		for j in range(nparameter):
			cor_matrix[i,j] = np.corrcoef((positions[:,i],positions[:,j]))[1,0]
	np.save('%scor_matrix' %(directory),cor_matrix)

	fig, axes = plt.subplots(nrows=nparameter, ncols=nparameter,figsize=(14.69,8.0), dpi=300)#,sharex=True, sharey=True)

	left  = 0.1  # the left side of the subplots of the figure
	right = 0.925    # the right side of the subplots of the figure
	bottom = 0.075   # the bottom of the subplots of the figure
	top = 0.97      # the top of the subplots of the figure
	wspace = 0.0   # the amount of width reserved for blank space between subplots
	hspace = 0.0   # the amount of height reserved for white space between subplots
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

	alpha=0.5
	alpha_more = 0.8
	alpha_less = 0.2
	lw = 2
	if set_scale:
		borders = []
	if use_scale:
		borders = np.load('prior_borders.npy')
	t = 0
	if plot_mcmc_space:
		for i in range(nparameter):
			for j in range(nparameter):
				axes[i,j].locator_params(nbins=4)
				if j==1:
					axes[i,j].locator_params(nbins=4)
				if i == j:
					counts, edges = np.histogram(positions[:,j], bins=10)
					max_count = float(np.max(counts))
					counts = np.divide(counts,max_count)
					axes[i,j].bar(left = edges[:-1], height = counts, width = edges[1]-edges[0], color = 'grey', alpha = alpha, linewidth = 0, edgecolor = 'blue')
					if use_scale:
						axes[i,j].set_xlim(borders[t][0])
						axes[i,j].plot( borders[t][2], borders[t][3], c="k", linestyle = '--', alpha=1, lw=lw )
						t += 1
					else:
						axes[i,j].set_xlim(min(positions[:,j]),max(positions[:,j]))
					axes[i,j].set_ylim(0,1.05)
					if j != 0:
						plt.setp(axes[i,j].get_yticklabels(), visible=False)

					if set_scale:
						borders.append([axes[i,j].get_xlim(),axes[i,j].get_ylim(),xgauss,ygauss])

					axes[i,j].vlines(np.percentile(positions[:,j],15.865),axes[i,j].get_ylim()[0],axes[i,j].get_ylim()[1], color = 'k',alpha=alpha,linewidth = lw,linestyle = 'dashed')    
					axes[i,j].vlines(np.percentile(positions[:,j],100-15.865),axes[i,j].get_ylim()[0],axes[i,j].get_ylim()[1], color = 'k',alpha=alpha,linewidth = lw,linestyle = 'dashed')  
					axes[i,j].vlines(np.percentile(positions[:,j],50),axes[i,j].get_ylim()[0],axes[i,j].get_ylim()[1], color = 'k',alpha=alpha, linewidth = lw)
					axes[i,j].text( 0.5, 1.03, r'$%.2f_{-%.2f}^{+%.2f}$'%(np.percentile(positions[:,j],50),np.percentile(positions[:,j],50)-np.percentile(positions[:,j],15.865),np.percentile(positions[:,j],100-15.865)-np.percentile(positions[:,j],50)),fontsize=text_size, ha="center" ,transform=axes[i,j].transAxes)

				if i>j:
					if j != 0:
						plt.setp(axes[i,j].get_yticklabels(), visible=False)
					
					#corner.hist2d(positions[:,j],positions[:,i], ax = axes[i,j],bins = 15 , levels=(1-np.exp(-0.5),1-np.exp(-2.0),1-np.exp(-4.5)))
					axes[i,j].plot(positions_max[:,j],positions_max[:,i],'kx',markersize = 10,mew=2.5)
					im = axes[i,j].scatter(positions[:,j],positions[:,i],c=posterior,cmap=m,vmin=vmin,vmax=vmax,edgecolor='None',s=45,marker='o',alpha=alpha_less,rasterized=True)#,vmin=0,vmax=vmax)

					if use_scale:
						axes[i,j].set_xlim(borders[t][0])
						axes[i,j].set_ylim(borders[t][1])
						axes[i,j].plot( borders[t][2], borders[t][3], "k",linestyle = '--', alpha=1, lw = lw )
						t += 1
					else:
						axes[i,j].set_xlim(min(positions[:,j]),max(positions[:,j]))
						axes[i,j].set_ylim(min(positions[:,i]),max(positions[:,i]))


					if set_scale:
						borders.append([axes[i,j].get_xlim(),axes[i,j].get_ylim(),xk,yk])

				if j>i:
					correlation_coefficient = np.corrcoef((positions[:,i],positions[:,j]))
					axes[i,j].text( 0.6, 0.5, "%.2f"%(correlation_coefficient[1,0]),fontsize=cor_text, ha="center" ,transform=axes[i,j].transAxes)	
					axes[i,j].axis('off')
				if i == nparameter-1:
					axes[i,j].set_xlabel(parameter_names[j])
				if j == 0:
					axes[i,j].set_ylabel(parameter_names[i])
		#if nparameter>3:
		#	axes[0,3].set_title('vmax = %.2f, obtained at %s' %(vmax,str(positions_max)))
		#	#axes[0,1].set_title('vmax = %.2f, obtained at %s and %d evals thrown out' %(vmax,str(positions_max),throw_out))
		
		fig.savefig('%sparameter_space_sorted.png' %(directory),dpi=300,bbox_inches='tight')
		plt.clf()

	if set_scale:
		np.save('%sprior_borders' %(directory), borders)
		np.save('prior_borders', borders)

