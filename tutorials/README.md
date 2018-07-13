## Chempy tutorial roadmap

You can load the notebooks using ``jupyter notebook`` and then change the content interactively, once you downloaded this repository.
You can also just have a look at them here on github.

In **tutorial 1,2** and **3** you will learn the basic concepts of galactic chemical evolution (GCE), like the solar abundances, stellar initial mass function (IMF), star formation rate (SFR) or nucleosynthetic yields.

**Tutorial 4** shows how different the enrichment from a simple stellar population (SSP) can be when using different yield sets and also different IMFs.

**Tutorial 5** illustrates the usage of *Chempy* to produce a complete one-zone chemical evolution and how the results depend on the IMF, SFR and yield set.
You will also see how to compare one-zone models with a distribution of observed stars. We will also plot the fractional contribution of each nucleosynthetic channel to the elements.

In **tutorial 6** the parameter inference using the MCMC is shown. You will learn how to run the MCMC for a specific set of parameters and how to plot the posterior distribution.
You will also see how the prediction correlate and that elements contain redundant information.
Finally the wildcard is introduced with which you can optimize *Chempy* chemical evolution parameters for your own stellar abundances.

**Tutorial 7** load the abundance tracks from the models of Paper 1 in order for you to be able to compare to it.

In **tutorial 8** it is outlined how to produce IMF-integrated and metallicity dependent SSP yield tables that you can plug into your own simulation to build on the flexible Chempy framework. The Yield table can use arbitrary time-steps, it can track the different enrichment processes. You can also change between net and gross yields.
