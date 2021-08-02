# Chempy
Flexible one-zone open-box chemical evolution modeling. Abundance fitting and stellar feedback calculation

## Recent Developments


[Oliver Philcox](https://github.com/oliverphilcox), in 2019 produced another Chempy [paper](https://ui.adsabs.harvard.edu/abs/2019arXiv190900812P/abstract) focusing on multi-star inference.

[Oliver Philcox](https://github.com/oliverphilcox), during his 2017 summer internship at MPIA coded a [NeuralNet add-on of Chempy](https://github.com/oliverphilcox/ChempyScoring) (together with a very nice jupyter tutorial), which is orders of magnitudes faster than the original version and used it to score different yield-tables from the literature which lead to this [publication](http://adsabs.harvard.edu/abs/2018ApJ...861...40P). He also included many more CC-SN yieldsets.

Kirsten Blancato in this [paper](https://ui.adsabs.harvard.edu/abs/2019ApJ...883...34B/abstract) used single zone models per star and got varying IMF high mass slopes across the chemical abundance space.

[Nathan Sandford](https://github.com/NathanSandford), produced an interactive and very instructive [Widget](https://hub.mybinder.org/user/nathansandford-chempy-widget-ibik9tdn/notebooks/chempy_widget.ipynb) which you can run in your browser to see the effect that the star formation history has on abundance patterns.

## Installation

```
pip install git+https://github.com/jan-rybizki/Chempy.git
```
Chempy should run with the latest python 2 and python 3 version.
Its dependencies are: [Numpy](http://numpy.scipy.org/), [SciPy](http://www.scipy.org/), [matplotlib](http://matplotlib.sourceforge.net/), [multiprocessing](https://docs.python.org/2/library/multiprocessing.html#module-multiprocessing) and [emcee](http://dan.iel.fm/emcee/current/) (for the MCMC), and [corner](http://corner.readthedocs.io/en/latest/) (for the MCMC plots). They are all pip installable and you can also get part of it with [Anaconda](https://www.continuum.io/downloads).

### Installation without admin rights:
You can install *Chempy* into a folder where you have write access:
```
pip install --install-option='--prefix=~/extra_package/' git+https://github.com/jan-rybizki/Chempy.git
```
Then you have to add the `site-packages/` folder which will be one of the newly created subfolders in `extra_package/` into the ```PYTHONPATH``` variable, e.g.:
```
export PYTHONPATH=~/extra_package/lib/python2.7/site-packages/:$PYTHONPATH
```
If you want this to be permanent, you can add the last line to your `.bashrc`.


## Authors
- Jan Rybizki (MPIA, rybizki@mpia.de)

## Collaborators
- Hans-Walter Rix (MPIA)
- Andreas Just (ZAH)
- Morgan Fouesneau (MPIA)

## Links
- <a href="http://arxiv.org/abs/1702.08729"><img src="http://img.shields.io/badge/arXiv-1702.08729-orange.svg?style=flat" alt="arxiv:1702.08729" /></a>
- <a href="http://ascl.net/1702.011"><img src="https://img.shields.io/badge/ascl-1702.011-blue.svg?colorB=262255" alt="ascl:1702.011" /></a>
- An early version of Chempy is presented in chapter 4 of my [phd thesis](http://nbn-resolving.de/urn:nbn:de:bsz:16-heidok-199349).

## Getting started
The jupyter [tutorial](https://github.com/jan-rybizki/Chempy/tree/master/tutorials) illustrates the basic usage of Chempy and basic concepts of galactic chemical evolution modeling. It can be inspected in the github repository or you can run it interactively on your local machine.

To run it interactively first clone the repository with
```
git clone https://github.com/jan-rybizki/Chempy.git
```
Then you can ```jupyter notebook``` from within the tutorial folder (it will run if you have installed *Chempy*). 
If you did not install Chempy you can still run the tutorial but need to point to the files in the Chempy folder. Basically you have to ```cd ../Chempy/``` and then replace each ```from Chempy import ...``` with ```from . import ...```.

You can also read the automatically generated [manual](https://chempy.readthedocs.io/en/latest/).

## Compare to Chempy data
If you want to compare your abundance model/data to Chempy paper one results, look at the [tutorial 7](https://github.com/jan-rybizki/Chempy/blob/master/tutorials/7-Acessing%20Chempy%20paper%201%20abundance%20tracks.ipynb) where the stored abundance tracks are loaded and plotted for one element.

## Extract yield tables for chemical evolution
If you want to use the flexible framework of Chempy to produce IMF integrated metallicity dependent yield tables for your SPH or other Chemical Evolution model you can use [tutorial 8](https://github.com/jan-rybizki/Chempy/blob/master/tutorials/8-Yield%20tables%20for%20SPH%20simulations%20and%20comparison%20to%20other%20tables.ipynb). You can use net or gross yields and also look at individual processes contribution to the overall SSP yield table.

## Attribution
Please cite the [paper](http://adsabs.harvard.edu/abs/2017A%26A...605A..59R) when using the code in your research.


\bibitem[Rybizki et al.(2017)]{2017A&A...605A..59R} Rybizki, J., Just, A., \& Rix, H.-W.\ 2017, \aap, 605, A59
