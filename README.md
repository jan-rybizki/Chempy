# Chempy
Flexible one-zone open-box chemical evolution modeling. Abundance fitting and stellar feedback calculation

## Installation

```
pip install git+https://github.com/jan-rybizki/Chempy.git@v0.1
```
Chempy should run with the latest python 2 and python 3 version.
Its dependencies are: [Numpy](http://numpy.scipy.org/), [SciPy](http://www.scipy.org/), [matplotlib](http://matplotlib.sourceforge.net/), [multiprocessing](https://docs.python.org/2/library/multiprocessing.html#module-multiprocessing) and [emcee](http://dan.iel.fm/emcee/current/) (for the MCMC), and [corner](http://corner.readthedocs.io/en/latest/) (for the MCMC plots). They are all pip installable and you can also get part of it with [Anaconda](https://www.continuum.io/downloads).

### Installation without admin rights:
You can install *Chempy* into a folder where you have write access:
```
pip install --install-option='--prefix=~/extra_package/' git+https://github.com/jan-rybizki/Chempy.git@v0.1
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
- [Paper on overleaf](https://www.overleaf.com/read/jypbdqqqrdnp), Arxiv link following
- An early version of Chempy is presented in chapter 4 of my [phd thesis](http://nbn-resolving.de/urn:nbn:de:bsz:16-heidok-199349).

## Getting started
The jupyter [tutorial](https://github.com/jan-rybizki/Chempy/tree/master/tutorials) illustrates the basic usage of Chempy and basic concepts of galactic chemical evolution modeling. It can be inspected in the github repository or you can run it interactively on your local machine.

To run it interactively first clone the repository with
```
git clone https://github.com/jan-rybizki/Chempy.git@v0.1
```
Then you can ```jupyter notebook``` from within the tutorial folder (it will run if you have installed *Chempy*). 
If you did not install Chempy you can still run the tutorial but need to point to the files in the Chempy folder. Basically you have to ```cd ../Chempy/``` and then replace each ```from Chempy import ...``` with ```from . import ...```.

You can also have a look at the *preliminary* [documentaion](http://www.mpia.de/homes/rybizki/html/index.html) which gives an overview over the Chempy classes and functions.
