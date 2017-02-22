.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Chempy - Flexible one-zone open-box chemical evolution modeling
===============================================================

This package is a elemental abundance fitting tool for stellar abundances.
It takes literature yield, integrates them over the IMF for simple stellar (SSP) populations
and uses successive SSPs to integrate the chemical evolution of a single zone ISM over time.

Package main content
~~~~~~~~~~~~~~~~~~~~

This package is mostly organized around 2 main classes:

* :class:`Chempy.weighted_yield.SSP` that handles the SSP enrichment.
 
* :class:`Chempy.cem_function.cem` that handles the time-integration. 


Additionally, and for convenience 

* the mcmc function can use cem to infer most likely parameters
  :class:`Chempy.wrapper.mcmc`

Contents
---------

.. toctree::
   :maxdepth: 1

   modules

Installation
------------

* Using pip:

.. code-block:: none
	
	pip install git+https://github.com/jan-rybizki/Chempy.git

* Manually:

.. code-block:: none

        git clone https://github.com/jan-rybizki/Chempy
        cd Chempy
        python setup.py intall



Quick Start
~~~~~~~~~~~

.. code-block:: python

	# Loading the default parameters
	from Chempy.parameter import ModelParameters
	a = ModelParameters()        

	# Evaluating the Chempy posterior at the maximum prior parameters
	from Chempy.cem_function import cem
	a.testing_output = True
	a.summary_pdf = True
	a.observational_constraints_index = ['gas_reservoir','sn_ratio','sol_norm']
	posterior, blobs = cem(a.p0,a)	

-----------------

References
----------

* Chempy paper


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

