.. pyphot documentation master file, created by
   sphinx-quickstart on Wed Oct  5 11:25:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Chempy
===============================================================

This package is a flexible elemental abundance fitting tool for stellar abundances.
It takes nucleosynthetic yield tables from literature, integrates them over the IMF for simple stellar (SSP) populations
and uses successive SSPs to integrate the chemical evolution of an open-box one-zone inter stellar medium (ISM) over time.

Package main content
~~~~~~~~~~~~~~~~~~~~

This package is mostly organized around:

* :class:`Chempy.wrapper.SSP_wrap` that handles the SSP enrichment.
 
* :py:func:`Chempy.wrapper.Chempy` that handles the time-integration. 


Additionally, and for convenience 

* :py:func:`Chempy.wrapper.mcmc` can be used to infer most likely parameters.

To get started with using Chempy it is strongly recommended to browse through the extensive tutorial_.

.. _tutorial: https://github.com/jan-rybizki/Chempy/tree/master/tutorials  

Contents
---------

.. toctree::
   :maxdepth: 3

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

Look into the very extensive jupyter tutorial_. Here is a short preview:

.. _tutorial: https://github.com/jan-rybizki/Chempy/tree/master/tutorials

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

