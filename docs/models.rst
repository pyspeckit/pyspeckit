Models
======
See :doc:`parameters` for information on how to restrict/modify model parameters.

.. module:: pyspeckit.spectrum.models

The generic SpectralModel class is a wrapper for model functions.  A model
should take in an X-axis and some number of parameters.  In order to declare a
SpectralModel, you give SpectralModel the function name and the number of
parameters it requires.  The rest of the options are optional, though parnames
& shortvarnames are strongly recommended.  If you do not specify fitunits,
your fitting code must deal with units internally.

Here are some examples of how to make your own fitters::

    hill5_fitter = model.SpectralModel(hill5_model, 5,
            parnames=['tau', 'v_lsr',  'v_infall',  'sigma', 'tpeak'], 
            parlimited=[(True,False),(False,False),(True,False),(True,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)],
            # specify the parameter names (TeX is OK)
            shortvarnames=("\\tau","v_{lsr}","v_{infall}","\\sigma","T_{peak}"), 
            fitunits='Hz' )

    gaussfitter = model.SpectralModel(gaussian, 3,
            parnames=['amplitude','shift','width'], 
            parlimited=[(False,False),(False,False),(True,False)], 
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=('A',r'\Delta x',r'\sigma'))

Then you can `register <registration.html>`_ these fitters.

Fitting
-------

Once you have a model defined, you can fit it using the
`pyspeckit.Spectrum.specfit` module.  Documents on fitting have not been
prepared yet, but you can learn most of the tricks by looking at the various
fitting examples and the :doc:`parameters` documentation.

See also :doc:`fitting`.


.. TODO::
    Implement the gaussian-hermite profile described here:
    `<http://pipelinesandarchives.blogspot.com/2012/09/fit1d-new-smurf-command-for-acsis-data.html>`_


Related Documents
-----------------

.. toctree::
   :maxdepth: 1

   Parameters <parameters>

Generic Models and Tools
------------------------

.. toctree::
   :maxdepth: 1

   Gaussian model <gaussian_model>
   lorentzian_model
   voigt_model
   LTE Molecule Model <lte_molecule_model>
   Hyperfine Line model <hyperfine_model>
   modelgrid
   template_model


Specific Models
---------------

.. toctree::
   :maxdepth: 1

   Ammonia Temperature and Hyperfine model <ammonia_model>
   Formaldehyde model <formaldehyde_model>
   HCN model <hcn_model>
   hill5infall_model
   n2hp_model
   hydrogen_model


API Documentation for Models
----------------------------

We include the API documentation for the generic model and fitter wrappers here.

.. automodule:: pyspeckit.spectrum.models.model
    :members:

.. automodule:: pyspeckit.spectrum.models.fitter
    :members:
