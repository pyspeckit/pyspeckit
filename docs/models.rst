Models
==================
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

.. TODO::
    Implement the gaussian-hermite profile described here:
    `<http://pipelinesandarchives.blogspot.com/2012/09/fit1d-new-smurf-command-for-acsis-data.html>`_

API Documentation for Models
----------------------------

.. automodule:: pyspeckit.spectrum.models.model
    :members:
.. automodule:: pyspeckit.spectrum.models.ammonia
    :members:
.. automodule:: pyspeckit.spectrum.models.fitter
    :members:
.. automodule:: pyspeckit.spectrum.models.formaldehyde
    :members:
.. automodule:: pyspeckit.spectrum.models.inherited_gaussfitter
    :members:
.. automodule:: pyspeckit.spectrum.models.hcn
    :members:
.. automodule:: pyspeckit.spectrum.models.hill5infall
    :members:
.. automodule:: pyspeckit.spectrum.models.hyperfine
    :members:
.. automodule:: pyspeckit.spectrum.models.inherited_lorentzian
    :members:
.. automodule:: pyspeckit.spectrum.models.modelgrid
    :members:
.. automodule:: pyspeckit.spectrum.models.n2hp
    :members:
.. automodule:: pyspeckit.spectrum.models.inherited_voigtfitter
    :members:
.. automodule:: pyspeckit.spectrum.models.hydrogen
    :members:

Model Documentation Table of Contents
-------------------------------------
.. toctree::
   :maxdepth: 3

   Parameters <parameters>

