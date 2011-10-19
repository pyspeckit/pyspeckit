Models
==================

.. module:: models
.. currentmodule:: pyspeckit.spectrum.models
.. automodule:: pyspeckit.spectrum.models
    :members:

The generic SpectralModel class is a wrapper for model functions.  A model
should take in an X-axis and some number of parameters.  In order to declare a
SpectralModel, you give SpectralModel the function name and the number of
parameters it requires.  The rest of the options are optional, though parnames
& shortvarnames are strongly recommended.  If you do not specify fitunits,
your fitting code must deal with units internally.::

    hill5_fitter = model.SpectralModel(hill5_model, 5,
            parnames=['tau', 'v_lsr',  'v_infall',  'sigma', 'tpeak'], 
            parlimited=[(True,False),(False,False),(True,False),(True,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=("\\tau","v_{lsr}","v_{infall}","\\sigma","T_{peak}"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )


.. autoclass:: pyspeckit.spectrum.models.model.SpectralModel
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:

.. autoclass:: pyspeckit.spectrum.models.hyperfine.hyperfinemodel
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:

.. autoclass:: pyspeckit.spectrum.models.fitter.SimpleFitter
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:

.. automodule:: pyspeckit.spectrum.models
    :members:

.. automodule:: pyspeckit.spectrum.models.ammonia
    :members:
.. automodule:: pyspeckit.spectrum.models.fitter
    :members:
.. automodule:: pyspeckit.spectrum.models.formaldehyde
    :members:
.. automodule:: pyspeckit.spectrum.models.gaussfitter
    :members:
.. automodule:: pyspeckit.spectrum.models.hcn
    :members:
.. automodule:: pyspeckit.spectrum.models.hill5infall
    :members:
.. automodule:: pyspeckit.spectrum.models.hyperfine
    :members:
.. automodule:: pyspeckit.spectrum.models.lorentzian
    :members:
.. automodule:: pyspeckit.spectrum.models.modelgrid
    :members:
.. automodule:: pyspeckit.spectrum.models.n2hp
    :members:
.. automodule:: pyspeckit.spectrum.models.voigtfitter
    :members:

