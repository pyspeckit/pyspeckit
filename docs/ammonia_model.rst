Ammonia Models
==============

The Ammonia modeling tools include a set of constants in the
``ammonia_constants`` module and the following ammonia modeling tools
listed below.

There is also an ammonia fitter wrapper; see :doc:`wrappers`.

Note that there are two modules described here: the multi-rotational-transition
fitter, which has its own set of custom functions, and a generic hyperfine-line
fitting module meant to fit a single metastable (or non-metastable) transition.

If you run into a problem, have a look at the `ammonia-tagged issues on github
<https://github.com/pyspeckit/pyspeckit/issues?q=label%3Aammonia+>`__

Note that the default ammonia spectrum, with only the CMB as a background, is
intentionally restricted to be positive to avoid nonphysical (absorption
against the CMB) parameters.  If you find your spectral fits sometimes
returning errors about a negative spectrum, you may need the 'restricted
excitation temperature' model.  In this model, the excitation temperature
is no longer a free parameter, but instead is set to trot + deltat, where
deltat is a positive value.  The only difference from the main ammonia model
is this restriction, which can result in more numerically stable results.
An example can be found `in the examples directory
<https://github.com/pyspeckit/pyspeckit/blob/master/examples/synthetic_nLTE_ammonia_spectrum_example_witherrorestimates.py>`__

.. automodule:: pyspeckit.spectrum.models.ammonia
    :members:

.. automodule:: pyspeckit.spectrum.models.ammonia_hf
    :members:
