Fitting a continuum model as a model
====================================

.. include:: <isogrk3.txt>

This example shows the initialization of a pyspeckit object from numpy arrays,
as in :doc:`example_fromscratch`, but it adds the twist of including a steep
continuum.

We fit the continuum using the polynomial continuum model, which gives access
to the error on the polynomial fit parameters.  No such parameters are
accessible via the `pyspeckit.Spectrum.baseline` tools because they use
`numpy.poly1d` to fit the data.

.. include:: example_continuum_fromscratch.py
   :literal:

