Baseline Fitting
----------------

There are a number of cool features in baselining that aren't well-described
below, partly due to Sphinx errors as of 12/22/2011.

``exclude`` and ``include`` allow you to specify which parts of the spectrum to use
for baseline fitting.  Enter values as pairs of coordinates.  

Excludefit makes use of an existing fit and excludes all points with signal above
a (very low) threshold when fitting the baseline.  Going back and forth between
``baseline(excludefit=True)`` and ``specfit()`` is a nice way to iteratively measure
the baseline & emission/absorption line components.

API
~~~

.. automodule:: pyspeckit.spectrum.baseline
.. autoclass:: Baseline
    :members:
    :undoc-members:
    :special-members:
