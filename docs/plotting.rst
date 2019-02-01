Basic Plotting Guide
====================

The plotting tool in pyspeckit is intended to make publication-quality plots
straightforward to produce.

For details on the various plotting tools, please see the `examples` and the
`plotter documentation <pyspeckit.spectrum.plotters.Plotter>`.

A few basic examples are shown in the snippet below, with comments describing
the various steps.  This example (https://github.com/jmangum/spectrumplot)
shows how to use pyspeckit with `spectral_cube
<https://spectral-cube.readthedocs.io/en/latest/index.html>`_
to extract and label spectra.  Other examples:


.. literalinclude:: basic_plot.py
   :language: python



Basic plot example:


.. figure:: images/basic_plot_example.png
        :alt: Basic plot example
        :figwidth: 800
        :width: 800


Basic plot example with a fit and an annotation (annotation is on by default):


.. figure:: images/basic_plot_example_withfit.png
        :alt: Basic plot example with a fit and an annotation (default)
        :figwidth: 800
        :width: 800


Basic plot example with a fit, but with no annotation:


.. figure:: images/basic_plot_example_withfit_no_annotation.png
        :alt: Basic plot example with a fit, but with no annotation
        :figwidth: 800
        :width: 800


Basic plot example with a second spectrum overlaid in green: 


.. figure:: images/basic_plot_example_with_second_spectrum_overlaid_in_green.png
        :alt: Basic plot example with a second spectrum overlaid in green
        :figwidth: 800
        :width: 800


Basic plot example with a second spectrum overlaid in green plus adjusted limits:


.. figure:: images/basic_plot_example_with_second_spectrum_overlaid_in_green_wider_limits.png
        :alt: Basic plot example with a second spectrum overlaid in green plus adjusted limits
        :figwidth: 800
        :width: 800


Basic plot example with a second spectrum offset and overlaid in red, again with adjusted limits:


.. figure:: images/basic_plot_example_with_second_spectrum_offset_overlaid_in_red.png
        :alt: Basic plot example with a second spectrum offset and overlaid in red, again with adjusted limits
        :figwidth: 800
        :width: 800



API Documentation for Plotting
------------------------------

We include the API documentation for the generic model and fitter wrappers here.

.. automodule:: pyspeckit.spectrum.plotters
    :members:
