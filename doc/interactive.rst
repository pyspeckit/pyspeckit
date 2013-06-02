A guide to interactive fitting
==============================

A step-by-step example of how to use the interactive fitter.

In short, we will do the following:

.. code-block:: python

    # 1. Load the spectrum
    sp = pyspeckit.Spectrum('hr2421.fit')

    # 2. Plot a particular spectral line
    sp.plotter(xmin=4700,xmax=5000)

    # 3. Need to fit the continuum first
    sp.baseline(interactive=True, subtract=False)

    # 4... (much work takes place interactively at this stage)

    # 5. Start up an interactive line-fitting session
    sp.specfit(interactive=True)

.. note:: 

    If you don't see a plot window after step #2 above, make sure you're using
    matplotlib in interactive mode.  This may require starting ipython as
    ``ipython --pylab``

.. figure:: images/interactive_example_hr2421_baseline_firstclick.png
	:alt: First click fitting a baseline
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_secondclick.png
	:alt: Second click fitting a baseline
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_secondclick_highlight.png
	:alt: The results of the second click
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_thirdclick.png
	:alt: Third click fitting a baseline
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_fourthclick.png
	:alt: Fourth click fitting a baseline 
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_fourthclick_highlight.png
	:alt: The results of the fourth click ("exclude")
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_baseline_fifthclick_fit.png
	:alt: Fifth click fitting a baseline (run the fit)
        :figwidth: 800
        :width: 800

This is where you start the line-fitter:

.. code-block:: python

    # Start up an interactive line-fitting session
    sp.specfit(interactive=True)

.. figure:: images/interactive_example_hr2421_firstclick.png
	:alt: First click fitting a spectral line
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_secondclick.png
	:alt: Second click fitting a spectral line
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_secondclick_highlight.png
	:alt: Results of the second click (highight the fit region)
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_thirdclick.png
	:alt: Third click fitting a spectral line
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_fourthclick.png
	:alt: Fourth click fitting a spectral line
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_gaussmodelguess.png
	:alt: Results of the fourth click: make a gaussian guess
        :figwidth: 800
        :width: 800

.. figure:: images/interactive_example_hr2421_fifthclick_fit.png
	:alt: Fifth click fitting a spectral line - do the fit
        :figwidth: 800
        :width: 800

