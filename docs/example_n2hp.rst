
.. include:: <isogrk3.txt>

Radio Fitting: N\ :sub:`2`\ H+ example
============================================
Example hyperfine line fitting for the N\ :sub:`2`\ H+ 1-0 line.


:: 

    import pyspeckit

    # Load the spectrum
    sp = pyspeckit.Spectrum('n2hp_opha_example.fits')

    # Register the fitter
    # The N2H+ fitter is 'built-in' but is not registered by default; this example
    # shows how to register a fitting procedure
    # 'multi' indicates that it is possible to fit multiple components and a
    # background will not automatically be fit 
    # 4 is the number of parameters in the model (excitation temperature,
    # optical depth, line center, and line width)
    sp.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter,4)

    # Run the fitter
    sp.specfit(fittype='n2hp_vtau',guesses=[15,2,4,0.2])

    # Plot the results
    sp.plotter()
    # Re-run the fitter (to get proper error bars) and show the individual fit components
    sp.specfit(fittype='n2hp_vtau', guesses=[15,2,4,0.2], show_hyperfine_components=True)

    # Save the figure (this step is just so that an image can be included on the web page)
    sp.plotter.savefig('n2hp_ophA_fit.png')

.. figure:: images/n2hp_ophA_fit.png
	:alt: Fit to the 15 hyperfine components of N2H+ simultaneously.  The (0)'s indicate that this is the 0'th velocity component being fit
        :figwidth: 650
        :width: 650
