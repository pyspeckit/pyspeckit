Fitting using a Template
------------------------
Pyspeckit allows you to use a spectral template as the model to fit to your spectrum.
See `pyspeckit.spectrum.models.template`.

If your model spectrum only requires a shift and a scale, it's easy to use:

.. code-block:: python

    from pyspeckit.spectrum.models.template import template_fitter

    template = pyspeckit.Spectrum('template_spectrum.fits')
    dataspec = pyspeckit.Spectrum("DataSpectrum.fits")

    # Create the fitter from the template spectrum and "Register" it
    template_fitter = template_fitter(template,xshift_units='angstroms')
    dataspec.Registry.add_fitter('template',template_fitter,2)

    # The fitted parameters are amplitude & xshift
    # perform the fit:
    dataspec.specfit(fittype='template',guesses=[1,0])

    # print the results
    print dataspec.specfit.parinfo
