Guide for IRAF users
====================
PySpecKit is similar in intent and implementation to IRAF's *splot* routine.
IRAF users will probably want to do most of their data reduction (i.e., turning
an image/cube into a 1D wavelength-calibrated spectrum) in IRAF, but will be
comfortable fitting lines and making publication-quality plots using PySpecKit.

Loading a Spectrum
------------------
If you have an IRAF spectrum, it is straightforward to load into PySpecKit::

    sp = pyspeckit.Spectrum('iraf_spectrum.ms.fits')
    sp.plotter()

Fitting Line Profiles
---------------------
.. note:: See :doc:`interactive` for a comprehensive graphical demonstration of these instructions.

In IRAF, a line profile is fitted using *k* to start the fitter, then *k*, *l*, or *v*
to perform the fit.

In PySpecKit, the continuum (*baseline*) and line profile are determined separately.

Instead of using a key twice to specify the continuum level, a continuum must
be fitted from the data.  This is done by pressing *b* to turn on the baseline fitter.
Click or press *1* to select baseline regions - they will be highlighted in green.
Press *3* to fit the baseline and display it as an orange line.

In PySpecKit, the interactive fitter is started by pressing *f* in the plot
window.  After pressing *f*, instructions will be provided in the terminal
window telling you which line profiles are implemented.  Select one of these,
or use a gaussian by default.

Select the line fitting region by pressing *1* on either side of the line.
Select the peak and full-width-half-maximum of the line by pressing *2* at each
of these locations in turn.  You may repeat this cycle indefinitely to fit
multiple profiles (comparable to IRAF*s *deblend* capability).  Then, press *3*
to perform the fit.

The fitted parameters can be accessed (as variables, or printed) through the
`Spectrum.specfit.parinfo` parameter.  They will also be displayed in the plot
legend.

 
..  PySpecKit shares some capabilities with `IRAF <http://iraf.net>`_, but IRAF is
..  a much more extensive tool suite designed to deal with images and 2 or 3
..  dimensional spectra.  PySpecKit is not that - there is no aperture extraction
..  toolkit, no way to trace stellar spectra, and no geometric transforms for images.

