Model Fitting
=============

General documentation for model fitting is below.  There are several multi-component model fitting tools that have been developed based on pyspeckit:

 * `Mike Chen's multi-vlsr <https://github.com/mcyc/multi_vlsr>`_
 * `Vlas Sokolov's pyspecnest <https://github.com/vlas-sokolov/pyspecnest>`_ and `multicube <https://github.com/vlas-sokolov/multicube>`_
 * `Jonny Henshaw's SCOUSE <https://github.com/jdhenshaw/scousepy>`_
 * `Jared Keown's astroclover <https://github.com/jakeown/astroclover>`_

Model comparison: log-likelihood, AIC, and BIC
----------------------------------------------

After a fit, the `Specfit` object provides the (Gaussian) log-likelihood and
the common information criteria as properties: ``logL``, ``AIC`` (Akaike
information criterion), ``AICc`` (AIC with the small-sample correction), and
``BIC`` (Bayesian information criterion).  These assume normally-distributed,
uncorrelated errors given by the error spectrum:

.. math::

    \ln L = -\frac{1}{2}\chi^2 - \sum_i \ln\left(\sqrt{2\pi}\,\sigma_i\right)

with :math:`\mathrm{AIC} = 2k - 2\ln L` and
:math:`\mathrm{BIC} = k\ln(n) - 2\ln L`, where :math:`k` is the number of
*free* parameters (fixed and tied parameters are not counted; see
``Specfit.n_free_parameters``) and :math:`n` is the number of pixels included
in the fit (``Specfit.npix_fitted``).  Lower values indicate a preferred
model, so these can be used to decide, e.g., how many velocity components a
spectrum contains::

    import numpy as np
    from pyspeckit import Spectrum

    # synthetic two-component spectrum with known noise
    xarr = np.linspace(-6, 6, 240)
    sigma = 0.1
    data = (1.5 * np.exp(-(xarr + 2.5)**2 / (2 * 0.5**2)) +
            1.0 * np.exp(-(xarr - 2.5)**2 / (2 * 0.4**2)) +
            np.random.randn(xarr.size) * sigma)
    sp = Spectrum(xarr=xarr, data=data, error=np.ones(xarr.size)*sigma)

    # one-component fit
    sp.specfit(fittype='gaussian', guesses=[1.2, -2.0, 0.6])
    aic_one, bic_one = sp.specfit.AIC, sp.specfit.BIC

    # two-component fit
    sp.specfit(fittype='gaussian',
               guesses=[1.2, -2.0, 0.6, 0.8, 2.0, 0.5])
    aic_two, bic_two = sp.specfit.AIC, sp.specfit.BIC

    print("Delta-AIC (1comp - 2comp) = {0}".format(aic_one - aic_two))
    print("Delta-BIC (1comp - 2comp) = {0}".format(bic_one - bic_two))
    # both are strongly positive: the two-component model is preferred

If no error spectrum is provided, the errors are estimated from the data and
the log-likelihood (and therefore the absolute value of the information
criteria) is only defined up to an additive constant; differences between
models evaluated on the same spectrum remain meaningful.

.. automodule:: pyspeckit.spectrum.fitters
.. autoclass:: Specfit
    :show-inheritance:
    :members:
    :inherited-members:
    :undoc-members:


.. toctree::
   fit_custom_model
