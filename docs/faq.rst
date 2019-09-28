Frequently Asked Questions
==========================

I see "UnitConversionError: 'km / s' (speed) and 'GHz' (frequency) are not convertible" errors when using the ammonia fitter
----------------------------------------------------------------------------------------------------------------------------

For versions of pyspeckit 0.16 (May 2015) and later, pyspeckit uses astropy's
units for the spectroscopic axis.  It therefore requires an ``equivalency`` to be
defined.

To create a ``SpectroscopicAxis`` with the appropriate equivalency defined, the
axis must have a reference frequency (``refX``) and a ``velocity_convention``.::


    >>> from astropy import units as u
    >>> from pyspeckit.spectrum.units import SpectroscopicAxis
    >>> import numpy as np
    >>> xarr = SpectroscopicAxis(np.linspace(22,24)*u.GHz,
                                 refX=23*u.GHz,
                                 velocity_convention='radio')

If you have loaded a spectrum from file and it doesn't contain the appropriate
metadata (usually a ``CTYPE`` in the header), you can set the ``refX`` and
``velocity_convention`` manually.  The options for ``velocity_convention``
are ``radio``, ``optical``, and ``relativisitic``.

For details on the meaning of the various velocity conventions, see `Frank
Ghigo's site <http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html>`_ or
`FITS Paper III`_,
especially
`table 4 <http://www.aanda.org/articles/aa/full/2006/05/aa3818-05/aa3818-05.html>`_.
For details on the accepted FITS ``CTYPE`` keywords, see `FITS Paper III`_.
In particular, their
`table 1 <http://www.aanda.org/articles/aa/full/2006/05/aa3818-05/aa3818-05.html>`_
specifies all valid spectral coordinate type codes.

.. _FITS Paper III: http://adsabs.harvard.edu/abs/2006A%26A...446..747G
