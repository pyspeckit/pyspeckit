PySpecKit Spectroscopic Toolkit
-------------------------------

* Documentation: `https://pyspeckit.readthedocs.io/en/latest/`_
* Source pages: `<https://github.com/pyspeckit/pyspeckit>`_

This is a code framework designed to allow for analysis of spectroscopic data
from a wide variety of astronomical instruments.  It is motivated by the lack
of general spectroscopic analysis tools applicable at multiple wavelengths
(compare to IRAF, SPLAT, etc. - these are wavelength-specific and/or do not
make user scripting easy).

Initial implementation focuses on optical and radio applications, e.g.
gaussian and voigt profile fitting, baseline/continuum fitting, and equivalent
width measurements.  However, the design is meant to be extensible.  We want
additional features to be trivial to implement.

In that vein, there is a growing set of spectral models implemented.  The
model and hyperfinemodel classes makes
model implementation quite straightforward, if not trivial (you still have to
plug in the right frequency offsets and relative line strengths).

Plotting is straightforward, as is usually the case with matplotlib-based
codes.  We have a few different methods of error bar plotting implemented, a
decent (and expanding) units class for pretty printing of spectroscopic units,
and different methods of fit plotting.

We're also looking for more users to give us more use cases!  Contact us
at pyspeckit@gmail.com or post `issues
<https://github.com/pyspeckit/pyspeckit/issues>`_


Requirements:
`astropy <http://www.astropy.org>`_
`matplotlib <http://matplotlib.org/>`_
`numpy <http://numpy.org/>`_

Optional:
`atpy <http://atpy.readthedocs.org/>`_
`scipy <http://www.scipy.org/>`_
`lmfit <https://github.com/lmfit/lmfit-py>`_
`lineid_plot <https://pythonhosted.org/lineid_plot/>`_
`spectral-cube <http://spectral-cube.readthedocs.io/>`_

Authors:
`Adam Ginsburg <adam.g.ginsburg@gmail.com>`_
`Jordan Mirocha <mirochaj@gmail.com>`_
(or both of us at pyspeckit@gmail.com)

Contributors: (see https://github.com/pyspeckit/pyspeckit/graphs/contributors)

 * Erik Rosolowsky (Ammonia models, RADEX-based models)
 * Vlas Sokolov (cube fitting, ammonia and N2H+ modeling)
 * Miguel de Val-Borro (CLASS file reading, python3 compatibility, bugfixes)
 * Brigitta Sipocz (internals & logistics)
 * Jaime Pineda (N2D+, N2H+ models, bugfixes)
 * Allison Youngblood (H2 model fit)
 * Taylor Hogge (Ammonia models)
 * Dinos Kousidis (ESO summer student - Astropy integration)
 * Mike Lum (EQW fitter)
 * Matt Craig, Erik Tollerud, Thomas Robitaille (minor - Astropy integration)


The PySpecKit logo uses the Voyager 1 image of Earth known as the "Pale Blue Dot"
[ `original source <http://visibleearth.nasa.gov/view_rec.php?id=601>`_ |  `reprocessed image <http://instructors.dwrl.utexas.edu/mcginnis/sites/instructors.cwrl.utexas.edu.mcginnis/files/pale_blue_dot2.jpg>`_ ]


.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.12490.svg
   :target: https://doi.org/10.5281/zenodo.12490

