PySpecKit Spectroscopic Toolkit
-------------------------------
`<http://pyspeckit.bitbucket.org>`_
and `<https://github.com/pyspeckit/pyspeckit>`_

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
<https://bitbucket.org/pyspeckit/pyspeckit.bitbucket.org/issues>`_


Requirements:
`astropy <http://www.astropy.org>`_
`matplotlib <http://matplotlib.sourceforge.net/>`_
`numpy <http://numpy.scipy.org/>`_

Optional:
`atpy <http://atpy.github.com/>`_
`scipy <http://www.scipy.org/>`_
`lmfit <https://github.com/lmfit/lmfit-py>`_
`lineid_plot <http://packages.python.org/lineid_plot/>`_
`spectral-cube <http://spectral-cube.readthedocs.org/>`_

Authors:
`Adam Ginsburg <adam.g.ginsburg@gmail.com>`_
`Jordan Mirocha <mirochaj@gmail.com>`_
(or both of us at pyspeckit@gmail.com)

The PySpecKit logo uses the Voyager 1 image of Earth known as the "Pale Blue Dot"
[ `original source <http://visibleearth.nasa.gov/view_rec.php?id=601>`_ |  `reprocessed image <http://instructors.cwrl.utexas.edu/mcginnis/sites/instructors.cwrl.utexas.edu.mcginnis/files/pale_blue_dot2.jpg>`_ ]


