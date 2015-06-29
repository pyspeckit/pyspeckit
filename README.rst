PySpecKit Spectroscopic Toolkit
-------------------------------

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/pyspeckit/pyspeckit
   :target: https://gitter.im/pyspeckit/pyspeckit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

 * Documentation: `<http://pyspeckit.readthedocs.org/en/latest/>`_
 * Homepage: `<pyspeckit.bitbucket.org>`_
 * Source pages: `<https://github.com/pyspeckit/pyspeckit>`_ and `<https://bitbucket.org/pyspeckit/pyspeckit>`_

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
<https://bitbucket.org/pyspeckit/pyspeckit/issues>`_


Requirements:
`matplotlib <http://matplotlib.sourceforge.net/>`_
`numpy <http://numpy.scipy.org/>`_

Preferably at least one of:
`pyfits <http://www.stsci.edu/resources/software_hardware/pyfits/Download>`_
`atpy <http://atpy.readthedocs.org/>`_
`asciitable <http://cxc.harvard.edu/contrib/asciitable/>`_

Optional:
`scipy <http://www.scipy.org/>`_
`lmfit-py <https://github.com/newville/lmfit-py>`_
`lineid_plot <http://packages.python.org/lineid_plot/>`_

Authors:
`Adam Ginsburg <adam.g.ginsburg@gmail.com>`_
`Jordan Mirocha <mirochaj@gmail.com>`_
(or both of us at pyspeckit@gmail.com)

The PySpecKit logo uses the Voyager 1 image of Earth known as the "Pale Blue Dot"
[ `original source <http://visibleearth.nasa.gov/view_rec.php?id=601>`_ |  `reprocessed image <http://instructors.dwrl.utexas.edu/mcginnis/sites/instructors.cwrl.utexas.edu.mcginnis/files/pale_blue_dot2.jpg>`_ ]


.. image:: https://zenodo.org/badge/6907/pyspeckit/pyspeckit.png
    :target: http://dx.doi.org/10.5281/zenodo.12490


