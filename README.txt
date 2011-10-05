Spectroscopic Toolkit
---------------------

This is a code framework designed to allow for analysis of spectroscopic data
from a wide variety of astronomical instruments.  It is motivated by the lack
of general spectroscopic analysis tools applicable at multiple wavelengths.

Initial implementation focuses on optical and radio applications, e.g.
gaussian and voigt profile fitting, baseline/continuum fitting, and equivalent
width measurements.  However, the design is meant to be extensible.  We want
additional features to be trivial to implement.

In that vein, there is a growing set of spectral models implemented, primarily
of hyperfine components of radio lines.  The "hyperfine" model class makes
model implementation quite straightforward, if not trivial (you still have to
plug in the right frequency offsets and relative line strengths).

Requirements:
matplotlib
numpy
mpfit (http://code.google.com/p/astrolibpy/source/browse/trunk/mpfit)

At least one of:
pyfits (http://www.stsci.edu/resources/software_hardware/pyfits/Download)
atpy (http://atpy.github.com/)
asciitable (http://cxc.harvard.edu/contrib/asciitable/)

Optional:
scipy 

A good model for our code is:
http://atpy.github.com/developers.html
from which we derived inspiration for the "register functions" capability.

Authors:
Adam Ginsburg (adam.g.ginsburg@gmail.com)
Jordan Mirocha (mirochaj@gmail.com)
