PySpecKit Projects
==================

A few projects appropriate for a Google Summer of Code, `ESA Summer of Code in
Space <http://sophia.estec.esa.int/socis2013/>`_, or similar are outlined
below.

Incorporate astropy.units into pyspeckit.units
----------------------------------------------
This project is at the core of both `pyspeckit` and `specutils
<https://github.com/astropy/specutils>`_.

The most important base functionality of a spectroscopic toolkit
is to be able to convert between different spectroscopic systems,
e.g. wavelength<->frequency<->velocity.  This can be achieved
using astropy_'s `unit equivalencies
<https://github.com/astropy/astropy/pull/1176>`_.  

The X-axis unit changes will be straightforward project that should require
about 2 weeks to complete.  The more complicated and interesting project is
creating Y-axis units (i.e., flux units) that appropriately adjust with changes
to the X-axis.  These would make use of other astropy_ unit equivalencies,
e.g.  the `spectral density
<https://github.com/astropy/astropy/blob/master/astropy/units/equivalencies.py#L44>`_
equivalency.

The end goal will be to have a `Spectrum` object that will live in `specutils`_
and be inherited by `pyspeckit`, which will provide the interface to modeling
and graphical tools.

Refactor and Expand the pyspeckit modeling tools
------------------------------------------------
Since the development of pyspeckit, there has been substantial progress on a
more general class of `modeling tools from astropy
<http://docs.astropy.org/en/latest/modeling/index.html>`_.
Pyspeckit already has a wide variety of data fitting and modeling tools that
can readily be modified to use the astropy modeling formalism.

Details of this project need to be worked out, but will include:

 * refactoring pyspeckit.models to use astropy.models
 * building a graphical interface to astropy.models

Expand the Unit Test suite
--------------------------

Pyspeckit is a complicated code suite, which has led to many bugs, particularly
in the UI.  An improved unit test suite would help prevent or remove these
bugs.  Such a project would start by breaking down the existing tests, which
are really end-to-end tests, into their component units.
