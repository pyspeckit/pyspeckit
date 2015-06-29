.. pyspeckit documentation master file, created by
   sphinx-quickstart on Thu Jun 30 16:32:03 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySpecKit
=========
An `extensible <registration.html>`_ spectroscopic analysis toolkit for astronomy.

If you're just getting started, see the :doc:`examples page <examples>`.

To cite pyspeckit, use http://adsabs.harvard.edu/abs/2011ascl.soft09001G.

Downloads
^^^^^^^^^

 * `March 2015 version <https://github.com/pyspeckit/pyspeckit/archive/pyspeckit_0.1.15.tar.gz>`_
 * `latest commit from bitbucket <https://bitbucket.org/pyspeckit/pyspeckit/get/master.tar.gz>`_ (see :doc:`install`) 
 * `latest commit from github <https://github.com/pyspeckit/pyspeckit/archive/master.zip>`_ (same as above)
 * `pypi entry <http://pypi.python.org/pypi/pyspeckit>`_.


Supported file types and their formats:

    * :doc:`fitsfiles` 
    * :doc:`txtfiles` 
    * :doc:`hdf5files` 

Guides / Getting Started
========================
* :doc:`guide_class`
    A simple getting started guide aimed at Gildas-CLASS users

* :doc:`guide_iraf`
    Intended for users of IRAF's `splot` interactive fitting routine.

Classes and API
===============
At the core, PySpecKit runs on a 'Spectroscopic Object' class called
:class:`Spectrum <pyspeckit.spectrum.classes.Spectrum>`.  Therefore everything
interesting about PySpecKit can be learned by digging into the properties of
this class.

* :doc:`spectrum` can read a variety of individual spectra types

  + :class:`~pyspeckit.spectrum.classes.Spectrum` The Spectrum class, which is the core
    of pyspeckit.  The ``__init__`` procedure opens a spectrum file.
  + :class:`~pyspeckit.spectrum.classes.Spectra` A group of ``Spectrum`` s.  Generally for
    when you have multiple wavelength observations you want to stitch together
    (e.g., two filterbanks on a heterodyne system, or the red/blue spectra from
    a multi-band spectrometer like the Double Imaging Spectrograph)
  + :class:`~pyspeckit.spectrum.classes.ObsBlock` An Observation Block - multiple
    spectra of different objects or different times covering the same
    wavelength range

* :doc:`cubes` is used to deal with data cubes and has functionality similar to
  `GAIA <http://astro.dur.ac.uk/~pdraper/gaia/gaia.html>`_ and `ds9
  <http://hea-www.harvard.edu/RD/ds9/>`_.

  + :class:`~pyspeckit.cubes.SpectralCube.Cube` A Cube of Spectra.  Has features to collapse
    the cube along the spectral axis and fit spectra to each element of the
    cube.  Is meant to replicate `Starlink's GAIA
    <http://astro.dur.ac.uk/~pdraper/gaia/gaia.html>`_ in some ways, but with
    less emphasis on speed and much greater emphasis on spectral line fitting.

Features
========

 * :doc:`baseline` describes baseline & continuum fitting.
 * :doc:`fitting` describes the general process of model fitting.
 * :doc:`measurements` is a toolkit for performing EQW, column, and other measurements...
 * :doc:`units` contains the all-important
   :class:`~pyspeckit.spectrum.units.SpectroscopicAxis` class that is used to deal with
   coordinate transformations
 * :doc:`registration` describes the extensible qualities of pyspeckit


.. toctree::
   :maxdepth: 3
   :hidden:

   Index <self>
   install
   Models & Fitting <models>
   Class Features <features>
   File Readers <readers>
   High-Level Wrappers <wrappers>
   Examples <examples>
   Interactive Use <interactive>

..   Base Classes <classes>
