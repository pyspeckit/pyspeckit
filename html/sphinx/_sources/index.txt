.. pyspeckit documentation master file, created by
   sphinx-quickstart on Thu Jun 30 16:32:03 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySpecKit
=====================================
An `extensible <registration.html>`_ spectroscopic analysis toolkit for astronomy.

If you're just getting started, see the `examples page
<../examples.html>`_.

`Download <../../get/tip.tar.gz>`_ or see :doc:`install`


Supported file types and their formats:

    * :doc:`fitsfiles` 
    * :doc:`txtfiles` 
    * :doc:`hdf5files` 


Classes
=====================================

* :doc:`spectrum` can read a variety of individual spectra types

  * :class:`~pyspeckit.spectrum.Spectrum` The Spectrum class, which is the core
    of pyspeckit.  The __init__ procedure opens a spectrum file.
  * :class:`~pyspeckit.spectrum.Spectra` A group of Spectrums.  Generally for
    when you have multiple wavelength observations you want to stitch together
    (e.g., two filterbanks on a heterodyne system, or the red/blue spectra from
    a multi-band spectrometer like the Double Imaging Spectrograph)
  * :class:`~pyspeckit.spectrum.ObsBlock` An Observation Block - multiple
    spectra of different objects or different times covering the same
    wavelength range

* :doc:`cubes` is used to deal with data cubes and has functionality similar to
  `GAIA <http://astro.dur.ac.uk/~pdraper/gaia/gaia.html>`_ and `ds9
  <http://hea-www.harvard.edu/RD/ds9/>`_.

  * :class:`~pyspeckit.cubes.Cube` A Cube of Spectra.  Has features to collapse
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
  :class:`~pyspeckit.units.SpectroscopicAxis` class that is used to deal with
  coordinate transformations
* :doc:`registration` describes the extensible qualities of pyspeckit


.. toctree::
   :maxdepth: 2
   :hidden:

   Index <self>
   install
   Models & Fitting <models>
   Base Classes <classes>
   Class Features <features>
   File Readers <readers>
   High-Level Wrappers <wrappers>
