CHANGES
*******

Version 0.1.19 (unreleased)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Version 0.1.18.1 (2015-12-10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    * MAJOR bugfix #125: Discovered a major issue introduced in #122 in which a
      typo in the radiative transfer equation results in an incorrect model

Version 0.1.18 (2015-12-09)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Major change: this is the first version of pyspeckit fully compatible with
python3.

    * Bugfixes: #92, #93 related to fiteach, parameter validation,
      spectroscopic axis array handling, in_range checks, and the ammonia fit
      and plot wrapper.  Also includes an update to the versioning scheme such
      that dev versions (like 0.1.18.dev) will include a git devstring when
      possible.
    * Bugfix #91: Offset of 0.13 km/s between ammonia frequency used and that
      in the LOVAS catalog
    * Bugfix #90: Allow fitting of multiple peaks interactively (again)
    * Enhancement #89: Allow data to be passed to Spectrum as a quantity
    * Bugfix #88: line_ids needed to handle quantities
    * Enhancement #85: Remove hard-coded T_A* label
    * Bugfix #82: Fix a loophole in which parameters were not checked for being
      in bounds when passed as `guesses`
    * Addition #99: fiteach example with NIR data cubes, specfit example using
      user-created model, and a guide to installing pyspeckit via GitHub.
    * Enhancement #100: Display EQW highlighting at appropriate location
    * Addition #102: New N2D+ model
    * Bugfix #105: Ammonia thin works (though it is still not recommended)
    * Bugfix #106: Baselines were dependent on the X-axis coordinate unit.
    * Enhancement #108, 110, 111, 112: Python 3 compatibility
    * Bugfix #113: Improved plotting & robustness of measure_fwhm, especially
      when baseline is not subtracted
    * Bugfix #115: py3ify steppify
    * Enhancement #117/#119: Better error messages when using
      ``SpectralCube.load_model`` and more robust figure (re-)opening (if you
      closed a figure created by spectral cube and then tried to plot again, it
      would result in a crash.  This behavior was never supported, so it is not
      a bugfix, but now the behavior should work)
    * Bugfix #118: Re-fitting within a given window sometimes failed, possibly
      due to how matplotlib handles event handlers.  Refactoring should make
      this much more robust.
    * Enhancement #120: astropy-helpers + travis-ci related infrastructure
      update
    * Bugfix #121: Use correct variable name in N2D+ degeneracy
    * Bugfix #122: Informative errors on file reading & exact RT equation in
      ``hyperfine`` model
    * Addition #124: LTE model

Version 0.1.17 (2015-07-14)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

    * Bugfixes: #67, #69, #71, #74, #75 related to fiteach
    * Bugfix for error bar plotting (PR #76, issue #78)
    * Documentation cleanup and enhancement (#77)

Release 0.1.16 (2015-05-21)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

    * Major refactor: use astropy's Quantity and units to replace pyspeckit's
      minimalist unit framework.  You may see deprecation warnings as a result!
      (https://github.com/pyspeckit/pyspeckit/pull/26)
    * The refactor led to many subsequent bugfixes: #61, #55, #51, and others
    * Removal of the `peakbgfit` default method: instead, the default is to treat
      all fits as multifits.  Changes came from #32, #25
    * New ammonia models (#28, #50): absorption against a background source and
      treating each line independently without a temperature connecting them

Release 0.1.15 (2014-11-09)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Bugfix: write_fit is part of Cube, not CubeStack
    * Bugfix: excludefit must occur after selectregion if fit_plotted_area is True
    * API change: For the fitter & baseliner, data selection is end-inclusive
                  if specified in world coordinates
    * Bugfix: numpy 1.8 added a "writeable" flag that broke units; that is now 
              corrected
              http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.flags.html
    * Baseline: add a `fit` method that is independent from `button2action` to
      make masking and fitting independent processes
    * Added `astropy <http://astropy.org>`_ as a dependency
    * Converted to astropy-helper template for setup
    * Changed 'units' -> 'unit' in Spectrum
    * Baseline: add spline fitting
    * Add Zenodo badge
    * Allow `parinfo` to be used in place of `guesses` in specfit

Release 0.1.14 (2013-09-10)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Bugfix: integral(direct=True) double-subtracted the baseline if it was
              already subtracted
    * New Feature: Models now include analytic integrals (only implemented for
                   Gaussian so far)
    * New Feature: hyperfine models added that allow varying amplitudes & widths

Release 0.1.13 (2013-03-04)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * Cleanup excess data files
    * Bugfixes in baseline fitting
    * Added astropy.models interface

Release 0.1.12
~~~~~~~~~~~~~~
    * New formaldehyde 218 GHz fitter
    * Allow parinfo to be input as "guess=" or "parinfo="

Release 0.1.11
~~~~~~~~~~~~~~
    * bugfix to EQW non-fitted ("empirical")

Release 0.1.10
~~~~~~~~~~~~~~
    * bugfix: unit conversion with reference wavelength
    * bugfix: interactive buttons "reconnected" each time 
    * new feature: voigt profile interactive guess now has 2 widths 

Release 0.1.9 
~~~~~~~~~~~~~
    * Added `lineid_plot <http://packages.python.org/lineid_plot/>`_ tools
    * Baseline can fit power laws
    * New TSPEC unmerged, IRAF fits readers
    * astropy.io.fits compatibility fixes
    * General bugfixes
    * Voigt Profile Fitter - bugfix, previously abused notation / misused widths

Release 0.1.8
~~~~~~~~~~~~~
    * BUGFIX RELEASE
    * Documentation updates
    * Speed boost for model grids
    * Some support for pymodelfit
    * added emcee and pymc generators

Release 0.1.7
~~~~~~~~~~~~~
    
    * Added cross-correlation 
    * Added (some) unit tests
    * New interactive slider widget for manipulating models (Spectrum.specfit.add_sliders())
    * Subtle but very significant bug-fix: parinfo is now a single uniform
      list, based on the ParinfoList class.
    * You can now fit based on what you see by using the 'use_window_limits=True' kwarg
      .. warning:: This changes the default behavior in interactive mode!
    * lmfit-py can now be used for fitting via the 'use_lmfit' kwarg
    * BUGFIX: SpectroscopicAxis can be converted between units even when scalar
    * velocity frames are read from FITS headers following the VELDEF conventions

Release 0.1.6 
~~~~~~~~~~~~~

    * H2CO fit / plot wrapper
    * bugfixes
    * setup.py no longer tries to write config files

Release 0.1.5 
~~~~~~~~~~~~~

    * Added GBT (GBTIDL SDFITS file) and ALFALFA (ALFALFA idlsave .src) readers
    * added extinction model (just a function, not a complete model yet)

Release 0.1.4 
~~~~~~~~~~~~~

    * removed setuptools & distribute (they fail)
    * added hydrogen.py to models
    * first release to inherit from astropy's Spectrum1D

Release 0.1.3 
~~~~~~~~~~~~~

    * some internal cleanup / refactoring
    * override slicing (__getitem__ features)
    * parallel moment & fitting in Cubes repaired

Release 0.1.2 
~~~~~~~~~~~~~

    * added MIT license, moved mpfit and parallel_map inside pyspeckit as
    * submodules

Release 0.1.1 
~~~~~~~~~~~~~

    * bugfixes and versioning work

Release 0.1.0 
~~~~~~~~~~~~~

    * Initial creation
