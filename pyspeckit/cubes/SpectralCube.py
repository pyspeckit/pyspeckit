"""
Cubes
=====

Tools to deal with spectroscopic data cubes.

Some features in Cubes require additional packages:

   * smoothing - requires agpy_\'s smooth and parallel_map routines
   * `pyregion <git://github.com/astropy/pyregion.git>`_


The 'grunt work' is performed by the :py:mod:`cubes` module


"""
from __future__ import print_function

import time
import sys
import traceback
import numpy as np
import types
import copy
import itertools
from ..specwarnings import warn,PyspeckitWarning

import astropy
from astropy.io import fits
from astropy import log
from astropy import wcs
from astropy import units
from astropy.utils.console import ProgressBar
from six import iteritems, string_types
from functools import wraps

# import parent package
from .. import spectrum
from ..spectrum import smooth
from ..spectrum.units import (generate_xarr, SpectroscopicAxis,
                              SpectroscopicAxes)
from ..parallel_map import parallel_map
from ..spectrum import history

# import local things
from . import mapplot
from . import cubes

def not_for_cubes(func):

    @wraps(func)
    def wrapper(*args):
        warn("This operation ({0}) operates on the spectrum selected "
             "from the cube, e.g. with `set_spectrum` or `set_apspec`"
             ", it does not operate on the whole cube.", PyspeckitWarning)
        return func(*args)
    return wrapper

class Cube(spectrum.Spectrum):

    def __init__(self, filename=None, cube=None, xarr=None, xunit=None,
                 errorcube=None, header=None, x0=0, y0=0,
                 maskmap=None,
                 **kwargs):
        """
        A pyspeckit Cube object.  Can be created from a FITS file on disk or
        from an array or a `spectral_cube.SpectralCube` object.  If an array
        is used to insantiate the cube, the `xarr` keyword must be given,
        specifying the X-axis units

        Parameters
        ----------
        filename : str, optional
            The name of a FITS file to open and read from.  Must be 3D
        cube : `np.ndarray`, `spectral_cube.SpectralCube`, or \
               `astropy.units.Quantity`
            The data from which to instantiate a Cube object.  If it is
            an array or an astropy Quantity (which is an array with attached
            units), the X-axis must be specified.  If this is given as a
            SpectralCube object, the X-axis and units should be handled
            automatically.
        xarr : `np.ndarray` or `astropy.units.Quantity`, optional
            The X-axis of the spectra from each cube.  This actually
            corresponds to axis 0, or what we normally refer to as the Z-axis
            of the cube, but it indicates the X-axis in a plot of intensity vs
            wavelength.  The units for this array are specified in the `xunit`
            keyword unless a `~astropy.units.Quantity` is given.
        xunit : str, optional
            The unit of the ``xarr`` array if ``xarr`` is given as a numpy
            array
        errorcube : `np.ndarray`, `spectral_cube.SpectralCube`,\
                    or `~astropy.units.Quantity`, optional
            A cube with the same shape as the input cube providing the 1-sigma
            error for each voxel.  This can be specified more efficiently as an
            error map for most use cases, but that approach has not yet been
            implemented.  However, you can pass a 2D error map to `fiteach`.
        header : `fits.Header` or dict, optional
            The header associated with the data.  Only needed if the cube is
            given as an array or a quantity.
        x0, y0 : int
            The initial spectrum to use.  The `Cube` object can be treated as
            a `pyspeckit.Spectrum` object, with all the associated tools
            (plotter, fitter) using the `set_spectrum` method to select a pixel
            from the cube to plot and fit.  However, it is generally more sensible
            to extract individual spectra and treat them separately using the
            `get_spectrum` method, so these keywords MAY BE DEPRECATED in the
            future.
        maskmap : `np.ndarray`, optional
            A boolean mask map, where ``True`` implies that the data are good.
            This will be used for both plotting using `mapplot` and fitting
            using `fiteach`.
        kwargs : dict
            passed to `load_fits`, which are passed to `SpectralCube.read`

        """

        if filename is not None:
            self.load_fits(filename, **kwargs)
            return
        else:
            if hasattr(cube, 'spectral_axis'):
                # Load from a SpectralCube instance
                self.cube = cube.hdu.data
                if (cube.unit in ('undefined', units.dimensionless_unscaled)
                    and 'BUNIT' in cube._meta):
                    self.unit = cube._meta['BUNIT']
                else:
                    self.unit = cube.unit
                log.debug("Self.unit: {0}".format(self.unit))
                if xarr is None:
                    if cube.spectral_axis.flags['OWNDATA']:
                        xarr = SpectroscopicAxis(cube.spectral_axis,
                                                 unit=cube.spectral_axis.unit,
                                                 refX=cube.wcs.wcs.restfrq,
                                                 refX_unit='Hz')
                    else:
                        xarr = SpectroscopicAxis(cube.spectral_axis.copy(),
                                                 unit=cube.spectral_axis.unit,
                                                 refX=cube.wcs.wcs.restfrq,
                                                 refX_unit='Hz')
                if header is None:
                    header = cube.header
            elif hasattr(cube, 'unit'):
                self.cube = cube.value
                self.unit = cube.unit
            else:
                self.cube = cube

            if hasattr(errorcube, 'spectral_axis'):
                # Load from a SpectralCube instance
                self.errorcube = errorcube.hdu.data
            elif hasattr(errorcube, 'unit'):
                self.errorcube = errorcube.value
            else:
                self.errorcube = errorcube
            if hasattr(xarr, 'flags'):
                log.debug("XARR flags: {0}".format(xarr.flags))
            self.xarr = generate_xarr(xarr, unit=xunit)
            if hasattr(xarr, 'flags'):
                log.debug("self.xarr flags: {0}".format(xarr.flags))
            self.header = header
            self.error = None
            if self.cube is not None:
                self.data = self.cube[:,int(y0),int(x0)]

        if not hasattr(self, '_unit'):
            self.unit = units.dimensionless_unscaled
        log.debug("Self.unit before header: {0}".format(self.unit))
        if self.header is not None:
            self.parse_header(self.header)
        else:
            log.debug("self.header is None: {0}".format(self.header))
            self.unit = 'undefined'
            self.header = fits.Header()
        log.debug("Self.unit after header: {0}".format(self.unit))

        if maskmap is not None:
            if maskmap.ndim != 2:
                raise ValueError("Mask map must be two-dimensional.")
            self.maskmap = maskmap
        else:
            self.maskmap = np.ones(self.cube.shape[1:],dtype='bool')

        if isinstance(filename,str):
            self.fileprefix = filename.rsplit('.', 1)[0]    # Everything prior to .fits or .txt
        else:
            self.fileprefix = "pyfitsHDU"

        self.plotter = spectrum.plotters.Plotter(self)
        self._register_fitters()
        self.specfit = spectrum.fitters.Specfit(self,Registry=self.Registry)
        self.baseline = spectrum.baseline.Baseline(self)
        self.speclines = spectrum.speclines
        # Initialize writers
        self.writer = {}
        for writer in spectrum.writers.writers:
            self.writer[writer] = spectrum.writers.writers[writer](self)

        # Special.  This needs to be modified to be more flexible; for now I need it to work for nh3
        self.plot_special = None
        self.plot_special_kwargs = {}
        self._modelcube = None

        if self.header:
            self.wcs = wcs.WCS(self.header)
            self.wcs.wcs.fix()
            self._spectral_axis_number = self.wcs.wcs.spec+1
            self._first_cel_axis_num = np.where(self.wcs.wcs.axis_types // 1000 == 2)[0][0]+1

            # TODO: Improve this!!!
            self.system = ('galactic'
                           if ('CTYPE{0}'.format(self._first_cel_axis_num)
                               in self.header and 'GLON' in
                               self.header['CTYPE{0}'.format(self._first_cel_axis_num)])
                           else 'celestial')
        else:
            self._spectral_axis_number = 2
            self._first_cel_axis_num = 0
            self.system = 'PIXEL'

        self.mapplot = mapplot.MapPlotter(self)

    def load_fits(self, fitsfile, **kwargs):
        try:
            from spectral_cube import SpectralCube
        except ImportError:
            raise ImportError("Could not import spectral_cube.  As of pyspeckit"
                              " 0.17, spectral_cube is required for cube reading. "
                              "It can be pip installed or acquired from "
                              "spectral-cube.rtfd.org.")
        mycube = SpectralCube.read(fitsfile, **kwargs)
        return self.load_spectral_cube(mycube)

    def load_spectral_cube(self, cube):
        """
        Load the cube from a spectral_cube.SpectralCube object
        """
        self.__init__(cube=cube)

    def __repr__(self):
        return (r'<Cube object over spectral range %6.5g :'
                ' %6.5g %s and flux range = [%2.1f, %2.1f]'
                ' %s with shape %r at %s>' %
                (self.xarr.min().value, self.xarr.max().value, self.xarr.unit,
                 self.data.min(), self.data.max(), self.unit, self.cube.shape,
                 str(hex(self.__hash__()))))


    def copy(self,deep=True):
        """
        Create a copy of the spectral cube with its own plotter, fitter, etc.
        Useful for, e.g., comparing smoothed to unsmoothed data
        """
        newcube = copy.copy(self)
        newcube.header = copy.copy(self.header)

        deep_attr_lst = ['xarr', 'data', 'cube', 'maskmap',
                         'error', 'errorcube']
        if deep:
            for attr in deep_attr_lst:
                setattr(newcube, attr, copy.copy(getattr(self, attr)))
            if hasattr(self, 'wcs'):
                newcube.wcs = self.wcs.deepcopy()
            newcube.header = self.header.copy()

        newcube.plotter = self.plotter.copy(parent=newcube)
        newcube._register_fitters()
        newcube.specfit = self.specfit.copy(parent=newcube)
        newcube.specfit.Spectrum.plotter = newcube.plotter
        newcube.baseline = self.baseline.copy(parent=newcube)
        newcube.baseline.Spectrum.plotter = newcube.plotter

        newcube.mapplot = self.mapplot.copy(parent=newcube)
        newcube.mapplot.Cube = newcube

        return newcube

    def _update_header_from_xarr(self):
        """Uses SpectroscopiAxis' _make_header method to update Cube header"""
        self.header['NAXIS3'] = self.xarr.size

        self.xarr._make_header()
        sp_naxis = self._spectral_axis_number

        # change keywords in xarr._make_header from, e.g., CRPIX1 to CRPIX3
        newhead = {(key.replace('1', str(sp_naxis))
                    if key.endswith('1') else key): val
                   for key, val in iteritems(self.xarr.wcshead)}

        for key, val in iteritems(newhead):
            if isinstance(val, units.Quantity):
                newhead[key] = val.value
            elif (isinstance(val, units.CompositeUnit)
                  or isinstance(val, units.Unit)):
                newhead[key] = val.to_string()
            log.debug("Updating header: {}: {}".format(key, val))

        self.header.update(newhead)

    def slice(self, start=None, stop=None, unit='pixel', preserve_fits=False,
              copy=True, update_header=False):
        """
        Slice a cube along the spectral axis
        (equivalent to "spectral_slab" from the spectral_cube package)

        Parameters
        ----------
        start : numpy.float or int
            start of slice
        stop : numpy.float or int
            stop of slice
        unit : str
            allowed values are any supported physical unit, 'pixel'
        update_header : bool
            modifies the header of the spectral cube according to the slice
        """

        x_in_units = self.xarr.as_unit(unit)
        start_ind = x_in_units.x_to_pix(start)
        stop_ind = x_in_units.x_to_pix(stop)
        if start_ind > stop_ind:
            start_ind, stop_ind = stop_ind, start_ind
        spectrum_slice = slice(start_ind,stop_ind)

        if not copy:
            raise NotImplementedError("Must copy when slicing a cube.")
        newcube = self.copy()
        newcube.cube = newcube.cube[spectrum_slice,:,:]
        if hasattr(newcube,'errcube'):
            newcube.errcube = newcube.errcube[spectrum_slice,:,:]
        newcube.data = newcube.data[spectrum_slice]
        if newcube.error is not None:
            newcube.error = newcube.error[spectrum_slice]
        newcube.xarr = newcube.xarr[spectrum_slice]

        # create new specfit / baseline instances (otherwise they'll be the wrong length)
        newcube._register_fitters()
        newcube.baseline = spectrum.baseline.Baseline(newcube)
        newcube.specfit = spectrum.fitters.Specfit(newcube,Registry=newcube.Registry)

        if preserve_fits:
            newcube.specfit.modelpars = self.specfit.modelpars
            newcube.specfit.parinfo = self.specfit.parinfo
            newcube.baseline.baselinepars = self.baseline.baselinepars
            newcube.baseline.order = self.baseline.order

        # modify the header in the new cube
        if update_header:
            newcube._update_header_from_xarr()
            # create a new wcs instance from the updated header
            newcube.wcs = wcs.WCS(newcube.header)
            newcube.wcs.wcs.fix()
            newcube._spectral_axis_number = newcube.wcs.wcs.spec + 1
            newcube._first_cel_axis_num = np.where(newcube.wcs.wcs.axis_types
                                                   // 1000 == 2)[0][0] + 1

        return newcube


    def __getitem__(self, indx):
        """
        If [] is used on a cube, slice on the cube and use
        the first dimension to slice on the xarr and the data
        """

        return Cube(xarr=self.xarr.__getitem__(indx[0]),
                    cube=self.cube[indx],
                    errorcube=self.errorcube[indx] if self.errorcube else None,
                    maskmap=self.maskmap[indx[1:]] if self.maskmap is not None else None,
                    header=self.header
                   )

    def set_spectrum(self, x, y):
        self.data = self.cube[:,int(y),int(x)]
        if self.errorcube is not None:
            self.error = self.errorcube[:,int(y),int(x)]

    def plot_spectrum(self, x, y, plot_fit=False, **kwargs):
        """
        Fill the .data array with a real spectrum and plot it
        """

        self.set_spectrum(x,y)

        if self.plot_special is None:
            self.plotter(**kwargs)
            if plot_fit:
                self.plot_fit(x,y)
            self.plotted_spectrum = self
        else:
            sp = self.get_spectrum(x,y)
            sp.plot_special = types.MethodType(self.plot_special, sp)
            combined_kwargs = dict(kwargs.items())
            combined_kwargs.update(self.plot_special_kwargs)
            self._spdict = sp.plot_special(**combined_kwargs)
            self.plotted_spectrum = sp
            self.plotter = sp.plotter
            self.plotter.refresh = lambda: [spi.plotter.refresh()
                                            for spi in self._spdict.values()]
            self.specfit.modelplot = [comp
                                      for spi in self._spdict.values()
                                      for comp in spi.specfit.modelplot]
            self.specfit._plotted_components = [comp
                                                for spi in self._spdict.values()
                                                for comp in spi.specfit._plotted_components]

    def plot_fit(self, x, y, silent=False, **kwargs):
        """
        If fiteach has been run, plot the best fit at the specified location

        Parameters
        ----------
        x : int
        y : int
            The x, y coordinates of the pixel (indices 2 and 1 respectively in
            numpy notation)
        """
        if not hasattr(self,'parcube'):
            if not silent:
                log.info("Must run fiteach before plotting a fit.  "
                         "If you want to fit a single spectrum, "
                         "use plot_spectrum() and specfit() directly.")
            return

        if self.plot_special is not None:
            # don't try to overplot a fit on a "special" plot
            # this is already handled in plot_spectrum
            return

        if not self.has_fit[int(y), int(x)]:
            # no fit to plot
            return

        self.specfit.modelpars = self.parcube[:,int(y),int(x)]
        if np.any(np.isnan(self.specfit.modelpars)):
            log.exception("Attempted to plot a model with NaN parameters.")
            return

        self.specfit.npeaks = self.specfit.fitter.npeaks
        self.specfit.model = self.specfit.fitter.n_modelfunc(self.specfit.modelpars,
                                                             **self.specfit.fitter.modelfunc_kwargs)(self.xarr)

        # set the parinfo values correctly for annotations
        self.specfit.parinfo.values = self.parcube[:,int(y),int(x)]
        self.specfit.parinfo.errors = self.errcube[:,int(y),int(x)]
        self.specfit.fitter.parinfo.values = self.parcube[:,int(y),int(x)]
        self.specfit.fitter.parinfo.errors = self.errcube[:,int(y),int(x)]
        #for pi,p,e in zip(self.specfit.parinfo,
        #                  self.specfit.modelpars,
        #                  self.errcube[:,int(y),int(x)]):
        #    try:
        #        pi['value'] = p
        #        pi['error'] = e
        #    except ValueError:
        #        # likely to happen for failed fits
        #        pass

        self.specfit.plot_fit(**kwargs)

    def plot_apspec(self, aperture, coordsys=None, reset_ylimits=True,
                    wunit='arcsec',
                    method='mean', **kwargs):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube coordinates)

        Parameters
        ----------
        aperture : list
            A list of aperture parameters, e.g.
             * For a circular aperture, len(ap)=3:
               + ``ap = [xcen,ycen,radius]``
             * For an elliptical aperture, len(ap)=5:
               + ``ap = [xcen,ycen,height,width,PA]``
        coordsys : None or str
            The coordinate system of the aperture (e.g., galactic, fk5, None
            for pixel)
        method : 'mean' or 'sum'
            Either average over parellel spectra or sum them.
        """


        if self.plot_special is None:
            self.set_apspec(aperture, coordsys=coordsys, method=method)
            self.plotter(reset_ylimits=reset_ylimits, **kwargs)
        else:
            #self.plot_special(reset_ylimits=reset_ylimits, **dict(kwargs.items()+self.plot_special_kwargs.items()))

            sp = self.get_apspec(aperture, coordsys=coordsys, wunit=wunit, method=method)
            sp.plot_special = types.MethodType(self.plot_special, sp)
            combined_kwargs = dict(kwargs.items())
            combined_kwargs.update(self.plot_special_kwargs)
            sp.plot_special(reset_ylimits=reset_ylimits, **combined_kwargs)

    def get_spectrum(self, x, y):
        """
        Very simple: get the spectrum at coordinates x,y

        (inherits fitter from self)

        Returns a SpectroscopicAxis instance
        """

        ct = 'CTYPE{0}'.format(self._first_cel_axis_num)
        header = cubes.speccen_header(fits.Header(cards=[(k,v) for k,v in
                                                         iteritems(self.header)
                                                         if k != 'HISTORY']),
                                      lon=x, lat=y, system=self.system,
                                      proj=(self.header[ct][-3:]
                                            if ct in self.header else
                                            'CAR'))

        sp = spectrum.Spectrum(xarr=self.xarr.copy(), data=self.cube[:,int(y),int(x)],
                               header=header, error=(self.errorcube[:,int(y),int(x)] if
                                                     self.errorcube is not None
                                                     else None),
                               unit=self.unit,
                               model_registry=self.Registry,
                              )

        sp.specfit = self.specfit.copy(parent=sp, registry=sp.Registry)
        # explicitly re-do this (test)
        sp.specfit.includemask = self.specfit.includemask.copy()
        sp.specfit.Spectrum = sp

        if hasattr(self, 'parcube'):
            if self.has_fit[int(y),int(x)]:
                # only set parameters if they're valid
                sp.specfit.modelpars = self.parcube[:,int(y),int(x)]
                if hasattr(self.specfit,'parinfo') and self.specfit.parinfo is not None:
                    # set the parinfo values correctly for annotations
                    for pi,p,e in zip(sp.specfit.parinfo, sp.specfit.modelpars, self.errcube[:,int(y),int(x)]):
                        try:
                            pi['value'] = p
                            pi['error'] = e
                        except ValueError:
                            pass

                if hasattr(self.specfit,'fitter') and self.specfit.fitter is not None:
                    sp.specfit.fitter.mpp = sp.specfit.modelpars # also for annotations (differs depending on which function... sigh... need to unify)
                    sp.specfit.npeaks = self.specfit.fitter.npeaks
                    sp.specfit.fitter.npeaks = len(sp.specfit.modelpars) / sp.specfit.fitter.npars
                    sp.specfit.fitter.parinfo = sp.specfit.parinfo
                    try:
                        sp.specfit.model = sp.specfit.fitter.n_modelfunc(sp.specfit.modelpars,
                                                                         **sp.specfit.fitter.modelfunc_kwargs)(sp.xarr)
                    except ValueError:
                        # possibly invalid model parameters, just skip
                        sp.specfit.model = np.zeros_like(sp.data)

        return sp

    def get_apspec(self, aperture, coordsys=None, method='mean', **kwargs):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube pixel coordinates)

        *aperture* [tuple or list] (x, y, radius)
            The aperture to use when extracting the data

        *coordsys* [ 'celestial' | 'galactic' | None]
            the coordinate system the aperture is specified in
            None indicates pixel coordinates (default)

        *wunit* [str]
            arcsec, arcmin, or degree


        """

        if coordsys is not None:
            wcs = self.mapplot.wcs
        else:
            wcs = None

        data = cubes.extract_aperture(self.cube, aperture,
                                      coordsys=coordsys,
                                      wcs=wcs,
                                      method=method,
                                      **kwargs)
        if self.errorcube is not None:
            error = cubes.extract_aperture(self.errorcube, aperture,
                                           coordsys=coordsys,
                                           wcs=self.mapplot.wcs,
                                           method='error', **kwargs)
        else:
            error = None

        ct = 'CTYPE{0}'.format(self._first_cel_axis_num)
        header = cubes.speccen_header(fits.Header(cards=[(k,v) for k,v in
                                                         iteritems(self.header)
                                                         if k != 'HISTORY']),
                                      lon=aperture[0],
                                      lat=aperture[1],
                                      system=self.system,
                                      proj=self.header[ct][-3:])
        if len(aperture) == 3:
            header['APRADIUS'] = aperture[2]
        if len(aperture) == 5:
            header['APMAJ'] = aperture[2]
            header['APMIN'] = aperture[3]
            header['APREFF'] = (aperture[2]*aperture[3])**0.5
            header['APPA'] = aperture[4]

        sp = spectrum.Spectrum(xarr=self.xarr.copy(),
                               data=data,
                               error=error,
                               header=header,
                               model_registry=self.Registry,
                              )

        sp.specfit = self.specfit.copy(parent=sp, registry=sp.Registry)

        return sp

    def set_apspec(self, aperture, coordsys=None, method='mean'):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube coordinates)
        """

        if coordsys is not None:
            self.data = cubes.extract_aperture(self.cube, aperture,
                                               coordsys=coordsys,
                                               wcs=self.mapplot.wcs,
                                               method=method)
        else:
            self.data = cubes.extract_aperture(self.cube, aperture,
                                               coordsys=None, method=method)

    def get_modelcube(self, update=False, multicore=1):
        """
        Return or generate a "model cube", which will have the same shape as
        the ``.cube`` but will have spectra generated from the fitted model.

        If the model cube does not yet exist, one will be generated

        Parameters
        ----------
        update : bool
            If the cube has already been computed, set this to ``True`` to
            recompute the model.
        multicore: int
            if >1, try to use multiprocessing via parallel_map to run on multiple cores
        """
        if self._modelcube is None or update:
            yy,xx = np.indices(self.parcube.shape[1:])
            nanvals = np.any(~np.isfinite(self.parcube),axis=0)
            isvalid = np.any(self.parcube, axis=0) & ~nanvals
            valid_pixels = zip(xx[isvalid], yy[isvalid])
            self._modelcube = np.full_like(self.cube, np.nan)

            def model_a_pixel(xy):
                x,y = int(xy[0]), int(xy[1])
                self._modelcube[:,y,x] = self.specfit.get_full_model(pars=self.parcube[:,y,x])
                return ((x,y), self._modelcube[:,y,x])

            if multicore > 1:
                sequence = [(x,y) for x,y in valid_pixels]
                result = parallel_map(model_a_pixel, sequence, numcores=multicore)
                merged_result = [core_result for core_result in result
                                 if core_result is not None]
                for mr in merged_result:
                    ((x,y), model) = mr
                    x = int(x)
                    y = int(y)
                    self._modelcube[:,y,x] = model
            else:
                # progressbar doesn't work with zip; I'm therefore giving up on
                # "efficiency" in memory by making a list here.
                for xy in ProgressBar(list(valid_pixels)):
                    model_a_pixel(xy)

        return self._modelcube


    def fiteach(self, errspec=None, errmap=None, guesses=(), verbose=True,
                verbose_level=1, quiet=True, signal_cut=3, usemomentcube=None,
                blank_value=0, integral=False, direct_integral=False,
                absorption=False, use_nearest_as_guess=False,
                use_neighbor_as_guess=False, start_from_point=(0,0),
                multicore=1, position_order=None, continuum_map=None,
                prevalidate_guesses=False, maskmap=None,
                skip_failed_fits=False,
                **fitkwargs):
        """
        Fit a spectrum to each valid pixel in the cube

        For guesses, priority is *use_nearest_as_guess*, *usemomentcube*,
        *guesses*, None

        Once you have successfully run this function, the results will be
        stored in the ``.parcube`` and ``.errcube`` attributes, which are each
        cubes of shape ``[npars, ny, nx]``, where npars is the number of fitted
        parameters and ``nx``, ``ny`` are the shape of the map.  ``errcube``
        contains the errors on the fitted parameters (1-sigma, as returned from
        the Levenberg-Marquardt fit's covariance matrix).  You can use the
        attribute ``has_fit``, which is a map of shape ``[ny,nx]`` to find
        which pixels have been successfully fit.


        Parameters
        ----------
        use_nearest_as_guess: bool
            Unless the fitted point is the first, it will find the nearest
            other point with a successful fit and use its best-fit parameters
            as the guess
        use_neighbor_as_guess: bool
            Set this keyword to use the average best-fit parameters from
            neighboring positions with successful fits as the guess
        start_from_point: tuple(int,int)
            Either start from the center or from a point defined by a tuple.
            Work outward from that starting point.
        position_order: ndarray[naxis=2]
            2D map of region with pixel values indicating the order in which
            to carry out the fitting.  Any type with increasing pixel values.
        guesses: tuple or ndarray[naxis=3]
            Either a tuple/list of guesses with len(guesses) = npars or a cube
            of guesses with shape [npars, ny, nx].
            NOT TRUE, but a good idea in principle:
            You can also use a dictionary of the form {(y,x): [list of length
            npars]}, where (y,x) specifies a pixel location. If the dictionary
            method is used, npars must be specified and it sets the length of
            the first parameter axis
        signal_cut: float
            Minimum signal-to-noise ratio to "cut" on (i.e., if peak in a given
            spectrum has s/n less than this value, ignore it)
        blank_value: float
            Value to replace non-fitted locations with.
        errmap: ndarray[naxis=2] or ndarray[naxis=3]
            A map of errors used for the individual pixels of the spectral
            cube. 2D errmap results in an equal weighting of each given
            spectrum, while a 3D array sets individual weights of each channel
        verbose: bool
        verbose_level: int
            Controls how much is output.
            0,1 - only changes frequency of updates in loop
            2 - print out messages when skipping pixels
            3 - print out messages when fitting pixels
            4 - specfit will be verbose
        multicore: int
            if >1, try to use multiprocessing via parallel_map to run on
            multiple cores
        continuum_map: np.ndarray
            Same shape as error map.  Subtract this from data before estimating
            noise.
        prevalidate_guesses: bool
            An extra check before fitting is run to make sure the guesses are
            all within the specified limits.  May be slow, so it is off by
            default.  It also should not be necessary, since careful checking
            is performed before each fit.
        maskmap : `np.ndarray`, optional
            A boolean mask map, where ``True`` implies that the data are good.
            This will be used for both plotting using `mapplot` and fitting
            using `fiteach`.  If ``None``, will use ``self.maskmap``.
        integral : bool
            If set, the integral of each spectral fit will be computed and
            stored in the attribute ``.integralmap``
        direct_integral : bool
            Return the integral of the *spectrum* (as opposed to the fitted
            model) over a range defined by the `integration_limits` if specified or
            `threshold` otherwise
        skip_failed_fits : bool
            Flag to forcibly skip failed fits that fail with "unknown error".
            Generally, you do not want this on, but this is the
            'finger-over-the-engine-light' approach that will allow these
            incomprehensible failures to go by and just ignore them.  Keep
            an eye on how many of these you get: if it's just one or two
            out of hundreds, then maybe those are just pathological cases
            that can be ignored.  If it's a significant fraction, you probably
            want to take a different approach.
        """
        if 'multifit' in fitkwargs:
            warn("The multifit keyword is no longer required.  All fits "
                 "allow for multiple components.", DeprecationWarning)

        if not hasattr(self.mapplot,'plane'):
            self.mapplot.makeplane()

        if maskmap is None:
            maskmap = self.maskmap

        yy,xx = np.indices(self.mapplot.plane.shape)
        if isinstance(self.mapplot.plane, np.ma.core.MaskedArray):
            OK = ((~self.mapplot.plane.mask) &
                  maskmap.astype('bool')).astype('bool')
        else:
            OK = (np.isfinite(self.mapplot.plane) &
                  maskmap.astype('bool')).astype('bool')

        # NAN guesses rule out the model too
        if hasattr(guesses,'shape') and guesses.shape[1:] == self.cube.shape[1:]:
            bad = np.isnan(guesses).sum(axis=0).astype('bool')
            OK &= (~bad)

        log.info("Fitting up to {0} spectra".format(OK.sum()))

        if start_from_point == 'center':
            start_from_point = (xx.max()/2., yy.max()/2.)
        if hasattr(position_order,'shape') and position_order.shape == self.cube.shape[1:]:
            sort_distance = np.argsort(position_order.flat)
        else:
            d_from_start = ((xx-start_from_point[1])**2 + (yy-start_from_point[0])**2)**0.5
            sort_distance = np.argsort(d_from_start.flat)

        if use_neighbor_as_guess or use_nearest_as_guess:
            distance = ((xx)**2 + (yy)**2)**0.5

        valid_pixels = list(zip(xx.flat[sort_distance][OK.flat[sort_distance]],
                                yy.flat[sort_distance][OK.flat[sort_distance]]))

        if len(valid_pixels) != len(set(valid_pixels)):
            raise ValueError("There are non-unique pixels in the 'valid pixel' list.  "
                             "This should not be possible and indicates a major error.")
        elif len(valid_pixels) == 0:
            raise ValueError("No valid pixels selected.")

        if start_from_point not in valid_pixels:
            raise ValueError("The starting fit position is not among the valid"
                             " pixels.  Check your selection criteria to make "
                             "sure you have not unintentionally excluded "
                             "this first fit pixel.")

        if verbose_level > 0:
            log.debug("Number of valid pixels: %i" % len(valid_pixels))

        guesses_are_moments = (isinstance(guesses, string_types) and
                                 guesses in ('moment','moments'))
        if guesses_are_moments or (usemomentcube and len(guesses)):
            if not hasattr(self, 'momentcube') and guesses_are_moments:
                self.momenteach()
            npars = self.momentcube.shape[0]
        else:
            npars = len(guesses)
            if npars == 0:
                raise ValueError("Parameter guesses are required.")

        self.parcube = np.zeros((npars,)+self.mapplot.plane.shape)
        self.errcube = np.zeros((npars,)+self.mapplot.plane.shape)
        if integral:
            self.integralmap = np.zeros((2,)+self.mapplot.plane.shape)

        # newly needed as of March 27, 2012.  Don't know why.
        if 'fittype' in fitkwargs:
            self.specfit.fittype = fitkwargs['fittype']
        self.specfit.fitter = self.specfit.Registry.multifitters[self.specfit.fittype]

        # TODO: VALIDATE THAT ALL GUESSES ARE WITHIN RANGE GIVEN THE
        # FITKWARG LIMITS

        # array to store whether pixels have fits
        self.has_fit = np.zeros(self.mapplot.plane.shape, dtype='bool')

        self._counter = 0

        self._tracebacks = {}

        t0 = time.time()

        def fit_a_pixel(iixy):
            ii,x,y = iixy
            sp = self.get_spectrum(x,y)

            # very annoying - cannot use min/max without checking type
            # maybe can use np.asarray here?
            # cannot use sp.data.mask because it can be a scalar boolean,
            # which does unpredictable things.
            if hasattr(sp.data, 'mask') and not isinstance(sp.data.mask, (bool,
                                                                          np.bool_)):
                sp.data[sp.data.mask] = np.nan
                sp.error[sp.data.mask] = np.nan
                sp.data = np.array(sp.data)
                sp.error = np.array(sp.error)

            if errspec is not None:
                sp.error = errspec
            elif errmap is not None:
                if self.errorcube is not None:
                    raise ValueError("Either the 'errmap' argument or"
                                     " self.errorcube attribute should be"
                                     " specified, but not both.")
                if errmap.shape == self.cube.shape[1:]:
                    sp.error = np.ones(sp.data.shape) * errmap[int(y),int(x)]
                elif errmap.shape == self.cube.shape:
                    sp.error = errmap[:, int(y), int(x)]
            elif self.errorcube is not None:
                sp.error = self.errorcube[:, int(y), int(x)]
            else:
                if ii==0:
                    # issue the warning only once (ii==0), but always issue
                    warn("Using data std() as error.  "
                         "If signal_cut is set, this can result in "
                         "some pixels not being fit.",
                         PyspeckitWarning)
                sp.error[:] = sp.data[sp.data==sp.data].std()
            if sp.error is None:
                raise TypeError("The Spectrum's error is unset.  This should "
                                "not be possible.  Please raise an Issue.")
            if signal_cut > 0 and not all(sp.error == 0):
                if continuum_map is not None:
                    with np.errstate(divide='raise'):
                        snr = (sp.data-continuum_map[int(y),int(x)]) / sp.error
                else:
                    with np.errstate(divide='raise'):
                        snr = sp.data / sp.error
                if absorption:
                    max_sn = np.nanmax(-1*snr)
                else:
                    max_sn = np.nanmax(snr)
                if max_sn < signal_cut:
                    if verbose_level > 1:
                        log.info("Skipped %4i,%4i (s/n=%0.2g)" % (x,y,max_sn))
                    return
                elif np.isnan(max_sn):
                    if verbose_level > 1:
                        log.info("Skipped %4i,%4i (s/n is nan; max(data)=%0.2g, min(error)=%0.2g)" %
                                 (x,y,np.nanmax(sp.data),np.nanmin(sp.error)))
                    return
                if verbose_level > 2:
                    log.info("Fitting %4i,%4i (s/n=%0.2g)" % (x,y,max_sn))
            else:
                max_sn = None
            sp.specfit.Registry = self.Registry # copy over fitter registry

            # Do some homework for local fits
            # Exclude out of bounds points
            xpatch, ypatch = get_neighbors(x,y,self.has_fit.shape)
            local_fits = self.has_fit[ypatch+y,xpatch+x]

            if use_nearest_as_guess and self.has_fit.sum() > 0:
                if verbose_level > 1 and ii == 0 or verbose_level > 4:
                    log.info("Using nearest fit as guess")
                rolled_distance = np.roll(np.roll(distance, x, 0), y, 1)
                # If there's no fit, set its distance to be unreasonably large
                # so it will be ignored by argmin
                nearest_ind = np.argmin(rolled_distance+1e10*(~self.has_fit))
                nearest_x, nearest_y = xx.flat[nearest_ind],yy.flat[nearest_ind]
                if np.all(np.isfinite(self.parcube[:,nearest_y,nearest_x])):
                    gg = self.parcube[:,nearest_y,nearest_x]
                else:
                    log.exception("Pixel {0},{1} had a fit including a NaN: {2}"
                                  " so it will not be used as a guess for {3},{4}"
                                  .format(nearest_x, nearest_y, self.parcube[:, nearest_y, nearest_x],
                                          x, y))
                    gg = guesses
            elif use_neighbor_as_guess and np.any(local_fits):
                # Array is N_guess X Nvalid_nbrs so averaging over
                # Axis=1 is the axis of all valid neighbors
                gg = np.mean(self.parcube[:,
                                          (ypatch+y)[local_fits],
                                          (xpatch+x)[local_fits]], axis=1)
                if np.any(~np.isfinite(gg)):
                    log.exception("Pixel {0},{1} neighbors had non-finite guess: {2}"
                                  .format(x, y, gg))
                    gg = guesses
            elif guesses_are_moments and usemomentcube is False:
                raise ValueError("usemomentcube must be set to True")
            elif guesses_are_moments or (usemomentcube and len(guesses)):
                if not guesses_are_moments and ii == 0:
                    log.warn("guesses will be ignored because usemomentcube "
                             "was set to True.", PyspeckitWarning)
                if verbose_level > 1 and ii == 0:
                    log.info("Using moment cube")
                gg = self.momentcube[:,int(y),int(x)]
            elif hasattr(guesses,'shape') and guesses.shape[1:] == self.cube.shape[1:]:
                if verbose_level > 1 and ii == 0:
                    log.info("Using input guess cube")
                gg = guesses[:,int(y),int(x)]
            elif isinstance(guesses, dict):
                if verbose_level > 1 and ii == 0:
                    log.info("Using input guess dict")
                gg = guesses[(int(y),int(x))]
            else:
                if verbose_level > 1 and ii == 0:
                    log.info("Using input guess")
                gg = guesses

            if np.all(np.isfinite(gg)):
                try:
                    with np.errstate(divide='raise'):
                        sp.specfit(guesses=gg, quiet=verbose_level<=3,
                                   verbose=verbose_level>3, **fitkwargs)
                    self.parcube[:,int(y),int(x)] = sp.specfit.modelpars
                    self.errcube[:,int(y),int(x)] = sp.specfit.modelerrs
                    if np.any(~np.isfinite(sp.specfit.modelpars)):
                        log.exception("Fit result included nan for pixel {0},{1}: "
                                      "{2}".format(x, y, sp.specfit.modelpars))
                        success = False
                        # this is basically a debug statement to try to get the
                        # code to crash here
                        raise KeyboardInterrupt
                    else:
                        success = True
                except Exception as ex:
                    exc_traceback = sys.exc_info()[2]
                    self._tracebacks[(ii,x,y)] = exc_traceback
                    log.exception("Fit number %i at %i,%i failed on error %s" % (ii,x,y, str(ex)))
                    log.exception("Failure was in file {0} at line {1}".format(
                        exc_traceback.tb_frame.f_code.co_filename,
                        exc_traceback.tb_lineno,))
                    traceback.print_tb(exc_traceback)
                    log.exception("Guesses were: {0}".format(str(gg)))
                    log.exception("Fitkwargs were: {0}".format(str(fitkwargs)))
                    success = False
                    if isinstance(ex, KeyboardInterrupt):
                        raise ex

                # keep this out of the 'try' statement
                if integral and success:
                    self.integralmap[:,int(y),int(x)] = sp.specfit.integral(direct=direct_integral,
                                                                            return_error=True)
                self.has_fit[int(y),int(x)] = success
            else:
                log.exception("Fit number {0} at {1},{2} had non-finite guesses {3}"
                              .format(ii, x, y, guesses))
                self.has_fit[int(y),int(x)] = False
                self.parcube[:,int(y),int(x)] = blank_value
                self.errcube[:,int(y),int(x)] = blank_value
                if integral:
                    self.integralmap[:,int(y),int(x)] = blank_value


            self._counter += 1
            if verbose:
                if ii % (min(10**(3-verbose_level),1)) == 0:
                    snmsg = " s/n=%5.1f" % (max_sn) if max_sn is not None else ""
                    npix = len(valid_pixels)
                    pct = 100 * (ii+1.0)/float(npix)
                    log.info("Finished fit %6i of %6i at (%4i,%4i)%s. Elapsed time is %0.1f seconds.  %%%01.f" %
                             (ii+1, npix, x, y, snmsg, time.time()-t0, pct))

            if sp.specfit.modelerrs is None:
                log.exception("Fit number %i at %i,%i failed with no specific error." % (ii,x,y))
                if hasattr(sp.specfit, 'mpfit_status'):
                    log.exception("mpfit status is {0}".format(sp.specfit.mpfit_status))
                log.exception("The problem is that the model errors were never set, "
                              "which implies that the fit simply failed to finish.")
                log.exception("The string representation of `sp.specfit.parinfo` is: {0}"
                              .format(sp.specfit.parinfo))
                log.exception("The string representation of `sp.specfit.fitter.parinfo` is: {0}"
                              .format(sp.specfit.fitter.parinfo))
                log.exception("modelpars is: {0}".format(str(sp.specfit.modelpars)))
                log.exception("cube modelpars are: {0}".format(str(self.parcube[:,int(y),int(x)])))
                log.exception("cube modelerrs are: {0}".format(str(self.errcube[:,int(y),int(x)])))
                log.exception("Guesses were: {0}".format(str(gg)))
                log.exception("Fitkwargs were: {0}".format(str(fitkwargs)))

                if skip_failed_fits:
                    # turn the flag into a count
                    log.exception("The fit never completed; something has gone wrong.  Failed fits = {0}".format(int(skip_failed_fits)))
                else:
                    raise TypeError("The fit never completed; something has gone wrong.")


            # blank out the errors (and possibly the values) wherever they are zero = assumed bad
            # this is done after the above exception to make sure we can inspect these values
            if blank_value != 0:
                self.parcube[self.parcube == 0] = blank_value
                self.errcube[self.parcube == 0] = blank_value

            if integral:
                return ((x,y), sp.specfit.modelpars, sp.specfit.modelerrs,
                        self.integralmap[:,int(y),int(x)])
            else:
                return ((x,y), sp.specfit.modelpars, sp.specfit.modelerrs)
        #### BEGIN TEST BLOCK ####
        # This test block is to make sure you don't run a 30 hour fitting
        # session that's just going to crash at the end.
        # try a first fit for exception-catching
        if len(start_from_point) == 2:
            try0 = fit_a_pixel((0,start_from_point[0], start_from_point[1]))
        else:
            try0 = fit_a_pixel((0,valid_pixels[0][0],valid_pixels[0][1]))
        try:
            len_guesses = len(self.momentcube) if (usemomentcube or
                                guesses_are_moments) else len(guesses)
            assert len(try0[1]) == len_guesses == len(self.parcube) == len(self.errcube)
            assert len(try0[2]) == len_guesses == len(self.parcube) == len(self.errcube)
        except TypeError as ex:
            if try0 is None:
                raise AssertionError("The first fitted pixel did not yield a "
                                     "fit. Please try starting from a "
                                     "different pixel.")
            else:
                raise ex
        except AssertionError:
            raise AssertionError("The first pixel had the wrong fit "
                                 "parameter shape.  This is probably "
                                 "a bug; please report it.")

        # This is a secondary test... I'm not sure it's necessary, but it
        # replicates what's inside the fit_a_pixel code and so should be a
        # useful sanity check
        x,y = valid_pixels[0]
        sp = self.get_spectrum(x,y)
        sp.specfit.Registry = self.Registry # copy over fitter registry
        # this reproduced code is needed because the functional wrapping
        # required for the multicore case prevents gg from being set earlier
        if usemomentcube or guesses_are_moments:
            gg = self.momentcube[:,int(y),int(x)]
        elif hasattr(guesses,'shape') and guesses.shape[1:] == self.cube.shape[1:]:
            gg = guesses[:,int(y),int(x)]
        else:
            gg = guesses

        # This is NOT in a try/except block because we want to raise the
        # exception here if an exception is going to happen
        sp.specfit(guesses=gg, **fitkwargs)

        if prevalidate_guesses:
            if guesses.ndim == 3:
                for ii,(x,y) in ProgressBar(tuple(enumerate(valid_pixels))):
                    pinf, _ = sp.specfit.fitter._make_parinfo(parvalues=guesses[:,int(y),int(x)], **fitkwargs)
                    sp.specfit._validate_parinfo(pinf, 'raise')
            else:
                pinf, _ = sp.specfit.fitter._make_parinfo(parvalues=guesses, **fitkwargs)
                sp.specfit._validate_parinfo(pinf, 'raise')

        #### END TEST BLOCK ####


        if multicore > 1:
            sequence = [(ii,x,y) for ii,(x,y) in tuple(enumerate(valid_pixels))]
            with np.errstate(divide='raise'):
                result = parallel_map(fit_a_pixel, sequence, numcores=multicore)
            self._result = result # backup - don't want to lose data in the case of a failure
            # a lot of ugly hacking to deal with the way parallel_map returns
            # its results needs TWO levels of None-filtering, because any
            # individual result can be None (I guess?) but apparently (and this
            # part I don't believe) any individual *fit* result can be None as
            # well (apparently the x,y pairs can also be None?)
            merged_result = [core_result for core_result in result if
                             core_result is not None]
            # for some reason, every other time I run this code, merged_result
            # ends up with a different intrinsic shape.  This is an attempt to
            # force it to maintain a sensible shape.
            try:
                if integral:
                    ((x,y), m1, m2, intgl) = merged_result[0]
                else:
                    ((x,y), m1, m2) = merged_result[0]
            except ValueError:
                if verbose > 1:
                    log.exception("ERROR: merged_result[0] is {0} which has the"
                                  " wrong shape".format(merged_result[0]))
                merged_result = itertools.chain.from_iterable(merged_result)
            for TEMP in merged_result:
                if TEMP is None:
                    # this shouldn't be possible, but it appears to happen
                    # anyway.  parallel_map is great, up to a limit that was
                    # reached long before this level of complexity
                    log.debug("Skipped a None entry: {0}".format(str(TEMP)))
                    continue
                try:
                    if integral:
                        ((x,y), modelpars, modelerrs, intgl) = TEMP
                    else:
                        ((x,y), modelpars, modelerrs) = TEMP
                except TypeError:
                    # implies that TEMP does not have the shape ((a,b),c,d)
                    # as above, shouldn't be possible, but it happens...
                    log.debug("Skipped a misshapen entry: {0}".format(str(TEMP)))
                    continue
                if ((len(modelpars) != len(modelerrs)) or
                    (len(modelpars) != len(self.parcube))):
                    raise ValueError("There was a serious problem; modelpar and"
                                     " error shape don't match that of the "
                                     "parameter cubes")
                if ((any([x is None for x in modelpars]) or
                     np.any(np.isnan(modelpars)) or
                     any([x is None for x in modelerrs]) or
                     np.any(np.isnan(modelerrs)))):
                    self.parcube[:,int(y),int(x)] = np.nan
                    self.errcube[:,int(y),int(x)] = np.nan
                    self.has_fit[int(y),int(x)] = False
                else:
                    self.parcube[:,int(y),int(x)] = modelpars
                    self.errcube[:,int(y),int(x)] = modelerrs
                    self.has_fit[int(y),int(x)] = max(modelpars) > 0
                if integral:
                    self.integralmap[:,int(y),int(x)] = intgl
        else:
            for ii,(x,y) in enumerate(valid_pixels):
                fit_a_pixel((ii,x,y))


        # March 27, 2014: This is EXTREMELY confusing.  This isn't in a loop...
        # make sure the fitter / fittype are set for the cube
        # this has to be done within the loop because skipped-over spectra
        # don't ever get their fittypes set
        self.specfit.fitter = sp.specfit.fitter
        self.specfit.fittype = sp.specfit.fittype
        self.specfit.parinfo = sp.specfit.parinfo

        if verbose:
            log.info("Finished final fit %i.  "
                     "Elapsed time was %0.1f seconds" % (len(valid_pixels), time.time()-t0))

        pars_are_finite = np.all(np.isfinite(self.parcube), axis=0)

        # if you see one of these exceptions, please try to produce a minimum
        # working example and report it as a bug.
        # all non-finite fit parameters should be has_fit=False
        assert np.all(~self.has_fit[~pars_are_finite]), "Non-finite parameters found in fits"


    def momenteach(self, verbose=True, verbose_level=1, multicore=1, **kwargs):
        """
        Return a cube of the moments of each pixel

        Parameters
        ----------
        multicore: int
            if >1, try to use multiprocessing via parallel_map to run on multiple cores
        """

        if not hasattr(self.mapplot,'plane'):
            self.mapplot.makeplane()

        if 'vheight' not in kwargs:
            kwargs['vheight'] = False

        yy,xx = np.indices(self.mapplot.plane.shape)
        if isinstance(self.mapplot.plane, np.ma.core.MaskedArray):
            OK = (~self.mapplot.plane.mask) * self.maskmap
        else:
            OK = np.isfinite(self.mapplot.plane) * self.maskmap
        valid_pixels = zip(xx[OK],yy[OK])

        # run the moment process to find out how many elements are in a moment
        _temp_moment = self.get_spectrum(xx[OK][0], yy[OK][0]).moments(**kwargs)

        self.momentcube = np.zeros((len(_temp_moment),)+self.mapplot.plane.shape)

        t0 = time.time()

        def moment_a_pixel(iixy):
            ii,x,y = iixy
            sp = self.get_spectrum(x,y)
            self.momentcube[:,int(y),int(x)] = sp.moments(**kwargs)
            if verbose:
                if ii % 10**(3-verbose_level) == 0:
                    log.info("Finished moment %i.  "
                             "Elapsed time is %0.1f seconds" % (ii, time.time()-t0))

            return ((x,y), self.momentcube[:,int(y),int(x)])

        if multicore > 1:
            sequence = [(ii,x,y) for ii,(x,y) in tuple(enumerate(valid_pixels))]
            result = parallel_map(moment_a_pixel, sequence, numcores=multicore)
            merged_result = [core_result.tolist()
                             for core_result in result
                             if core_result is not None]
            for TEMP in merged_result:
                ((x,y), moments) = TEMP
                self.momentcube[:,int(y),int(x)] = moments
        else:
            for ii,(x,y) in enumerate(valid_pixels):
                moment_a_pixel((ii,x,y))

        if verbose:
            log.info("Finished final moment %i.  "
                     "Elapsed time was %0.1f seconds" % (OK.sum(), time.time()-t0))

    def show_moment(self, momentnumber, **kwargs):
        """
        If moments have been computed, display them in the mapplot window
        """

        if not hasattr(self,'momentcube'):
            raise ValueError("Compute moments first")

        self.mapplot.plane = self.momentcube[momentnumber,:,:].squeeze()

        self.mapplot(estimator=None, **kwargs)

    def show_fit_param(self, parnumber, **kwargs):
        """
        If pars have been computed, display them in the mapplot window

        Parameters
        ----------
        parnumber : int
            The index of the parameter in the parameter cube
        """

        if not hasattr(self,'parcube'):
            raise ValueError("Compute fit parameters first")

        self.mapplot.plane = self.parcube[parnumber,:,:].squeeze()

        self.mapplot(estimator=None, **kwargs)


    def load_model_fit(self, fitsfilename, npars, npeaks=1, fittype=None,
                       _temp_fit_loc=(0,0)):
        """
        Load a parameter + error cube into the ``.parcube`` and ``.errcube``
        attributes.  The models can then be examined and plotted using
        ``.mapplot`` as if you had run ``.fiteach``.

        Parameters
        ----------
        fitsfilename : str
            The filename containing the parameter cube written with `write_fit`
        npars : int
            The number of parameters in the model fit for a single spectrum
        npeaks : int
            The number of independent peaks fit toward each spectrum
        fittype : str, optional
            The name of the fittype, e.g. 'gaussian' or 'voigt', from the
            pyspeckit fitter registry.  This is optional; it should have
            been written to the FITS header and will be read from there if
            it is not specified
        _temp_fit_loc : tuple (int,int)
            The initial spectrum to use to generate components of the class.
            This should not need to be changed.
        """
        try:
            import astropy.io.fits as pyfits
        except ImportError:
            import pyfits

        cubefile = pyfits.open(fitsfilename,ignore_missing_end=True)
        cube = cubefile[0].data

        if cube.shape[0] != npars * npeaks * 2:
            raise ValueError("The cube shape is not correct.  The cube has "
                             "first dimension = {0}, but it should be {1}. "
                             "The keyword npars = number of parameters per "
                             "model component, and npeaks = number of "
                             "independent peaks.  You gave npars={2} and "
                             "npeaks={3}".format(cube.shape[0], npars*npeaks*2,
                                                 npars, npeaks))

        # grab a spectrum and fit it however badly you want
        # this is just to __init__ the relevant data structures
        if fittype is None:
            if cubefile[0].header.get('FITTYPE'):
                fittype = cubefile[0].header.get('FITTYPE')
            else:
                raise KeyError("Must specify FITTYPE or include it in cube header.")

        self.parcube = cube[:npars*npeaks,:,:]
        self.errcube = cube[npars*npeaks:npars*npeaks*2,:,:]

        if np.any(np.all(self.parcube == 0, axis=(1,2))):
            # there are some slices where all parameters are zero, we should
            # ignore this when establishing whether there's a fit (some
            # parameters, like fortho, can be locked to zero)
            self.has_fit = np.all((np.isfinite(self.parcube)), axis=0)
        else:
            self.has_fit = np.all((self.parcube != 0) &
                                  (np.isfinite(self.parcube)), axis=0)

        nanvals = ~np.isfinite(self.parcube)
        nanvals_flat = np.any(nanvals, axis=0)
        if np.any(nanvals):
            warn("NaN or infinite values encountered in parameter cube.",
                 PyspeckitWarning)


        # make sure params are within limits
        fitter = self.specfit.Registry.multifitters[fittype]
        guesses,throwaway = fitter._make_parinfo(npeaks=npeaks)

        try:
            x,y = _temp_fit_loc
            sp = self.get_spectrum(x,y)
            guesses.values = self.parcube[:,int(y),int(x)]
            sp.specfit(fittype=fittype, guesses=guesses.values)
        except Exception as ex1:
            try:
                OKmask = np.any(self.parcube, axis=0) & ~nanvals_flat
                whereOK = np.where(OKmask)
                x,y = whereOK[1][0],whereOK[0][0]
                sp = self.get_spectrum(x,y)
                guesses.values = self.parcube[:,int(y),int(x)]
                sp.specfit(fittype=fittype, guesses=guesses.values)
            except Exception as ex2:
                log.error("Fitting the pixel at location {0} failed with error: {1}.  "
                          "Re-trying at location {2} failed with error {3}.  "
                          "Try setting _temp_fit_loc to a valid pixel".format(_temp_fit_loc, ex1,
                                                                              (x,y), ex2))

        self.specfit.fitter = sp.specfit.fitter
        self.specfit.fittype = sp.specfit.fittype
        self.specfit.parinfo = sp.specfit.parinfo

    def smooth(self,factor,**kwargs):
        """
        Smooth the spectrum by factor `factor`.

        Documentation from the :mod:`cubes.spectral_smooth` module:

        """
        factor = round(factor)
        self.cube = cubes.spectral_smooth(self.cube,factor,**kwargs)
        self.xarr = self.xarr[::factor]
        if hasattr(self,'data'):
            self.data = smooth.smooth(self.data,factor,**kwargs)
        if len(self.xarr) != self.cube.shape[0]:
            raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
        if self.errorcube is not None:
            self.errorcube = cubes.spectral_smooth(self.errorcube,factor,**kwargs)

        self._smooth_header(factor)

    __doc__ += "cubes.spectral_smooth doc: \n" + cubes.spectral_smooth.__doc__

    def _smooth_header(self,factor):
        """
        Internal - correct the FITS header parameters when smoothing
        """
        if self.header.get('CDELT3') is not None and self.header.get('CRPIX3') is not None:
            self.header['CDELT3'] = self.header.get('CDELT3') * float(factor)
            self.header['CRPIX3'] = self.header.get('CRPIX3') / float(factor)

            history.write_history(self.header,"SMOOTH: Smoothed and downsampled spectrum by factor %i" % (factor))
            history.write_history(self.header,"SMOOTH: Changed CRPIX3 from %f to %f" % (self.header.get('CRPIX3')*float(factor),self.header.get('CRPIX3')))
            history.write_history(self.header,"SMOOTH: Changed CDELT3 from %f to %f" % (self.header.get('CRPIX3')/float(factor),self.header.get('CRPIX3')))


    def write_fit(self, fitcubefilename, overwrite=False):
        """
        Write out a fit cube containing the ``.parcube`` and ``.errcube`` using
        the information in the fit's parinfo to set the header keywords.  The
        ``PLANE#`` keywords will be used to indicate the content of each plane
        in the data cube written to the FITS file.  All of the fitted
        parameters will be written first, followed by all of the errors on
        those parameters.  So, for example, if you have fitted a single
        gaussian to each pixel, the dimensions of the saved cube will be ``[6,
        ny, nx]``, and they will be the amplitude, centroid, width, error on
        amplitude, error on centroid, and error on width, respectively.

        To load such a file back in for plotting purposes, see
        `SpectralCube.load_model_fit`.

        Parameters
        ----------
        fitcubefilename: string
            Filename to write to
        overwrite: bool
            Overwrite file if it exists?
        """

        try:
            import astropy.io.fits as pyfits
        except ImportError:
            import pyfits

        try:
            fitcubefile = pyfits.PrimaryHDU(data=np.concatenate([self.parcube,self.errcube]), header=self.header)
            fitcubefile.header['FITTYPE'] = self.specfit.fittype

            for ii,par in enumerate(self.specfit.parinfo):
                kw = "PLANE%i" % ii
                parname = par['parname'].strip('0123456789')
                fitcubefile.header[kw] = parname
            # set error parameters
            for jj,par in enumerate(self.specfit.parinfo):
                kw = "PLANE%i" % (ii+jj+1)
                parname = "e"+par['parname'].strip('0123456789')
                fitcubefile.header[kw] = parname

            # overwrite the WCS
            fitcubefile.header['CDELT3'] = 1
            fitcubefile.header['CTYPE3'] = 'FITPAR'
            fitcubefile.header['CRVAL3'] = 0
            fitcubefile.header['CRPIX3'] = 1
        except AttributeError:
            log.exception("Make sure you run the cube fitter first.")
            return

        if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
            fitcubefile.writeto(fitcubefilename, overwrite=overwrite)
        else:
            fitcubefile.writeto(fitcubefilename, clobber=overwrite)

    def write_cube(self):
        raise NotImplementedError


class CubeStack(Cube):
    """
    The Cube equivalent of Spectra: for stitching multiple cubes with the same
    spatial grid but different frequencies together
    """

    def __init__(self, cubelist, xunit='GHz', x0=0, y0=0, maskmap=None, **kwargs):
        """
        Initialize the Cube.  Accepts FITS files.

        x0,y0 - initial spectrum to use (defaults to lower-left corner)
        """

        log.info("Creating Cube Stack")
        cubelist = list(cubelist)
        for ii,cube in enumerate(cubelist):
            if type(cube) is str:
                cube = Cube(cube)
                cubelist[ii] = cube
            if cube.xarr.unit != xunit:
                # convert all inputs to same (non-velocity) unit
                cube.xarr.convert_to_unit(xunit, **kwargs)

        self.cubelist = cubelist

        log.info("Concatenating data")
        self.xarr = SpectroscopicAxes([sp.xarr for sp in cubelist])
        self.cube = np.ma.concatenate([icube.cube for icube in cubelist])

        if np.any([icube.errorcube is not None for icube in cubelist]):
            if all([icube.errorcube is not None for icube in cubelist]):
                self.errorcube = np.ma.concatenate([icube.errorcube for icube in cubelist])
            else:
                raise ValueError("Mismatched error cubes.")
        else:
            self.errorcube = None

        if hasattr(self.cube,'mask'):
            try:
                if self.cube.mask in (False,np.bool_(False)):
                    # mask causes major problems internally for numpy...
                    self.cube = np.array(self.cube)
            except ValueError:
                # this means that self.cube.mask is an array;
                # techically that's alright
                pass
        self._sort()
        self.data = self.cube[:,int(y0),int(x0)]
        self.error = self.errorcube[:,int(y0),int(x0)] if self.errorcube is not None else None

        self.header = cubelist[0].header.copy()
        for cube in cubelist:
            for key,value in cube.header.items():
                if key in ['HISTORY', 'COMMENT']:
                    continue
                self.header[key] = value

        if self.header:
            self.wcs = wcs.WCS(self.header)
            self.wcs.wcs.fix()
            self._spectral_axis_number = self.wcs.wcs.spec+1
            self._first_cel_axis_num = np.where(self.wcs.wcs.axis_types // 1000 == 2)[0][0]+1

            # TODO: Improve this!!!
            self.system = ('galactic'
                           if ('CTYPE{0}'.format(self._first_cel_axis_num)
                               in self.header and 'GLON' in
                               self.header['CTYPE{0}'.format(self._first_cel_axis_num)])
                           else 'celestial')
        else:
            self._spectral_axis_number = 3
            self._first_cel_axis_num = 1
            self.system = 'PIXEL'


        self.unit = cubelist[0].unit
        for cube in cubelist:
            if cube.unit != self.unit:
                raise ValueError("Mismatched units "
                                 "{0} and {1}".format(cube.unit, self.unit))

        self.fileprefix = cubelist[0].fileprefix # first is the best?

        if maskmap is not None:
            self.maskmap = maskmap
        else:
            self.maskmap = np.ones(self.cube.shape[1:],dtype='bool')

        self._register_fitters()
        self.plotter = spectrum.plotters.Plotter(self)
        self.specfit = spectrum.fitters.Specfit(self,Registry=self.Registry)
        self.baseline = spectrum.baseline.Baseline(self)
        self.speclines = spectrum.speclines
        # Initialize writers TO DO: DO WRITERS WORK FOR CUBES?
        self.writer = {}
        for writer in spectrum.writers.writers:
            self.writer[writer] = spectrum.writers.writers[writer](self)

        # Special.  This needs to be modified to be more flexible; for now I need it to work for nh3
        self.plot_special = None
        self.plot_special_kwargs = {}
        self._modelcube = None

        self.mapplot = mapplot.MapPlotter(self)

    def _sort(self):
        """ Sort the data in order of increasing X (could be decreasing, but
        must be monotonic for plotting reasons) """

        indices = self.xarr.argsort()
        self.xarr = self.xarr[indices]
        self.cube = self.cube[indices,:,:]
        if self.errorcube is not None:
            self.errorcube = self.errorcube[indices,:,:]

def get_neighbors(x, y, shape):
    """
    Find the 9 nearest neighbors, excluding self and any out of bounds points
    """
    ysh, xsh = shape
    xpyp = [(ii,jj)
            for ii,jj in itertools.product((-1,0,1),
                                           (-1,0,1))
            if (ii+x < xsh) and (ii+x >= 0)
            and (jj+y < ysh) and (jj+y >= 0)
            and not (ii==0 and jj==0)]
    xpatch, ypatch = zip(*xpyp)

    return np.array(xpatch, dtype='int'), np.array(ypatch, dtype='int')

def test_get_neighbors():
    xp,yp = get_neighbors(0,0,[10,10])
    assert set(xp) == {0,1}
    assert set(yp) == {0,1}

    xp,yp = get_neighbors(0,1,[10,10])
    assert set(xp) == {0,1}
    assert set(yp) == {-1,0,1}

    xp,yp = get_neighbors(5,6,[10,10])
    assert set(xp) == {-1,0,1}
    assert set(yp) == {-1,0,1}

    xp,yp = get_neighbors(9,9,[10,10])
    assert set(xp) == {0,-1}
    assert set(yp) == {0,-1}

    xp,yp = get_neighbors(9,8,[10,10])
    assert set(xp) == {-1,0}
    assert set(yp) == {-1,0,1}
