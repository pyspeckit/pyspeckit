"""
=====
Cubes
=====

Tools to deal with spectroscopic data cubes.  


Many features in Cubes require additional packages:
    smoothing - requires `agpy <https://agpy.googlecode.com/svn/trunk>`_'s smooth and parallel_map routines
    `coords <https://www.stsci.edu/svn/ssb/astrolib/trunk/coords>`_
    `pyregion <git://github.com/leejjoon/pyregion.git>`_
    `pywcs <git://github.com/astropy/astropy.git>`_
    
The 'grunt work' is performed by the :py:mod:`cubes` module

"""
# import parent package
import pyspeckit
from pyspeckit import spectrum
# import local things
import mapplot
import readers
import time
import numpy as np

class Cube(spectrum.Spectrum):

    def __init__(self, filename, x0=0, y0=0, maskfilename=None, maskmap=None, **kwargs):
        """
        Initialize the Cube.  Accepts files in the following formats:
            - .fits

        x0,y0 - initial spectrum to use (defaults to lower-left corner)
        """

        try: 
            self.cube,self.xarr,self.header,self.fitsfile = readers.open_3d_fits(filename, **kwargs)
            self.data = self.cube[:,y0,x0]
            self.error = None
            self.parse_header(self.header)
        except TypeError as inst:
            print "Failed to read fits file: wrong TYPE."
            print inst
            raise inst

        if maskmap is not None:
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
        for writer in spectrum.writers.writers: self.writer[writer] = spectrum.writers.writers[writer](self)
        self.mapplot = mapplot.MapPlotter(self)

    def plot_spectrum(self, x, y, **kwargs):
        """
        Fill the .data array with a real spectrum and plot it
        """

        self.data = self.cube[:,y,x]

        self.plotter(**kwargs)

    def plot_fit(self, x, y, silent=False, **kwargs):
        """
        If fiteach has been run, plot the best fit
        """
        if not hasattr(self,'parcube'):
            if not silent: print "Must run fiteach before plotting a fit.  If you want to fit a single spectrum, use plot_spectrum() and specfit() directly."
            return

        self.data = self.cube[:,y,x]

        self.specfit.modelpars = self.parcube[:,y,x]
        self.specfit.model = self.specfit.fitter.n_modelfunc(self.specfit.modelpars)(self.xarr)

        self.specfit.plot_fit(**kwargs)

    def plot_apspec(self, aperture, coordsys=None, reset_ylimits=True, **kwargs):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube coordinates)
        """

        self.set_apspec(aperture, coordsys=coordsys)
        self.plotter(reset_ylimits=reset_ylimits, **kwargs)

    def get_spectrum(self, x, y):
        """
        Very simple: get the spectrum at coordinates x,y

        Returns a SpectroscopicAxis instance
        """
        return pyspeckit.Spectrum( xarr=self.xarr.copy(), data = self.cube[:,y,x],
                header=self.header)

    def get_apspec(self, aperture, coordsys=None):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube coordinates)
        """

        import cubes
        if coordsys is not None:
            return pyspeckit.Spectrum(xarr=self.xarr.copy(), data=cubes.extract_aperture( self.cube, aperture , coordsys=coordsys , wcs=self.mapplot.wcs ), header=self.header)
        else:
            return pyspeckit.Spectrum(xarr=self.xarr.copy(), data=cubes.extract_aperture( self.cube, aperture , coordsys=None), header=self.header)

    def set_apspec(self, aperture, coordsys=None):
        """
        Extract an aperture using cubes.extract_aperture
        (defaults to Cube coordinates)
        """

        import cubes
        if coordsys is not None:
            self.data = cubes.extract_aperture( self.cube, aperture , coordsys=coordsys , wcs=self.mapplot.wcs )
        else:
            self.data = cubes.extract_aperture( self.cube, aperture , coordsys=None)

    def fiteach(self, errspec=None, errmap=None, guesses=(), verbose=True,
            verbose_level=1, quiet=True, signal_cut=3, usemomentcube=False,
            blank_value=0, integral=True, direct=False, absorption=False,
            use_nearest_as_guess=False, start_from_point=(0,0), **fitkwargs):
        """
        Fit a spectrum to each valid pixel in the cube

        For guesses, priority is *use_nearest_as_guess*, *usemomentcube*,
        *guesses*, None

        *use_nearest_as_guess* [ False ] 
            Unless the fitted point is the first, it will find the nearest
            other point with a successful fit and use its best-fit parameters
            as the guess

        *start_from_point* [ 'center', (x,y) ]
            Either start from the center or from a point defined by a tuple.
            Work outward from that starting point.  

        *guesses* [ tuple, ndarray[naxis=3] ]
            Either a tuple/list of guesses with len(guesses) = npars or a cube
            of guesses with shape [npars, ny, nx]

        *signal_cut* [ float ]
            Minimum signal-to-noise ratio to "cut" on (i.e., if peak in a given
            spectrum has s/n less than this value, ignore it)

        *blank_value* [ float ]
            Value to replace non-fitted locations with.  A good alternative is
            numpy.nan

        *verbose*       [ bool ]
        *verbose_level* [ int ]
            Controls how much is output.
            0,1 - only changes frequency of updates in loop
            2 - print out messages when skipping pixels
            3 - print out messages when fitting pixels
            4 - specfit will be verbose 

        """

        if not hasattr(self.mapplot,'plane'):
            print "Must mapplot before fitting."
            return

        xx,yy = np.indices(self.mapplot.plane.shape)
        if isinstance(self.mapplot.plane, np.ma.core.MaskedArray): 
            OK = (True-self.mapplot.plane.mask) * self.maskmap
        else:
            OK = np.isfinite(self.mapplot.plane) * self.maskmap

        distance = ((xx)**2 + (yy)**2)**0.5
        if start_from_point == 'center':
            start_from_point = (xx.max()/2., yy.max/2.)
        d_from_start = np.roll( np.roll( distance, start_from_point[0], 0), start_from_point[1], 1)
        sort_distance = np.argsort(d_from_start.flat)

        valid_pixels = zip(xx.flat[sort_distance][OK.flat[sort_distance]],yy.flat[sort_distance][OK.flat[sort_distance]])

        if usemomentcube:
            npars = self.momentcube.shape[0]
        else:
            npars = len(guesses)

        self.parcube = np.zeros((npars,)+self.mapplot.plane.shape)
        self.errcube = np.zeros((npars,)+self.mapplot.plane.shape) 
        if integral: self.integralmap = np.zeros((2,)+self.mapplot.plane.shape)

        # array to store whether pixels have fits
        has_fit = np.zeros(self.mapplot.plane.shape, dtype='bool')

        t0 = time.time()

        for ii,(x,y) in enumerate(valid_pixels):
            sp = self.get_spectrum(x,y)
            if errspec is not None:
                sp.error = errspec
            elif errmap is not None:
                sp.error = np.ones(sp.data.shape) * errmap[y,x]
            else:
                if verbose_level > 1 and ii==0: print "WARNING: using data std() as error."
                sp.error[:] = sp.data.std()
            if sp.error is not None and signal_cut > 0:
                if absorption:
                    max_sn = (-1*sp.data / sp.error).max()
                else:
                    max_sn = (sp.data / sp.error).max()
                if max_sn < signal_cut:
                    if verbose_level > 1:
                        print "Skipped %4i,%4i (s/n=%0.2g)" % (x,y,max_sn)
                    continue
                if verbose_level > 2:
                    print "Fitting %4i,%4i (s/n=%0.2g)" % (x,y,max_sn)
            sp.specfit.Registry = self.Registry # copy over fitter registry
            
            if use_nearest_as_guess and has_fit.sum() > 0:
                if verbose_level > 1 and ii == 0: print "Using nearest fit as guess"
                d = np.roll( np.roll( distance, x, 0), y, 1)
                # If there's no fit, set its distance to be unreasonably large
                nearest_ind = np.argmin(d+1e10*(True-has_fit))
                nearest_x, nearest_y = xx.flat[nearest_ind],yy.flat[nearest_ind]
                gg = self.parcube[:,nearest_y,nearest_x]
            if usemomentcube:
                if verbose_level > 1 and ii == 0: print "Using moment cube"
                gg = self.momentcube[:,y,x]
            elif hasattr(guesses,'shape') and guesses.shape[1:] == self.cube.shape[1:]:
                if verbose_level > 1 and ii == 0: print "Using input guess cube"
                gg = guesses[:,y,x]
            else:
                if verbose_level > 1 and ii == 0: print "Using input guess"
                gg = guesses

            sp.specfit(guesses=gg, quiet=True, verbose=verbose_level>3, **fitkwargs)
            self.parcube[:,y,x] = sp.specfit.modelpars
            self.errcube[:,y,x] = sp.specfit.modelerrs
            if integral: self.integralmap[:,y,x] = sp.specfit.integral(direct=direct,return_error=True)
            has_fit[y,x] = True
        
            if blank_value != 0:
                self.errcube[self.parcube == 0] = blank_value
                self.parcube[self.parcube == 0] = blank_value

            if verbose:
                if ii % (min(10**(3-verbose_level),1)) == 0:
                    print "Finished fit %i.  Elapsed time is %0.1f seconds" % (ii, time.time()-t0)

        # make sure the fitter / fittype are set for the cube
        self.specfit.fitter = sp.specfit.fitter
        self.specfit.fittype = sp.specfit.fittype

        if verbose:
            print "Finished final fit %i.  Elapsed time was %0.1f seconds" % (ii, time.time()-t0)


    def momenteach(self, verbose=True, verbose_level=1, **kwargs):
        """
        Return a cube of the moments of each pixel
        """

        if not hasattr(self.mapplot,'plane'):
            print "Must mapplot before moment-fitting."
            return

        xx,yy = np.indices(self.mapplot.plane.shape)
        if isinstance(self.mapplot.plane, np.ma.core.MaskedArray): 
            OK = (True-self.mapplot.plane.mask) * self.maskmap
        else:
            OK = np.isfinite(self.mapplot.plane) * self.maskmap
        valid_pixels = zip(xx[OK],yy[OK])

        _temp_moment = self.get_spectrum(yy[OK][0],xx[OK][0]).moments(**kwargs)

        self.momentcube = np.zeros((len(_temp_moment),)+self.mapplot.plane.shape)

        t0 = time.time()

        for ii,(x,y) in enumerate(valid_pixels):
            sp = self.get_spectrum(x,y)
            self.momentcube[:,y,x] = sp.moments(**kwargs)
            if verbose:
                if ii % 10**(3-verbose_level) == 0:
                    print "Finished moment %i.  Elapsed time is %0.1f seconds" % (ii, time.time()-t0)

        if verbose:
            print "Finished final moment %i.  Elapsed time was %0.1f seconds" % (ii, time.time()-t0)

    def show_moment(self, momentnumber, **kwargs):
        """
        If moments have been computed, display them in the mapplot window
        """

        if not hasattr(self,'momentcube'):
            print "Compute moments first"
            return

        self.mapplot.plane = self.momentcube[momentnumber,:,:].squeeze()

        self.mapplot(estimator=None, **kwargs)

    def show_fit_param(self, parnumber, **kwargs):
        """
        If pars have been computed, display them in the mapplot window
        """

        if not hasattr(self,'parcube'):
            print "Compute fit parameters first"
            return

        self.mapplot.plane = self.parcube[parnumber,:,:].squeeze()

        self.mapplot(estimator=None, **kwargs)



class CubeStack(Cube):
    """
    The Cube equivalent of Spectra: for stitching multiple cubes with the same
    spatial grid but different frequencies together
    """

    def __init__(self, cubelist, xunits='GHz', x0=0, y0=0, maskmap=None, **kwargs):
        """
        Initialize the Cube.  Accepts files in the following formats:
            - .fits

        x0,y0 - initial spectrum to use (defaults to lower-left corner)
        """

        print "Creating Cube Stack"
        cubelist = list(cubelist)
        for ii,cube in enumerate(cubelist):
            if type(cube) is str:
                cube = Cube(cube)
                cubelist[ii] = cube
            if cube.xarr.units != xunits:
                # convert all inputs to same (non-velocity) units
                cube.xarr.convert_to_unit(xunits, **kwargs)

        self.cubelist = cubelist

        print "Concatenating data"
        self.xarr = spectrum.units.SpectroscopicAxes([sp.xarr for sp in cubelist])
        self.cube = np.ma.concatenate([cube.cube for cube in cubelist])
        if self.cube.mask in (False,np.bool_(False)):
            # mask causes major problems internally for numpy...
            self.cube = np.array(self.cube)
        self._sort()
        self.data = self.cube[:,y0,x0]
        self.error = None

        self.header = cubelist[0].header
        for cube in cubelist:
            for key,value in cube.header.items():
                self.header.update(key,value)

        
        self.units = cubelist[0].units
        for cube in cubelist: 
            if cube.units != self.units: 
                raise ValueError("Mismatched units")

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
        for writer in spectrum.writers.writers: self.writer[writer] = spectrum.writers.writers[writer](self)
        self.mapplot = mapplot.MapPlotter(self)

    def _sort(self):
        """ Sort the data in order of increasing X (could be decreasing, but
        must be monotonic for plotting reasons) """

        indices = self.xarr.argsort()
        self.xarr = self.xarr[indices]
        self.cube = self.cube[indices,:,:]
