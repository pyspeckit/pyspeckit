"""
~~~~~~~~
cubes.py
~~~~~~~~

From `agpy <http://code.google.com/p/agpy/source/browse/trunk/agpy/cubes.py>`_,
contains functions to perform various transformations on data cubes and their
headers.


"""
from __future__ import print_function
from six.moves import xrange
from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array,mean,nansum
from math import acos,atan2,tan
import numpy
import numpy as np
import copy
import os
import astropy.io.fits as fits
import astropy.wcs as pywcs
import tempfile
import warnings
import astropy
from astropy import coordinates
from astropy import log
try:
    from AG_fft_tools import smooth
    smoothOK = True
except ImportError:
    smoothOK = False
try:
    from scipy.interpolate import UnivariateSpline
    scipyOK = True
except ImportError:
    scipyOK = False

from . import posang # agpy code
from ..parallel_map import parallel_map

dtor = pi/180.0


def blfunc_generator(x=None, polyorder=None, splineorder=None,
                     sampling=1):
    """
    Generate a function that will fit a baseline (polynomial or spline) to a
    data set.  Either ``splineorder`` or ``polyorder`` must be set

    Parameters
    ----------
    x : np.ndarray or None
        The X-axis of the fitted array.  Will be set to
        ``np.arange(len(data))`` if not specified
    polyorder : None or int
        The polynomial order.
    splineorder : None or int
    sampling : int
        The sampling rate to use for the data.  Can set to higher numbers to
        effectively downsample the data before fitting
    """
    def blfunc(args, x=x):
        yfit,yreal = args
        if hasattr(yfit,'mask'):
            mask = ~yfit.mask
        else:
            mask = np.isfinite(yfit)

        if x is None:
            x = np.arange(yfit.size, dtype=yfit.dtype)

        ngood = np.count_nonzero(mask)
        if polyorder is not None:
            if ngood < polyorder:
                return yreal
            else:
                endpoint = ngood - (ngood % sampling)
                y = np.mean([yfit[mask][ii:endpoint:sampling]
                             for ii in range(sampling)], axis=0)
                polypars = np.polyfit(x[mask][sampling/2:endpoint:sampling],
                                      y, polyorder)
                return yreal-np.polyval(polypars, x).astype(yreal.dtype)

        elif splineorder is not None and scipyOK:
            if splineorder < 1 or splineorder > 4:
                raise ValueError("Spline order must be in {1,2,3,4}")
            elif ngood <= splineorder:
                return yreal
            else:
                log.debug("splinesampling: {0}  "
                          "splineorder: {1}".format(sampling, splineorder))
                endpoint = ngood - (ngood % sampling)
                y = np.mean([yfit[mask][ii:endpoint:sampling]
                             for ii in range(sampling)], axis=0)
                if len(y) <= splineorder:
                    raise ValueError("Sampling is too sparse.  Use finer sampling or "
                                     "decrease the spline order.")
                spl = UnivariateSpline(x[mask][sampling/2:endpoint:sampling],
                                       y,
                                       k=splineorder,
                                       s=0)
                return yreal-spl(x)
        else:
            raise ValueError("Must provide polyorder or splineorder")

    return blfunc


def baseline_cube(cube, polyorder=None, cubemask=None, splineorder=None,
                  numcores=None, sampling=1):
    """
    Given a cube, fit a polynomial to each spectrum

    Parameters
    ----------
    cube: np.ndarray
        An ndarray with ndim = 3, and the first dimension is the spectral axis
    polyorder: int
        Order of the polynomial to fit and subtract
    cubemask: boolean ndarray
        Mask to apply to cube.  Values that are True will be ignored when
        fitting.
    numcores : None or int
        Number of cores to use for parallelization.  If None, will be set to
        the number of available cores.
    """
    x = np.arange(cube.shape[0], dtype=cube.dtype)
    #polyfitfunc = lambda y: np.polyfit(x, y, polyorder)
    blfunc = blfunc_generator(x=x,
                              splineorder=splineorder,
                              polyorder=polyorder,
                              sampling=sampling)

    reshaped_cube = cube.reshape(cube.shape[0], cube.shape[1]*cube.shape[2]).T

    if cubemask is None:
        log.debug("No mask defined.")
        fit_cube = reshaped_cube
    else:
        if cubemask.dtype != 'bool':
            raise TypeError("Cube mask *must* be a boolean array.")
        if cubemask.shape != cube.shape:
            raise ValueError("Mask shape does not match cube shape")
        log.debug("Masking cube with shape {0} "
                  "with mask of shape {1}".format(cube.shape, cubemask.shape))
        masked_cube = cube.copy()
        masked_cube[cubemask] = np.nan
        fit_cube = masked_cube.reshape(cube.shape[0], cube.shape[1]*cube.shape[2]).T


    baselined = np.array(parallel_map(blfunc, zip(fit_cube,reshaped_cube), numcores=numcores))
    blcube = baselined.T.reshape(cube.shape)
    return blcube



def flatten_header(header,delete=False):
    """
    Attempt to turn an N-dimensional fits header into a 2-dimensional header
    Turns all CRPIX[>2] etc. into new keywords with suffix 'A'

    header must be a fits.Header instance
    """

    if not isinstance(header,fits.Header):
        raise Exception("flatten_header requires a fits.Header instance")

    newheader = header.copy()

    for key in newheader.keys():
        try:
            if delete and int(key[-1]) >= 3 and key[:2] in ['CD','CR','CT','CU','NA']:
                newheader.pop(key)
            elif (int(key[-1]) >= 3 or int(key[2])>=3) and key[:2] in ['CD','CR','CT','CU','NA','PC']:
                newheader.rename_keyword(key,'A'+key,force=True)
            if delete and (int(key[4]) >= 3 or int(key[7]) >= 3) and key[:2]=='PC' and key in newheader:
                newheader.pop(key)
        except ValueError:
            # if key[-1] is not an int
            pass
        except IndexError:
            # if len(key) < 2
            pass
    newheader['NAXIS'] = 2
    if header.get('WCSAXES'):
        newheader['WCSAXES'] = 2

    return newheader

def speccen_header(header, lon=None, lat=None, proj='TAN', system='celestial',
                   spectral_axis=3, celestial_axes=[1,2]):
    """
    Turn a cube header into a spectrum header, retaining RA/Dec vals where possible
    (speccen is like flatten; spec-ify would be better but, specify?  nah)

    Assumes 3rd axis is velocity
    """
    newheader = header.copy()
    new_spectral_axis = 1
    newheader['CRVAL{0}'.format(new_spectral_axis)] = header.get('CRVAL{0}'.format(spectral_axis))
    newheader['CRPIX{0}'.format(new_spectral_axis)] = header.get('CRPIX{0}'.format(spectral_axis))
    if 'CD{0}_{0}'.format(new_spectral_axis) in header:
        newheader.rename_keyword('CD{0}_{0}'.format(new_spectral_axis),
                                 'OLDCD{0}_{0}'.format(new_spectral_axis))
    elif 'CDELT{0}'.format(new_spectral_axis) in header:
        newheader.rename_keyword('CDELT{0}'.format(new_spectral_axis),'OLDCDEL{0}'.format(new_spectral_axis))
    if 'CD{0}_{0}'.format(spectral_axis) in header:
        newheader['CDELT{0}'.format(new_spectral_axis)] = header.get('CD{0}_{0}'.format(spectral_axis))
    elif 'CDELT{0}'.format(spectral_axis) in header:
        newheader['CDELT{0}'.format(new_spectral_axis)] = header.get('CDELT{0}'.format(spectral_axis))
    newheader['CTYPE{0}'.format(new_spectral_axis)] = 'VRAD'
    if header.get('CUNIT{0}'.format(spectral_axis)):
        newheader['CUNIT{0}'.format(new_spectral_axis)] = header.get('CUNIT{0}'.format(spectral_axis))
    else:
        print("Assuming CUNIT3 is km/s in speccen_header")
        newheader['CUNIT{0}'.format(new_spectral_axis)] = 'km/s'
    newheader['CRPIX2'] = 1
    newheader['CRPIX{0}'.format(spectral_axis)] = 1
    if system == 'celestial':
        c2 = 'RA---'
        c3 = 'DEC--'
    elif system == 'galactic':
        c2 = 'GLON-'
        c3 = 'GLAT-'
    elif system == 'PIXEL':
        c2 = 'PIX--'
        c3 = 'PIX--'
    newheader['CTYPE2'] = c2+proj
    newheader['CTYPE{0}'.format(spectral_axis)] = c3+proj

    if lon is not None:
        newheader['CRVAL2'] = lon
    if lat is not None:
        newheader['CRVAL{0}'.format(spectral_axis)] = lat

    if 'CD2_2' in header:
        newheader.rename_keyword('CD2_2','OLDCD2_2')
    if 'CD{0}_{0}'.format(spectral_axis) in header:
        newheader.rename_keyword('CD{0}_{0}'.format(spectral_axis),
                                 'OLDCD{0}_{0}'.format(spectral_axis))
    if 'CROTA2' in header:
        newheader.rename_keyword('CROTA2','OLDCROT2')

    return newheader

def extract_aperture(cube, ap, r_mask=False, wcs=None,
                     coordsys='galactic', wunit='arcsec', debug=False,
                     method='mean'):
    """
    Extract an aperture from a data cube.  E.g. to acquire a spectrum
    of an outflow that is extended.

    Cube should have shape [z,y,x], e.g.
    cube = fits.getdata('datacube.fits')

    Apertures are specified in PIXEL units with an origin of 0,0 (NOT the 1,1
    fits standard!) unless wcs and coordsys are specified

    Parameters
    ----------
    ap : list
        For a circular aperture, len(ap)=3:
            ap = [xcen,ycen,radius]
        For an elliptical aperture, len(ap)=5:
            ap = [xcen,ycen,height,width,PA]
    wcs : wcs
        a pywcs.WCS instance associated with the data cube
    coordsys : str
        the coordinate system the aperture is specified in.
        Options are 'celestial' and 'galactic'.  Default is 'galactic'
    wunit : str
        units of width/height.  default 'arcsec', options 'arcmin' and 'degree'
    method : str
        'mean' or 'sum' (average over spectra, or sum them)
        or 'error' for sqrt(sum-of-squares / n)

    Other Parameters
    ----------------
    r_mask : bool
    return mask in addition to spectrum (for error checking?)
    """
    warnings.warn("SpectralCube can do what subimage_integ does much more easily!",
                  DeprecationWarning)

    if wcs is not None and coordsys is not None:
        if debug:
            print("Converting aperture ",ap,)
        ap = aper_world2pix(ap,wcs,coordsys=coordsys,wunit=wunit)
        if debug:
            print(" to ",ap)

    if len(ap) == 3:
        sh = cube.shape
        yind,xind = indices(sh[1:3]) # recall that python indices are backwards
        dis = sqrt((xind-ap[0])**2+(yind-ap[1])**2)
        mask = dis < ap[2]
    elif len(ap) == 5:
        yinds,xinds = indices(cube.shape[1:3])
        th = (ap[4])*dtor
        xindr = (xinds-ap[0])*cos(th) + (yinds-ap[1])*sin(th)
        yindr = (xinds-ap[0])*-sin(th) + (yinds-ap[1])*cos(th)
        ratio = max(ap[2:4])/min(ap[2:4])
        mask = ((xindr*ratio)**2 + yindr**2)**0.5 < max(ap[2:4])
    else:
        raise Exception("Wrong number of parameters.  Need either 3 parameters "
                        "for a circular aperture or 5 parameters for an "
                        "elliptical aperture.")

    npixinmask = mask.sum()
    if method == 'mean':
        specsum = nansum(cube[:, mask], axis=1)
        spec = specsum / npixinmask
    elif method == 'error':
        specsum = nansum(cube[:, mask]**2, axis=1)
        spec = (specsum)**0.5 / npixinmask
    else:
        specsum = nansum(cube[:, mask], axis=1)

    if r_mask:
        return spec,mask
    else:
        return spec

def integ(file,vrange,xcen=None,xwidth=None,ycen=None,ywidth=None,**kwargs):
    """
    wrapper of subimage_integ that defaults to using the full image
    """
    if isinstance(file,fits.PrimaryHDU):
        header = file.header
        cube = file.data
    elif isinstance(file,fits.HDUList):
        header = file[0].header
        cube = file[0].data
    else:
        file = fits.open(file)
        header = file[0].header
        cube = file[0].data

    if None in [xcen,xwidth,ycen,ywidth]:
        xcen = header['NAXIS1'] / 2
        xwidth = xcen + header['NAXIS1'] % 2
        ycen = header['NAXIS2'] / 2
        ywidth = ycen + header['NAXIS2'] % 2

    return subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,header=header,**kwargs)

def subimage_integ(cube, xcen, xwidth, ycen, ywidth, vrange, header=None,
                   average=mean, dvmult=False, return_HDU=False,
                   units="pixels", zunits=None):
    """
    Returns a sub-image from a data cube integrated over the specified velocity range

    NOTE: With `spectral_cube <spectral-cube.rtfd.org>`_, subcube features can
    be easily applied with the `.subcube` method, and integration is handled
    separately.

    Parameters
    ----------
    cube : np.ndarray
        A 3-dimensional numpy array with dimensions (velocity, y, x)
    xcen,ycen : float
        The center in the X,Y-dimension.  See `units` below for unit information
    xwidth,ywidth : float
        The width in the X,Y-dimension.  See `units` below for unit information
        xwidth and ywidth are "radius" values, i.e. half the length that will be extracted
    vrange : (float,float)
        The velocity range to integrate over.  See `zunits` below for unit information
    header : `astropy.io.fits.Header` or None
        If specified, will allow the use of WCS units
    average : function
        The function to apply when 'integrating' over the subcube
    dvmult : bool
        If dvmult is set, multiply the average by DV (this is useful if you set
        average=sum and dvmul=True to get an integrated value, e.g. K km/s or
        Jy km/s)
    return_hdu : bool
        If specified, will return an HDU object, otherwise will return the
        array and header
    units : 'pixels' or 'wcs'
        If 'pixels', all units (xcen, ycen, xwidth, ywidth) will be in pixels.
        If 'wcs', the values will be converted from WCS units to pixel units
        using the WCS specified by the `header`
    zunits : 'pixels' or 'wcs' or None
        If None, will be set to be the same as `units`

    Returns
    -------
    subim, hdu : tuple
        A tuple (integrated array, header) if ``return_hdu`` is ``False``, or an HDU if
        it is True
    """

    if header:
        flathead = flatten_header(header.copy())
        wcs = pywcs.WCS(header=flathead)
        if header.get('CD3_3'): CD3 = header.get('CD3_3')
        else: CD3 = header.get('CDELT3')

    if units=="pixels":
        xlo = int( max([xcen-xwidth,0])              )
        ylo = int( max([ycen-ywidth,0])              )
        xhi = int( min([xcen+xwidth,cube.shape[2]])  )
        yhi = int( min([ycen+ywidth,cube.shape[1]])  )
    elif units=="wcs" and header:
        newxcen,newycen = wcs.wcs_world2pix(xcen,ycen,0)
        try:
            newxwid,newywid = xwidth / abs(wcs.wcs.cd[0,0]), ywidth / abs(wcs.wcs.cd[1,1])
        except AttributeError:
            newxwid,newywid = xwidth / abs(wcs.wcs.cdelt[0]), ywidth / abs(wcs.wcs.cdelt[1])
        xlo = int( max([newxcen-newxwid,0]) )
        ylo = int( max([newycen-newywid,0]) )
        xhi = int( min([newxcen+newxwid,cube.shape[2]]) )
        yhi = int( min([newycen+newywid,cube.shape[1]]) )
    else:
        print("Can only use wcs if you pass a header.")

    if zunits is None:
        zunits = units
    if zunits == 'pixels':
        zrange = vrange
    if zunits == 'wcs':
        zrange = ( array(vrange)-header.get('CRVAL3') ) / CD3 - 1 + header.get('CRPIX3')

    subim = average(cube[zrange[0]:zrange[1],ylo:yhi,xlo:xhi],axis=0)
    if dvmult and CD3: subim *= CD3
    elif dvmult:
        print("Error: could not multiply by dv; CD3=",CD3)

    if header is None:
        return subim
    else:
        # Cannot set crval2 != 0 for Galactic coordinates: therefore, probably
        # wrong approach in general
        #crv1,crv2 = wcs.wcs_pix2world(xlo,ylo,0)

        #try:
        #    flathead['CRVAL1'] = crv1[0]
        #    flathead['CRVAL2'] = crv2[0]
        #except IndexError:
        #    flathead['CRVAL1'] = crv1.item() # np 0-d arrays are not scalar
        #    flathead['CRVAL2'] = crv2.item() # np 0-d arrays are not scalar

        # xlo, ylo have been forced to integers already above
        flathead['CRPIX1'] = flathead['CRPIX1'] - xlo
        flathead['CRPIX2'] = flathead['CRPIX2'] - ylo

        if return_HDU:
            return fits.PrimaryHDU(data=subim,header=flathead)
        else:
            return subim,flathead

def subcube(cube, xcen, xwidth, ycen, ywidth, header=None,
        dvmult=False, return_HDU=False, units="pixels",
        widthunits="pixels"):
    """
    Crops a data cube

    All units assumed to be pixel units

    cube has dimensions (velocity, y, x)

    xwidth and ywidth are "radius" values, i.e. half the length that will be extracted

    if dvmult is set, multiple the average by DV (this is useful if you set
    average=sum and dvmul=True to get an integrated value)

    """

    if header:
        newheader = header.copy()
        flathead = flatten_header(header.copy())
        wcs = pywcs.WCS(header=flathead)

    if widthunits == "pixels":
        newxwid, newywid = xwidth, ywidth
    elif widthunits == "wcs":
        try:
            newxwid,newywid = xwidth / abs(wcs.wcs.cd[0,0]), ywidth / abs(wcs.wcs.cd[1,1])
        except AttributeError:
            newxwid,newywid = xwidth / abs(wcs.wcs.cdelt[0]), ywidth / abs(wcs.wcs.cdelt[1])
    else:
        raise Exception("widthunits must be either 'wcs' or 'pixels'")

    if units=="pixels":
        newxcen,newycen = xcen,ycen
    elif units=="wcs" and header:
        newxcen,newycen = wcs.wcs_world2pix(xcen,ycen,0)
    else:
        raise Exception("units must be either 'wcs' or 'pixels'")

    x1 = int( numpy.floor( max([newxcen-newxwid,0]) ) )
    y1 = int( numpy.floor( max([newycen-newywid,0]) ) )
    x2 = int( numpy.ceil( min([newxcen+newxwid,cube.shape[2]]) ) )
    y2 = int( numpy.ceil( min([newycen+newywid,cube.shape[1]]) ) )

    xhi = max(x1,x2)
    xlo = min(x1,x2)
    yhi = max(y1,y2)
    ylo = min(y1,y2)

    subim = cube[:,ylo:yhi,xlo:xhi]

    if return_HDU:

        xmid_sky,ymid_sky = wcs.wcs_pix2world(xlo+xwidth,ylo+ywidth,0)

        try:
            newheader['CRVAL1'] = xmid_sky[0]
            newheader['CRVAL2'] = ymid_sky[0]
        except IndexError:
            newheader['CRVAL1'] = float(xmid_sky)
            newheader['CRVAL2'] = float(ymid_sky)
        newheader['CRPIX1'] = 1+xwidth
        newheader['CRPIX2'] = 1+ywidth

        newHDU =  fits.PrimaryHDU(data=subim,header=newheader)
        if newHDU.header.get('NAXIS1') == 0 or newHDU.header.get('NAXIS2') == 0:
            raise Exception("Cube has been cropped to 0 in one dimension")

        return newHDU
    else:
        return subim

def aper_world2pix(ap,wcs,coordsys='galactic',wunit='arcsec'):
    """
    Converts an elliptical aperture (x,y,width,height,PA) from
    WCS to pixel coordinates given an input wcs (an instance
    of the pywcs.WCS class).  Must be a 2D WCS header.


    """
    convopt = {'arcsec':3600.0,'arcmin':60.0,'degree':1.0}
    try:
        conv = convopt[wunit]
    except:
        raise Exception("Must specify wunit='arcsec','arcmin', or 'degree'")

    if len(wcs.wcs.cdelt) != 2:
        raise Exception("WCS header is not strictly 2-dimensional.  Look for 3D keywords.")
    if '' in wcs.wcs.ctype:
        raise Exception("WCS header has no CTYPE.")

    if coordsys.lower() == 'galactic':
        pos = coordinates.SkyCoord(ap[0],ap[1],unit=('deg','deg'), frame='galactic')
    elif coordsys.lower() in ('radec','fk5','icrs','celestial'):
        pos = coordinates.SkyCoord(ap[0],ap[1],unit=('deg','deg'), frame='fk5')

    if wcs.wcs.ctype[0][:2] == 'RA':
        ra,dec = pos.icrs.ra.deg,pos.icrs.dec.deg
    elif wcs.wcs.ctype[0][:4] == 'GLON':
        ra,dec = pos.galactic.l.deg,pos.galactic.b.deg
    else:
        raise Exception("WCS CTYPE has no match.")
    # workaround for a broken wcs.wcs_sky2pix
    try:
        radif = (wcs.wcs.crval[0]-ra)*dtor
        gamma = acos(cos(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)*cos(radif)+sin(dec*dtor)*sin(wcs.wcs.crval[1]*dtor)) / dtor
        theta = atan2( sin(radif) , ( tan(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)-sin(wcs.wcs.crval[1]*dtor)*cos(radif) ) )
        x = -gamma * sin(theta) / wcs.wcs.cd[0,0] + wcs.wcs.crpix[0]
        y = gamma * cos(theta) / wcs.wcs.cd[1,1] + wcs.wcs.crpix[1]
    except:
        radif = (wcs.wcs.crval[0]-ra)*dtor
        gamma = acos(cos(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)*cos(radif)+sin(dec*dtor)*sin(wcs.wcs.crval[1]*dtor)) / dtor
        theta = atan2( sin(radif) , ( tan(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)-sin(wcs.wcs.crval[1]*dtor)*cos(radif) ) )
        x = -gamma * sin(theta) / wcs.wcs.cdelt[0] + wcs.wcs.crpix[0]
        y = gamma * cos(theta) / wcs.wcs.cdelt[1] + wcs.wcs.crpix[1]

    #print "DEBUG: x,y from math (vectors): ",x,y
    #x,y = wcs.wcs_world2pix(ra,dec,0)  # convert WCS coordinate to pixel coordinate (0 is origin, do not use fits convention)
    #print "DEBUG: x,y from wcs: ",x,y
    try:
        x=x[0] - 1 # change from FITS to python convention
        y=y[0] - 1 # change from FITS to python convention
        #print "DEBUG: x,y from math: ",x,y
    except:
        pass
    # cd is default, cdelt is backup
    if len(ap) > 3:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cd[0,0])
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cdelt[0])
        apold = copy.copy(ap)
        if len(ap) == 5:
            PA = ap[4]
            ap = [x,y,width,height,PA]
        else:
            ap = [x,y,width,height]
    elif len(ap) == 3:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
        apold = copy.copy(ap)
        ap = [x,y,width]
    else:
        raise TypeError("Aperture length is incorrect.")

    return ap


def getspec(lon,lat,rad,cube,header,r_fits=True,inherit=True,wunit='arcsec'):
    """
    Given a longitude, latitude, aperture radius (arcsec), and a cube file,
    return a .fits file or a spectrum.

    Parameters
    ----------
    lon: float
    lat: float
        longitude and latitude center of a circular aperture in WCS coordinates
        must be in coordinate system of the file
    rad: float
        radius (default degrees) of aperture
    """

    convopt = {'arcsec':1.0,'arcmin':60.0,'degree':3600.0}

    flathead = flatten_header(header)
    wcs = pywcs.WCS(flathead)
    if wcs.wcs.ctype[0][:2] == 'RA':
      coordsys='celestial'
    elif wcs.wcs.ctype[0][:4] == 'GLON':
      coordsys='galactic'
    spec = extract_aperture(cube,[lon,lat,rad],wcs=wcs,
            coordsys=coordsys,wunit=wunit)

    if nansum(spec) == 0:
        print("Total of extracted spectrum was zero. lon,lat,rad: ",lon,lat,rad)
        #import pdb; pdb.set_trace()

    if r_fits:
        if inherit:
            newhead = header.copy()
        else:
            newhead = fits.Header()
        try:
            newhead['CD1_1'] = header['CD3_3']
        except KeyError:
            newhead['CD1_1'] = header['CDELT3']
        newhead['CRPIX1'] = header['CRPIX3']
        newhead['CRVAL1'] = header['CRVAL3']
        try:
            newhead['CTYPE1'] = header['CTYPE3']
        except KeyError:
            newhead['CTYPE1'] = "VRAD"
        try:
            newhead['CUNIT1'] = header['CUNIT3']
        except KeyError:
            print("Header did not contain CUNIT3 keyword.  Defaulting to km/s")
            newhead['CUNIT1'] = "km/s"
        newhead['BUNIT'] = header['BUNIT']
        newhead['APGLON'] = lon
        newhead['APGLAT'] = lat
        newhead['APRAD'] = (rad*convopt[wunit],'arcseconds') # radius in arcsec
        newfile = fits.PrimaryHDU(data=spec,header=newhead)
        return newfile
    else:
        return spec

def getspec_reg(cubefilename,region,**kwargs):
    """
    Aperture extraction from a cube using a pyregion circle region

    The region must be in the same coordinate system as the cube header

    .. warning:: The second argument of getspec_reg requires a pyregion region list,
        and therefore this code depends on `pyregion`_.
    """

    ds9tocoords = {'fk5':'celestial','galactic':'galactic','icrs':'celestial'}

    if region.name != 'circle':
        raise Exception("Only circular apertures are implemented so far")

    l,b,r = region.coord_list
    #pos = coords.Position([l,b],system=ds9tocoords[region.coord_format])
    if isinstance(cubefilename,fits.HDUList):
        cubefile = cubefilename
    else:
        cubefile = fits.open(cubefilename)
    header = cubefile[0].header
    cube = cubefile[0].data
    if len(cube.shape) == 4: cube = cube[0,:,:,:]

    sp = getspec(l,b,r,cube,header,wunit='degree',**kwargs)

    return sp

def coords_in_image(fitsfile,lon,lat,system='galactic'):
    """
    Determine whether the coordinates are inside the image
    """
    if not isinstance(fitsfile,fits.HDUList):
        fitsfile = fits.open(fitsfile)

    wcs = pywcs.WCS(flatten_header(fitsfile[0].header))

    if 'RA' in wcs.wcs.ctype[0]:
        pos = coordinates.Position((lon,lat),system=system)
        lon,lat = pos.j2000()
    if 'GLON' in wcs.wcs.ctype[0]:
        pos = coordinates.Position((lon,lat),system=system)
        lon,lat = pos.galactic()

    x,y = wcs.wcs_world2pix(lon,lat,0)
    #DEBUG print x,y,wcs.naxis1,wcs.naxis2
    if (0 < x < wcs.naxis1) and (0 < y < wcs.naxis2):
        return True
    else:
        return False

def spectral_smooth(cube, smooth_factor, downsample=True, parallel=True,
                    numcores=None, **kwargs):
    """
    Smooth the cube along the spectral direction
    """

    yy,xx = numpy.indices(cube.shape[1:])

    if downsample:
        newshape = cube[::smooth_factor,:,:].shape
    else:
        newshape = cube.shape

    # need to make the cube "flat" along dims 1&2 for iteration in the "map"
    flatshape = (cube.shape[0],cube.shape[1]*cube.shape[2])

    Ssmooth = lambda x: smooth.smooth(x, smooth_factor, downsample=downsample, **kwargs)
    if parallel:
        newcube = numpy.array(parallel_map(Ssmooth, cube.reshape(flatshape).T, numcores=numcores)).T.reshape(newshape)
    else:
        newcube = numpy.array(map(Ssmooth, cube.reshape(flatshape).T)).T.reshape(newshape)

    #naive, non-optimal version
    # for (x,y) in zip(xx.flat,yy.flat):
    #     newcube[:,y,x] = smooth.smooth(cube[:,y,x], smooth_factor,
    #             downsample=downsample, **kwargs)

    return newcube

def plane_smooth(cube,cubedim=0,parallel=True,numcores=None,**kwargs):
    """
    parallel-map the smooth function

    Parameters
    ----------
    parallel: bool
        defaults True.  Set to false if you want serial (for debug purposes?)
    numcores: int
        pass to parallel_map (None = use all available)
    """
    if not smoothOK:
        return

    if cubedim != 0:
        cube = cube.swapaxes(0,cubedim)

    cubelist = [cube[ii,:,:] for ii in xrange(cube.shape[0])]

    Psmooth = lambda C: smooth.smooth(C,**kwargs)

    if parallel:
        smoothcube = array(parallel_map(Psmooth,cubelist,numcores=numcores))
    else:
        smoothcube = array(map(Psmooth,cubelist))

    if cubedim != 0:
        smoothcube = smoothcube.swapaxes(0,cubedim)

    return smoothcube


try:
    import montage

    def rotcrop_cube(x1, y1, x2, y2, cubename, outname, xwidth=25, ywidth=25,
                     in_system='galactic',  out_system='equatorial',
                     overwrite=True, newheader=None, xcen=None, ycen=None):
        """
        Crop a data cube and then rotate it with montage

        """

        cubefile = fits.open(cubename)

        if xcen is None and ycen is None:
            pos1 = coordinates.Position([x1,y1],system=in_system)
            pos2 = coordinates.Position([x2,y2],system=in_system)

            if cubefile[0].header.get('CTYPE1')[:2] == 'RA':
                x1,y1 = pos1.j2000()
                x2,y2 = pos2.j2000()
                coord_system = 'celestial'
            elif cubefile[0].header.get('CTYPE1')[:4] == 'GLON':
                x1,y1 = pos1.galactic()
                x2,y2 = pos2.galactic()
                coord_system = 'galactic'

            xcen = (x1+x2)/2.0
            ycen = (y1+y2)/2.0
            print(xcen,ycen,xwidth,ywidth,coord_system)
        else:
            coord_system = in_system

        sc = subcube(cubefile[0].data, xcen, xwidth, ycen, ywidth,
                widthunits='pixels', units="wcs", header=cubefile[0].header,
                return_HDU=True)
        # note: there should be no security risk here because fits' writeto
        # will not overwrite by default
        tempcube = tempfile.mktemp(suffix='.fits')
        sc.writeto(tempcube)

        pa = posang.posang(x1,y1,x2,y2,system=coord_system) - 90

        if newheader is None:
            newheader = sc.header.copy()
            cd11 = newheader.get('CDELT1') if newheader.get('CDELT1') else newheader.get('CD1_1')
            cd22 = newheader.get('CDELT2') if newheader.get('CDELT2') else newheader.get('CD2_2')
            cd12 = newheader.get('CD1_2') if newheader.get('CD1_2') else 0.0
            cd21 = newheader.get('CD2_1') if newheader.get('CD2_1') else 0.0
            cdelt = numpy.sqrt(cd11**2+cd12**2)

            tempheader = tempfile.mktemp(suffix='.hdr')
            ycensign = "+" if numpy.sign(ycen) >= 0 else "-"
            montage.mHdr("%s %1s%s" % (xcen, ycensign, numpy.abs(ycen)), xwidth*cdelt,
                    tempheader, system=out_system, height=ywidth*cdelt,
                    pix_size=cdelt*3600.0, rotation=pa)
            os.system("sed -i bck '/END/d' %s" % (tempheader))
            newheader2 = fits.Header()
            newheader2.fromTxtFile(tempheader)
            #newheader2.fromtextfile(tempheader)
            for key in ('CRPIX3','CRVAL3','CDELT3','CD3_3','CUNIT3','WCSTYPE3','CTYPE3'):
                if newheader.get(key):
                    newheader2[key] = newheader.get(key)
            if newheader.get('CD3_3') and newheader2.get('CDELT3') is None:
                newheader2['CDELT3'] = newheader.get('CD3_3')
            if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
                newheader2.toTxtFile(tempheader,overwrite=True)
            else:
                newheader2.toTxtFile(tempheader,clobber=True)
            #if newheader2.get('CDELT3') is None:
            #    raise Exception("No CD3_3 or CDELT3 in header.")
        else:
            if isinstance(newheader,str):
                newheader2 = fits.Header()
                newheader2.fromTxtFile(newheader)
            tempheader = tempfile.mktemp(suffix='.hdr')
            if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
                newheader2.toTxtFile(tempheader,overwrite=True)
            else:
                newheader2.toTxtFile(tempheader,clobber=True)


        montage.wrappers.reproject_cube(tempcube,outname,header=tempheader,clobber=overwrite)
        #print "\n",outname
        #os.system('imhead %s | grep CDELT' % outname)

        # AWFUL hack because montage removes CDELT3
        tempcube = fits.open(outname)
        tempcube.header = newheader2
        #if tempcube.header.get('CDELT3') is None:
        #    raise Exception("No CD3_3 or CDELT3 in header.")
        #print tempcube.header.get('CDELT3')
        if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
            tempcube.writeto(outname,overwrite=True)
        else:
            tempcube.writeto(outname,clobber=True)
        #print tempcube.get('CDELT3')
        #print "\n",outname
        #os.system('imhead %s | grep CDELT' % outname)


        return

    def resample_cube(cubefilename, header):
        inhdr = fits.getheader(cubefilename)

except:
    pass
