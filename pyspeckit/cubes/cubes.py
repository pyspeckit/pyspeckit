"""
~~~~~~~~
cubes.py
~~~~~~~~

From `agpy <http://code.google.com/p/agpy/source/browse/trunk/agpy/cubes.py>`_,
contains functions to perform various transformations on data cubes and their
headers.  

########

"""
from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array,mean,sum,nansum
from math import acos,atan2,tan
import numpy
import copy
import os
import pyfits
import tempfile
import posang # agpy code
import pyspeckit
from pyspeckit.specwarnings import warn
try:
    from AG_fft_tools.convolve import smooth
    from pyspeckit.parallel_map import parallel_map
    smoothOK = True
except ImportError:
    smoothOK = False
try:
    import coords
except ImportError:
    warn( "cubes.py requires coords for aper_world2pix and coords_in_image" )
try:
    import pyregion
except ImportError:
    warn( "cubes.py requires pyregion for getspec_reg" )
try:
    import astropy.wcs as pywcs
except ImportError:
    try:
        import pywcs
    except ImportError:
        warn( "cubes.py requires astropy.wcs or pywcs for some subimage_integ,aper_wordl2pix,getspec, and coords_in_image" )

dtor = pi/180.0

def flatten_header(header):
    """
    Attempt to turn an N-dimensional fits header into a 2-dimensional header
    Turns all CRPIX[>2] etc. into new keywords with suffix 'A'

    header must be a pyfits.Header instance
    """

    if not isinstance(header,pyfits.Header):
        raise Exception("flatten_header requires a pyfits.Header instance")

    newheader = header.copy()

    for key in newheader.keys():
        try:
            if int(key[-1]) >= 3 and key[:2] in ['CD','CR','CT','CU','NA']:
                newheader.rename_key(key,'A'+key,force=True)
        except ValueError:
            # if key[-1] is not an int
            pass
        except IndexError:
            # if len(key) < 2
            pass
    newheader.update('NAXIS',2)

    return newheader

def speccen_header(header,lon=None,lat=None):
    """
    Turn a cube header into a spectrum header, retaining RA/Dec vals where possible
    (speccen is like flatten; spec-ify would be better but, specify?  nah)

    Assumes 3rd axis is velocity
    """
    newheader = header.copy()
    newheader.update('CRVAL1',header.get('CRVAL3'))
    newheader.update('CRPIX1',header.get('CRPIX3'))
    if 'CD1_1' in header: newheader.rename_key('CD1_1','OLDCD1_1')
    elif 'CDELT1' in header: newheader.rename_key('CDELT1','OLDCDEL1')
    if 'CD3_3' in header: newheader.update('CDELT1',header.get('CD3_3'))
    elif 'CDELT3' in header: newheader.update('CDELT1',header.get('CDELT3'))
    newheader.update('CTYPE1','VRAD')
    if header.get('CUNIT3'): newheader.update('CUNIT1',header.get('CUNIT3'))
    else: 
        print "Assuming CUNIT3 is km/s in speccen_header"
        newheader.update('CUNIT1','km/s')
    newheader.update('CRPIX2',1)
    newheader.update('CTYPE2','RA---TAN')
    newheader.update('CRPIX3',1)
    newheader.update('CTYPE3','DEC--TAN')

    if lon is not None: newheader.update('CRVAL2',lon)
    if lat is not None: newheader.update('CRVAL3',lat)

    if 'CD2_2' in header: newheader.rename_key('CD2_2','OLDCD2_2')
    if 'CD3_3' in header: newheader.rename_key('CD3_3','OLDCD3_3')

    return newheader

def extract_aperture(cube,ap,r_mask=False,wcs=None,coordsys='galactic',wunit='arcsec',debug=False):
    """
    Extract an aperture from a data cube.  E.g. to acquire a spectrum
    of an outflow that is extended.

    Cube should have shape [z,y,x], e.g. 
    cube = pyfits.getdata('datacube.fits')

    Apertures are specified in PIXEL units with an origin of 0,0 (NOT the 1,1
    fits standard!) unless wcs and coordsys are specified
    
    INPUTS:
        wcs - a pywcs.WCS instance associated with the data cube
        coordsys - the coordinate system the aperture is specified in.
            Options are 'celestial' and 'galactic'.  Default is 'galactic'
        wunit - units of width/height.  default 'arcsec', options 'arcmin' and 'degree'

    For a circular aperture, len(ap)=3:
        ap = [xcen,ycen,radius]

    For an elliptical aperture, len(ap)=5:
        ap = [xcen,ycen,height,width,PA]

    Optional inputs:
        r_mask - return mask in addition to spectrum (for error checking?)
    """

    if wcs is not None and coordsys is not None:
        if debug: print "Converting aperture ",ap,
        ap = aper_world2pix(ap,wcs,coordsys=coordsys,wunit=wunit)
        if debug: print " to ",ap

    if len(ap) == 3:
        sh = cube.shape
        yind,xind = indices(sh[1:3]) # recall that python indices are backwards
        dis = sqrt((xind-ap[0])**2+(yind-ap[1])**2)
        mask = dis < ap[2]
    elif len(ap) == 5:
        yinds,xinds = indices(cube.shape[1:3])
        th = (ap[4])*dtor
        xindr = (xinds-ap[0])*cos(th)  + (yinds-ap[1])*sin(th)
        yindr = (xinds-ap[0])*-sin(th) + (yinds-ap[1])*cos(th)
        ratio = max(ap[2:4])/min(ap[2:4])
        mask = sqrt( (xindr*ratio)**2 + yindr**2) < max(ap[2:4])
    else:
        raise Exception("Wrong number of parameters.  Need either 3 parameters "+
                "for a circular aperture or 5 parameters for an elliptical "+ 
                "aperture.")

    npixinmask = mask.sum()
    mask3d = repeat(mask[newaxis,:,:],cube.shape[0],axis=0)
    spec = nansum(nansum((cube*mask3d),axis=2),axis=1) / npixinmask

    if r_mask:
        return spec,mask
    else:
        return spec

def integ(file,vrange,xcen=None,xwidth=None,ycen=None,ywidth=None,**kwargs):
    """
    wrapper of subimage_integ that defaults to using the full image
    """
    if isinstance(file,pyfits.PrimaryHDU):
        header = file.header
        cube = file.data
    elif isinstance(file,pyfits.HDUList):
        header = file[0].header
        cube = file[0].data
    else:
        file = pyfits.open(file)
        header = file[0].header
        cube = file[0].data

    if None in [xcen,xwidth,ycen,ywidth]:
        xcen = header['NAXIS1'] / 2
        xwidth = xcen
        ycen = header['NAXIS2'] / 2
        ywidth = ycen

    return subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,header=header,**kwargs)

def subimage_integ(cube, xcen, xwidth, ycen, ywidth, vrange, header=None,
        average=mean, dvmult=False, return_HDU=False, units="pixels",
        zunits=None):
    """
    Returns a sub-image from a data cube integrated over the specified velocity range

    All units assumed to be pixel units

    cube has dimensions (velocity, y, x)

    xwidth and ywidth are "radius" values, i.e. half the length that will be extracted

    if dvmult is set, multiple the average by DV (this is useful if you set
    average=sum and dvmul=True to get an integrated value)

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
        newxcen,newycen = wcs.wcs_sky2pix(xcen,ycen,0)
        try:
            newxwid,newywid = xwidth / abs(wcs.wcs.cd[0,0]), ywidth / abs(wcs.wcs.cd[1,1])
        except AttributeError:
            newxwid,newywid = xwidth / abs(wcs.wcs.cdelt[0]), ywidth / abs(wcs.wcs.cdelt[1])
        xlo = int( max([newxcen-newxwid,0]) )
        ylo = int( max([newycen-newywid,0]) )
        xhi = int( min([newxcen+newxwid,cube.shape[2]]) )
        yhi = int( min([newycen+newywid,cube.shape[1]]) )
    else:
        print "Can only use wcs if you pass a header."

    if zunits is None:
        zunits = units
    if zunits == 'pixels':
        zrange = vrange
    if zunits == 'wcs':
        zrange = ( array(vrange)-header.get('CRVAL3') ) / CD3 - 1 + header.get('CRPIX3')

    subim = average(cube[zrange[0]:zrange[1],ylo:yhi,xlo:xhi],axis=0)
    if dvmult and CD3: subim *= CD3
    elif dvmult: print "Error: could not multiply by dv; CD3=",CD3

    if header is None:
        return subim
    else:
        crv1,crv2 = wcs.wcs_pix2sky(xlo,ylo,0)

        flathead['CRVAL1'] = crv1[0]
        flathead['CRVAL2'] = crv2[0]
        flathead['CRPIX1'] = 1
        flathead['CRPIX2'] = 1
        
        if return_HDU:
            return pyfits.PrimaryHDU(data=subim,header=flathead)
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
        newxcen,newycen = wcs.wcs_sky2pix(xcen,ycen,0)
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

        xmid_sky,ymid_sky = wcs.wcs_pix2sky(xlo+xwidth,ylo+ywidth,0)

        newheader.update('CRVAL1',xmid_sky[0])
        newheader.update('CRVAL2',ymid_sky[0])
        newheader.update('CRPIX1',1+xwidth)
        newheader.update('CRPIX2',1+ywidth)
        
        newHDU =  pyfits.PrimaryHDU(data=subim,header=newheader)
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
    pos = coords.Position((ap[0],ap[1]),system=coordsys)
    if wcs.wcs.ctype[0][:2] == 'RA':
        ra,dec = pos.j2000()
        corrfactor = cos(dec*dtor)
    elif wcs.wcs.ctype[0][:4] == 'GLON':
        ra,dec = pos.galactic()
        corrfactor=1
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
    #x,y = wcs.wcs_sky2pix(ra,dec,0)  # convert WCS coordinate to pixel coordinate (0 is origin, do not use fits convention)
    #print "DEBUG: x,y from wcs: ",x,y
    try:
        x=x[0] - 1 # change from FITS to python convention
        y=y[0] - 1 # change from FITS to python convention
        #print "DEBUG: x,y from math: ",x,y
    except:
        pass
    # cd is default, cdelt is backup
    if len(ap) == 5:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cd[0,0])
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cdelt[0])
        PA = ap[4] 
        apold = copy.copy(ap)
        ap = [x,y,width,height,PA]
    elif len(ap) == 3:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
        apold = copy.copy(ap)
        ap = [x,y,width]

    return ap


def getspec(lon,lat,rad,cube,header,r_fits=True,inherit=True,wunit='arcsec'):
    """
    Given a longitude, latitude, aperture radius (arcsec), and a cube file,
    return a .fits file or a spectrum.
    
    lon,lat - longitude and latitude center of a circular aperture in WCS coordinates
                must be in coordinate system of the file
    rad     - radius (default degrees) of aperture
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
        print "Total of extracted spectrum was zero. lon,lat,rad: ",lon,lat,rad #  Tracing to find your problem."
        #import pdb; pdb.set_trace()

    if r_fits:
        if inherit:
            newhead = header.copy()
        else:
            newhead = pyfits.Header()
        try:
            newhead.update('CD1_1',header['CD3_3'])
        except KeyError:
            newhead.update('CD1_1',header['CDELT3'])
        newhead.update('CRPIX1',header['CRPIX3'])
        newhead.update('CRVAL1',header['CRVAL3'])
        try:
            newhead.update('CTYPE1',header['CTYPE3'])
        except KeyError:
            newhead.update('CTYPE1',"VRAD")
        try:
            newhead.update('CUNIT1',header['CUNIT3'])
        except KeyError:
            print "Header did not contain CUNIT3 keyword.  Defaulting to km/s"
            newhead.update('CUNIT1',"km/s")
        newhead.update('BUNIT',header['BUNIT'])
        newhead.update('APGLON',lon)
        newhead.update('APGLAT',lat)
        newhead.update('APRAD',rad*convopt[wunit],comment='arcseconds') # radius in arcsec
        newfile = pyfits.PrimaryHDU(data=spec,header=newhead)
        return newfile
    else:
        return spec

def getspec_reg(cubefilename,region,**kwargs):
    """
    Aperture extraction from a cube using a pyregion circle region

    The region must be in the same coordinate system as the cube header
    """

    ds9tocoords = {'fk5':'celestial','galactic':'galactic','icrs':'celestial'}

    if region.name != 'circle':
        raise Exception("Only circular apertures are implemented so far")

    l,b,r = region.coord_list
    #pos = coords.Position([l,b],system=ds9tocoords[region.coord_format])
    if isinstance(cubefilename,pyfits.HDUList):
        cubefile = cubefilename
    else:
        cubefile = pyfits.open(cubefilename)
    header = cubefile[0].header
    cube = cubefile[0].data
    if len(cube.shape) == 4: cube = cube[0,:,:,:]

    sp = getspec(l,b,r,cube,header,wunit='degree',**kwargs)

    return sp

def coords_in_image(fitsfile,lon,lat,system='galactic'):
    """
    Determine whether the coordinates are inside the image
    """
    if not isinstance(fitsfile,pyfits.HDUList):
        fitsfile = pyfits.open(fitsfile)

    wcs = pywcs.WCS(flatten_header(fitsfile[0].header))

    if 'RA' in wcs.wcs.ctype[0]:
        pos = coords.Position((lon,lat),system=system)
        lon,lat = pos.j2000()
    if 'GLON' in wcs.wcs.ctype[0]:
        pos = coords.Position((lon,lat),system=system)
        lon,lat = pos.galactic()

    x,y = wcs.wcs_sky2pix(lon,lat,0)
    #DEBUG print x,y,wcs.naxis1,wcs.naxis2
    if (0 < x < wcs.naxis1) and (0 < y < wcs.naxis2):
        return True
    else:
        return False

def spectral_smooth(cube, smooth_factor, **kwargs):
    """
    Smooth the cube along the spectral direction
    """

    newcube = numpy.zeros(cube.shape)
    yy,xx = numpy.indices(cube.shape[1:])

    for (x,y) in zip(xx.flat,yy.flat):
        newcube[:,y,x] = pyspeckit.smooth.smooth(cube[:,y,x], smooth_factor, **kwargs)

    return newcube

def plane_smooth(cube,cubedim=0,parallel=True,numcores=None,**kwargs):
    """
    parallel-map the smooth function

    parallel - defaults True.  Set to false if you want serial (for debug
        purposes?)
    numcores - pass to parallel_map (None = use all available)
    """
    if not smoothOK:
        return

    if cubedim != 0:
        cube = cube.swapaxes(0,cubedim)

    cubelist = [cube[ii,:,:] for ii in xrange(cube.shape[0])]

    Psmooth = lambda C: smooth(C,**kwargs)

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
            in_system='galactic',  out_system='equatorial', clobber=True):
        """
        Crop a data cube and then rotate it with montage

        """

        cubefile = pyfits.open(cubename)

        pos1 = coords.Position([x1,y1],system=in_system)
        pos2 = coords.Position([x2,y2],system=in_system)

        if cubefile[0].header.get('CTYPE1')[:2] == 'RA':
            x1,y1 = pos1.j2000()
            x2,y2 = pos2.j2000()
            coord_system = 'celestial'
        elif  cubefile[0].header.get('CTYPE1')[:2] == 'GLON':
            x1,y1 = pos1.galactic()
            x2,y2 = pos2.galactic()
            coord_system = 'galactic'

        xcen = (x1+x2)/2.0
        ycen = (y1+y2)/2.0

        sc = subcube(cubefile[0].data, xcen, xwidth, ycen, ywidth, 
                widthunits='pixels', units="wcs", header=cubefile[0].header,
                return_HDU=True)
        # note: there should be no security risk here because pyfits' writeto
        # will not overwrite by default
        tempcube = tempfile.mktemp(suffix='.fits')
        sc.writeto(tempcube)
        
        pa = posang.posang(x1,y1,x2,y2,system=coord_system) - 90

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
        newheader2 = pyfits.Header()
        newheader2.fromTxtFile(tempheader)
        for key in ('CRPIX3','CRVAL3','CDELT3','CD3_3','CUNIT3','WCSTYPE3','CTYPE3'):
            if newheader.get(key):
                newheader2.update(key,newheader.get(key))
        if newheader.get('CD3_3') and newheader2.get('CDELT3') is None:
            newheader2.update('CDELT3',newheader.get('CD3_3'))
        newheader2.toTxtFile(tempheader,clobber=True)
        #if newheader2.get('CDELT3') is None:
        #    raise Exception("No CD3_3 or CDELT3 in header.")

        montage.wrappers.reproject_cube(tempcube,outname,header=tempheader,clobber=clobber)
        #print "\n",outname
        #os.system('imhead %s | grep CDELT' % outname)

        # AWFUL hack because montage removes CDELT3
        tempcube = pyfits.open(outname)
        tempcube.header = newheader2
        #if tempcube.header.get('CDELT3') is None:
        #    raise Exception("No CD3_3 or CDELT3 in header.")
        #print tempcube.header.get('CDELT3')
        tempcube.writeto(outname,clobber=True)
        #print tempcube.get('CDELT3')
        #print "\n",outname
        #os.system('imhead %s | grep CDELT' % outname)

        #print "\nnewheader2"
        #print newheader2.ascard
        #print
        
        return

except:
    pass
