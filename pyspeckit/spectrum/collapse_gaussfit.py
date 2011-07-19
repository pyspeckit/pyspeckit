try:
    import scipy
    from scipy import optimize,sqrt
    from scipy.optimize import leastsq
    #from scipy.stats.stats import nanmedian,nanmean,_nanmedian
except ImportError:
    print "Scipy cold not be loaded.  Collapse_gaussfit may fail"
import numpy
from numpy import vectorize,zeros,exp,median,where,asarray,array,nonzero,ma,arange,square
import matplotlib
#matplotlib.use('Agg')
from pylab import indices,figure,clf,savefig,plot,legend,text,axes,title
import pickle
import pyfits
import time
from mad import MAD
from ratosexagesimal import ratos,dectos

def nanmedian(arr):
    """ nanmedian - this version is NOT capable of broadcasting (operating along axes) """
    return median(arr[arr==arr])
def nanmean(arr):
    """ nanmean - this version is NOT capable of broadcasting (operating along axes) """
    return (arr[arr==arr]).mean()

# read in file
# filename = sys.argv[1]
# fitsfile = pyfits.open(filename)
# cube = fitsfile[0].data

# def gaussian(dx,sigma):
#     return lambda x: exp( - (x-dx)**2 / sigma**2 )
# def return_param(xarr,param):
#     errorfunction = lambda p:gaussian(*p)(*indices(xarr.shape))-xarr
#     pars, cov, infodict, errmsg, success = optimize.leastsq(errorfunction, [len(xarr)/2.,1], full_output=1)
#     print errmsg
#     if param == 'width':
#         return pars[1]
#     elif param == 'center':
#         return pars[0]
#     else:
#         return 
def gaussian(dx,sigma,a):
    return lambda x: a*exp( - (x-dx)**2 / sigma**2 )
def double_gaussian(dx1,dx2,sigma1,sigma2,a1,a2):
    return lambda x: a1*exp( - (x-dx1)**2 / sigma1**2 ) +  a2*exp( - (x-dx2)**2 / sigma2**2 )
def triple_gaussian(dx1,dx2,dx3,sigma1,sigma2,sigma3,a1,a2,a3):
    return lambda x: abs(a1)*exp( - (x-dx1)**2 / sigma1**2 ) +  abs(a2)*exp( - (x-dx2)**2 / sigma2**2 ) +  abs(a3)*exp( - (x-dx3)**2 / sigma3**2 )
def n_gaussian(dx,sigma,a):
    def g(x):
        v = zeros(len(x))
        for i in range(len(dx)):
            v += a[i] * exp( - ( x - dx[i] )**2 / sigma[i]**2 )
        return v
    return g
def gerr(xarr):
    return lambda p:xarr-gaussian(*p)(*indices(xarr.shape))
def double_gerr(xarr):
    return lambda p:xarr-double_gaussian(*p)(*indices(xarr.shape))
def triple_gerr(xarr):
    return lambda p:xarr-triple_gaussian(*p)(*indices(xarr.shape))
def return_param(xarr,params=None):
    if params == None:
        params = [xarr.argmax(),5,xarr.max()]
    pars, cov, infodict, errmsg, success = optimize.leastsq(gerr(xarr), params, full_output=1)
    return pars
def return_double_param(xarr,params=None):
    if params == None:
        params = [xarr.argmax(),xarr.argmax()+3,4.2,2.3,xarr.max(),xarr.max()/2]
    pars, cov, infodict, errmsg, success = optimize.leastsq(double_gerr(xarr), params, full_output=1)
    return pars
def return_triple_param(xarr,params=None):
    """
    input parameters: center[1-3],width[1-3],amplitude[1-3]
    """
    if params == None:
        params = [xarr.argmax(),xarr.argmax()+3,xarr.argmax(),4.2,2.3,10,xarr.max(),xarr.max()/2.,xarr.max()/5.]
    pars, cov, infodict, errmsg, success = optimize.leastsq(triple_gerr(xarr), params, full_output=1)
    return pars



def adaptive_collapse_gaussfit(cube,axis=2,nsig=3,nrsig=4,prefix='interesting',
        vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x,doplot=True):
    """
    Attempts to fit one or two Gaussians to each spectrum in a data cube and returns the parameters of the fits.
    Adaptively determines where to fit two Gaussian components based on residuals.  Will fit 3 gaussians if a
    two-gaussian fit is not better than a certain threshold (specified by nsig), and those fits will be output
    to images with filename prefix+(coordinate).png.  The 3-gaussian fit parameters will not be returned because
    the automated fitting is very unlikely to get that part right.

    inputs:
    cube - a data cube with two spatial and one spectral dimensions
    axis - the axis of the spectral dimension
    nsig - number of sigma over the mean residual to trigger double-gaussian fitting
           also, cutoff to do any fitting at all
    prefix - the prefix (including directory name) of the output images from 3-gaussian fitting
    doplot - option to turn off plotting of triple-gaussian fits

    vconv,xtora,ytodec - functions to convert the axes from pixel coordinates to ra/dec/velocity coordinates

    returns:
    width_arr1,width_arr2,chi2_arr,offset_arr1,offset_arr2,amp_arr1,amp_arr2
    The Gaussian widths, line centers (in pixel units), amplitudes, and the chi-squared value, not in that order
    These returns are identical to the returns from double_gaussian, but all components will be zero for the second
    gaussian in the case of a single-gaussian fit

    the triple gaussian is guessed to be the double gaussian plus a broad, low-amplitude gaussian.  Ideally this should
    fit outflows reasonably well, but who knows if it really will.
    Another option is to fit a negative-amplitude gaussian to account for self-absorption

    """
    std_coll = cube.std(axis=axis)       # standard deviation of each spectrum
#    mad_coll = MAD(cube,axis=axis)
    mean_std = median(std_coll.ravel())  # median standard deviation (to reject high-signal spectra that have high std)
    if axis > 0:                         # force spectral axis to first axis
        cube = cube.swapaxes(0,axis)
    width_arr = zeros(cube.shape[1:])   # define gaussian param arrays
    width_arr1 = zeros(cube.shape[1:])   # define gaussian param arrays
    width_arr2 = zeros(cube.shape[1:])   # define gaussian param arrays
    amp_arr = zeros(cube.shape[1:])     # define gaussian param arrays
    amp_arr1 = zeros(cube.shape[1:])     # define gaussian param arrays
    amp_arr2 = zeros(cube.shape[1:])     # define gaussian param arrays
    chi2_arr = zeros(cube.shape[1:])     # define gaussian param arrays
    resid_arr = zeros(cube.shape[1:])    # define gaussian param arrays
    offset_arr = zeros(cube.shape[1:])  # define gaussian param arrays 
    offset_arr1 = zeros(cube.shape[1:])  # define gaussian param arrays 
    offset_arr2 = zeros(cube.shape[1:])  # define gaussian param arrays 
    ncarr = (cube.max(axis=0) > mean_std*nsig)        # cutoff: don't fit no-signal spectra
    starttime = time.time()              # timing for output
    print cube.shape
    print "Fitting a total of %i spectra with peak signal above %f" % (ncarr.sum(),mean_std*nsig)
    for i in xrange(cube.shape[1]):      # Loop over all elements for 
        t0 = time.time()
        nspec = (cube[:,i,:].max(axis=0) > mean_std*nsig).sum()
        print "Working on row %d with %d spectra to fit" % (i,nspec) ,
        for j in xrange(cube.shape[2]):
            if cube[:,i,j].max() > mean_std*nsig:
#            if cube[:,i,j].max() > MAD(cube[:,i,j]):            
                pars = return_param(cube[:,i,j])
                width_arr[i,j] = pars[1]
                width_arr1[i,j] = pars[1]
                amp_arr[i,j] = pars[2]
                amp_arr1[i,j] = pars[2]
#                chi2_arr[i,j] = sum(( gerr(cube[:,i,j])(pars) )**2) 
                resid_arr[i,j] = (gerr(cube[:,i,j])(pars)).sum()
                offset_arr[i,j] = pars[0]
                offset_arr1[i,j] = pars[0]
            else: 
                width_arr1[i,j] = numpy.nan
                chi2_arr[i,j] = numpy.nan
                resid_arr[i,j] = numpy.nan
                offset_arr1[i,j] = numpy.nan
        dt = time.time()-t0
        if nspec > 0:
            print "in %f seconds (average: %f)" % (dt,dt/float(nspec))
        else: 
            print 
    chi2_arr = resid_arr**2
    resids = ma.masked_where(numpy.isnan(chi2_arr),chi2_arr) # hide bad values
#    residcut = (resids.mean() + (resids.std() * nrsig) )  # Old versino - used standard deviation and mean
    residcut = (nanmedian(chi2_arr.ravel()) + (MAD(chi2_arr.ravel()) * nrsig) ) # New version: set cutoff by median + nrsig * MAD
    to_refit = (resids > residcut).astype('bool')
#    to_refit[numpy.isnan(to_refit)] = 0
    inds = array(nonzero(to_refit)).transpose()
    dgc,tgc = 0,0
    print "Refitting a total of %i spectra with peak residual above %f" % (to_refit.sum(),residcut)
    f=open("%s_triples.txt" % prefix,'w')
#    vconv = lambda x: (x-p3+1)*dv+v0    # convert to velocity frame
    vind = vconv(arange(cube[:,0,0].shape[0]))
    xind = arange(cube[:,0,0].shape[0])
    for ind in inds:
        i,j = ind
        doublepars = return_double_param(cube[:,i,j])
        old_chi2 = chi2_arr[i,j]
        new_chi2 = sum(square( double_gerr(cube[:,i,j])(doublepars) )) 
        if new_chi2 < old_chi2: # if 2 gaussians is an improvement, use it!
            chi2_arr[i,j] = new_chi2
            width_arr1[i,j] = doublepars[2]
            width_arr2[i,j] = doublepars[3]
            amp_arr1[i,j] = doublepars[4]
            amp_arr2[i,j] = doublepars[5]
            offset_arr1[i,j] = doublepars[0]
            offset_arr2[i,j] = doublepars[1]
            ncarr[i,j] += 1
        if new_chi2 > residcut: # Even if double was better, see if a triple might be better yet [but don't store it in the params arrays!]
            print >>f,"Triple-gaussian fitting at %i,%i (%i'th double, %i'th triple)" % (i,j,dgc,tgc)
            if tgc % 100 == 0:
                print "Triple-gaussian fitting at %i,%i (%i'th double, %i'th triple)" % (i,j,dgc,tgc)
            tgc += 1
            tpguess = [doublepars[0],doublepars[1],(doublepars[0]+doublepars[1])/2.,doublepars[2],doublepars[3],doublepars[2]*5.,doublepars[4],doublepars[5],doublepars[4]/5.]
            triplepars = return_triple_param(cube[:,i,j],params=tpguess)
            pars = [offset_arr[i,j],width_arr[i,j],amp_arr[i,j]]
            if doplot: # if you don't, there's really no point in fitting at all...
                ax = axes([.05,.05,.7,.9])
                plot(vind,cube[:,i,j],color='black',linestyle='steps',linewidth='.5')
                plot(vind,gaussian(*pars)(xind),'r-.',label="Single %f" % ( (gerr(cube[:,i,j])(pars)).sum() ) )
                plot(vind,double_gaussian(*doublepars)(xind),'g--',label="Double %f" % ( (double_gerr(cube[:,i,j])(doublepars)).sum() ))
                plot(vind,triple_gaussian(*triplepars)(xind),'b:',label="Triple %f" % ( (triple_gerr(cube[:,i,j])(triplepars)).sum() ),linewidth=2)
                pars[0] = vconv(pars[0])
                text(1.05,.8,"c1 %3.2f w1 %3.2f a1 %3.2f" % tuple(pars),transform=ax.transAxes,size='smaller')
                dp = [ vconv(doublepars[0]) , doublepars[2], doublepars[4], vconv(doublepars[1]), doublepars[3], doublepars[5] ]
                text(1.05,.6,"c1 %3.2f w1 %3.2f a1 %3.2f\nc2 %3.2f w2 %3.2f a2 %3.2f" % tuple(dp),transform=ax.transAxes,size='smaller')
                tp = [ vconv(triplepars[0]) , triplepars[3], triplepars[6], vconv(triplepars[1]), triplepars[4], triplepars[7], vconv(triplepars[2]), triplepars[5], triplepars[8]  ]
                text(1.05,.4,"c1 %3.2f w1 %3.2f a1 %3.2f\nc2 %3.2f w2 %3.2f a2 %3.2f\nc3 %3.2f w3 %3.2f a3 %3.2f" % tuple(tp),transform=ax.transAxes,size='smaller')
                title("Spectrum at %s %s" % (ratos(xtora(i)),dectos(ytodec(j))) ) 
                legend(loc='best')
                savefig("%s_%s.%s.png" % (prefix,i,j))
                clf()
            ncarr[i,j] += 1
            print >>f,triplepars
        dgc += 1

    f.close()
    print "Total time %f seconds for %i double and %i triple gaussians" % (time.time()-starttime,dgc,tgc)

    return width_arr1,width_arr2,chi2_arr,offset_arr1,offset_arr2,amp_arr1,amp_arr2,ncarr

def collapse_gaussfit(cube,axis=2):
    std_coll = cube.std(axis=axis)
    mean_std = median(std_coll.ravel())
    if axis > 0:
        cube = cube.swapaxes(0,axis)
    width_arr = zeros(cube.shape[1:])
    amp_arr = zeros(cube.shape[1:])
    chi2_arr = zeros(cube.shape[1:])
    offset_arr = zeros(cube.shape[1:])
    starttime = time.time()
    print cube.shape
    print "Fitting a total of %i spectra with peak signal above %f" % ((cube.max(axis=0) > mean_std).sum(),mean_std)
    for i in xrange(cube.shape[1]):
        t0 = time.time()
        nspec = (cube[:,i,:].max(axis=0) > mean_std).sum()
        print "Working on row %d with %d spectra to fit" % (i,nspec) ,
        for j in xrange(cube.shape[2]):
            if cube[:,i,j].max() > mean_std:
                pars = return_param(cube[:,i,j])
                width_arr[i,j] = pars[1]
                chi2_arr[i,j] = sum(( gerr(cube[:,i,j])(pars) )**2) 
                offset_arr[i,j] = pars[0]
                amp_arr[i,j] = pars[2]
            else: 
                width_arr[i,j] = numpy.nan
                chi2_arr[i,j] = numpy.nan
                offset_arr[i,j] = numpy.nan
                amp_arr[i,j] = numpy.nan
        dt = time.time()-t0
        if nspec > 0:
            print "in %f seconds (average: %f)" % (dt,dt/float(nspec))
    print "Total time %f seconds" % (time.time()-starttime)

    return width_arr,offset_arr,amp_arr,chi2_arr

# next step: find 2-gaussian fits
def collapse_double_gaussfit(cube,axis=2):
    std_coll = cube.std(axis=axis)
    mean_std = median(std_coll.ravel())
    if axis > 0:
        cube = cube.swapaxes(0,axis)
    width_arr1 = zeros(cube.shape[1:])
    width_arr2 = zeros(cube.shape[1:])
    amp_arr1 = zeros(cube.shape[1:])
    amp_arr2 = zeros(cube.shape[1:])
    chi2_arr = zeros(cube.shape[1:])
    offset_arr1 = zeros(cube.shape[1:])
    offset_arr2 = zeros(cube.shape[1:])
    starttime = time.time()
    print cube.shape
    print "Fitting a total of %i spectra with peak signal above %f" % ((cube.max(axis=0) > mean_std).sum(),mean_std)
    for i in xrange(cube.shape[1]):
        t0 = time.time()
        nspec = (cube[:,i,:].max(axis=0) > mean_std).sum()
        print "Working on row %d with %d spectra to fit" % (i,nspec) ,
        for j in xrange(cube.shape[2]):
            if cube[:,i,j].max() > mean_std:
                pars = return_double_param(cube[:,i,j])
                width_arr1[i,j] = pars[2]
                width_arr2[i,j] = pars[3]
                amp_arr1[i,j] = pars[4]
                amp_arr2[i,j] = pars[5]
                chi2_arr[i,j] = sum(( double_gerr(cube[:,i,j])(pars) )**2) 
                offset_arr1[i,j] = pars[0]
                offset_arr2[i,j] = pars[1]
            else: 
                width_arr1[i,j] = numpy.nan
                width_arr2[i,j] = numpy.nan
                chi2_arr[i,j] = numpy.nan
                offset_arr1[i,j] = numpy.nan
                offset_arr2[i,j] = numpy.nan
        dt = time.time()-t0
        if nspec > 0:
            print "in %f seconds (average: %f)" % (dt,dt/float(nspec))
    print "Total time %f seconds" % (time.time()-starttime)

    return width_arr1,width_arr2,chi2_arr,offset_arr1,offset_arr2,amp_arr1,amp_arr2

def wrap_collapse_gauss(filename,outprefix,redo='no'):
    """

    redo - if not equal to 'no', then...
    if collapse_gaussfit succeeded (to the extent that the .pysav files were written),
    but some part of the file writing or successive procedures failed, re-do those 
    procedures without redoing the whole collapse
    """
    fitsfile = pyfits.open(filename)
    dv,v0,p3 = fitsfile[0].header['CD3_3'],fitsfile[0].header['CRVAL3'],fitsfile[0].header['CRPIX3']

    cube = fitsfile[0].data
    cube = where(numpy.isnan(cube),0,cube)

    if redo=='no':
        doubleB = asarray(collapse_double_gaussfit(cube,axis=0))
        doubleB[numpy.isnan(doubleB)] = 0
        pickle.dump(doubleB,open('%s_doubleB.pysav' % outprefix,'w'))
    else:
        doubleB = pickle.load(open('%s_doubleB.pysav' % outprefix,'r'))
    db = doubleB
    gcd = double_gaussian(db[3],db[4],db[0],db[1],db[5],db[6])(indices(cube.shape)[0])
    fitsfile[0].data = gcd
    fitsfile.writeto('%s_doublegausscube.fits' % outprefix,clobber=True)

    gcd[numpy.isnan(gcd)] = 0
    doubleResids = cube-gcd
    fitsfile[0].data = doubleResids
    fitsfile.writeto('%s_doublegaussresids.fits' % outprefix,clobber=True)


    #doubleB[4] = (doubleB[4]-v0) / dv + p3-1
    #doubleB[3] = (doubleB[3]-v0) / dv + p3-1
    doubleB[4] = (doubleB[4]-p3+1) * dv + v0
    doubleB[3] = (doubleB[3]-p3+1) * dv + v0
    fitsfile[0].data = asarray(doubleB)
    fitsfile.writeto('%s_doublegausspars.fits' % outprefix,clobber=True)

    if redo=='no':
        singleB = asarray(collapse_gaussfit(cube,axis=0))
        pickle.dump(singleB,open('%s_singleB.pysav' % outprefix,'w'))
    else:
        singleB = pickle.load(open('%s_singleB.pysav' % outprefix,'r'))
    gc = gaussian(singleB[1],singleB[0],singleB[2])(indices(cube.shape)[0])
    singleB[1] = (singleB[1]-p3+1) * dv + v0
    fitsfile[0].data = gc
    fitsfile.writeto('%s_singlegausscube.fits' % outprefix,clobber=True)

    gc[numpy.isnan(gc)]=0
    singleResids = cube-gc
    fitsfile[0].data = singleResids
    fitsfile.writeto('%s_singlegaussresids.fits' % outprefix,clobber=True)
    fitsfile[0].data = asarray(singleB)
    fitsfile.writeto('%s_singlegausspars.fits' % outprefix,clobber=True)

    fitsfile[0].header.__delitem__('CD3_3')
    fitsfile[0].header.__delitem__('CRVAL3')
    fitsfile[0].header.__delitem__('CRPIX3')
    fitsfile[0].header.__delitem__('CUNIT3')
    fitsfile[0].header.__delitem__('CTYPE3')

    doubleResids[numpy.isnan(doubleResids)] = 0
    totalDResids = doubleResids.sum(axis=0)
    fitsfile[0].data = totalDResids
    fitsfile.writeto('%s_doublegauss_totalresids.fits' % outprefix,clobber=True)

    singleResids[numpy.isnan(singleResids)] = 0
    totalSResids = singleResids.sum(axis=0)
    fitsfile[0].data = totalSResids
    fitsfile.writeto('%s_singlegauss_totalresids.fits' % outprefix,clobber=True)

    return singleB,doubleB


def wrap_collapse_adaptive(filename,outprefix,redo='no',nsig=5,nrsig=2,doplot=True):
    """
    redo - if not equal to 'no', then...
    if collapse_gaussfit succeeded (to the extent that the .pysav files were written),
    but some part of the file writing or successive procedures failed, re-do those 
    procedures without redoing the whole collapse
    """
    fitsfile = pyfits.open(filename)
    dv,v0,p3 = fitsfile[0].header['CD3_3'],fitsfile[0].header['CRVAL3'],fitsfile[0].header['CRPIX3']
    dr,r0,p1 = fitsfile[0].header['CD1_1'],fitsfile[0].header['CRVAL1'],fitsfile[0].header['CRPIX1']
    dd,d0,p2 = fitsfile[0].header['CD2_2'],fitsfile[0].header['CRVAL2'],fitsfile[0].header['CRPIX2']
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0

    cube = fitsfile[0].data
    cube = where(numpy.isnan(cube),0,cube)

    if redo=='no':
        adaptB = asarray(adaptive_collapse_gaussfit(cube,axis=0,prefix=outprefix+'_triple',
            nsig=nsig,nrsig=nrsig,vconv=vconv,xtora=xtora,ytodec=ytodec,doplot=doplot))
        adaptB[numpy.isnan(adaptB)] = 0
        pickle.dump(adaptB,open('%s_adaptB.pysav' % outprefix,'w'))
    else:
        adaptB = pickle.load(open('%s_adaptB.pysav' % outprefix,'r'))
    db = adaptB
    gcd = double_gaussian(db[3],db[4],db[0],db[1],db[5],db[6])(indices(cube.shape)[0])
    fitsfile[0].data = gcd
    fitsfile.writeto('%s_adaptgausscube.fits' % outprefix,clobber=True)

    gcd[numpy.isnan(gcd)] = 0
    adaptResids = cube-gcd
    fitsfile[0].data = adaptResids
    fitsfile.writeto('%s_adaptgaussresids.fits' % outprefix,clobber=True)


    #adaptB[4] = (adaptB[4]-v0) / dv + p3-1
    #adaptB[3] = (adaptB[3]-v0) / dv + p3-1
    adaptB[4] = (adaptB[4]-p3+1) * dv + v0
    adaptB[3] = (adaptB[3]-p3+1) * dv + v0
    fitsfile[0].data = asarray(adaptB)
    fitsfile.writeto('%s_adaptgausspars.fits' % outprefix,clobber=True)

    fitsfile[0].header.__delitem__('CD3_3')
    fitsfile[0].header.__delitem__('CRVAL3')
    fitsfile[0].header.__delitem__('CRPIX3')
    fitsfile[0].header.__delitem__('CUNIT3')
    fitsfile[0].header.__delitem__('CTYPE3')

    adaptResids[numpy.isnan(adaptResids)] = 0
    totalDResids = adaptResids.sum(axis=0)
    fitsfile[0].data = totalDResids
    fitsfile.writeto('%s_adaptgauss_totalresids.fits' % outprefix,clobber=True)

    return adaptB

