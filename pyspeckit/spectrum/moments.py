import numpy

def moments(Xax, data, vheight=True, estimator=numpy.median, negamp=None,
        veryverbose=False,  **kwargs):
    """Returns (height, amplitude, x, width_x)
    the gaussian parameters of a 1D distribution by calculating its
    moments.  Depending on the input parameters, will only output 
    a subset of the above.
    "height" is the background level
    "amplitude" is the maximum (or minimum) of the data after background subtraction
    "x" is the first moment
    "width_x" is the second moment
    
    If using masked arrays, pass estimator=numpy.ma.median
    'estimator' is used to measure the background level (height)

    negamp can be used to force the peak negative (True), positive (False),
    or it will be "autodetected" (negamp=None)
    """

    Xax = numpy.array(Xax)

    dx = numpy.mean(numpy.diff(Xax)) # assume a regular grid
    integral = (data*dx).sum()
    height = estimator(data)
    
    # try to figure out whether pos or neg based on the minimum width of the pos/neg peaks
    Lpeakintegral = integral - height*len(Xax)*dx - (data[data>height]*dx).sum()
    Lamplitude = data.min()-height
    Lwidth_x = 0.5*(numpy.abs(Lpeakintegral / Lamplitude))
    Hpeakintegral = integral - height*len(Xax)*dx - (data[data<height]*dx).sum()
    Hamplitude = data.max()-height
    Hwidth_x = 0.5*(numpy.abs(Hpeakintegral / Hamplitude))
    Lstddev = Xax[data<data.mean()].std()
    Hstddev = Xax[data>data.mean()].std()
    #print "Lstddev: %10.3g  Hstddev: %10.3g" % (Lstddev,Hstddev)
    #print "Lwidth_x: %10.3g  Hwidth_x: %10.3g" % (Lwidth_x,Hwidth_x)

    if negamp: # can force the guess to be negative
        xcen,amplitude,width_x = Xax[numpy.argmin(data)],Lamplitude,Lwidth_x
    elif negamp is None:
        if Hstddev < Lstddev: 
            xcen,amplitude,width_x, = Xax[numpy.argmax(data)],Hamplitude,Hwidth_x
        else:                                                                   
            xcen,amplitude,width_x, = Xax[numpy.argmin(data)],Lamplitude,Lwidth_x
    else:  # if negamp==False, make positive
        xcen,amplitude,width_x = Xax[numpy.argmax(data)],Hamplitude,Hwidth_x

    if veryverbose:
        print "negamp: %s  amp,width,cen Lower: %g, %g   Upper: %g, %g  Center: %g" %\
                (negamp,Lamplitude,Lwidth_x,Hamplitude,Hwidth_x,xcen)
    mylist = [amplitude,xcen,width_x]
    if negamp and amplitude > 0 and veryverbose: 
        print "WARNING: likely fit failure.  negamp=True, but amplitude > 0"
    if negamp is False and amplitude < 0 and veryverbose: 
        print "WARNING: likely fit failure.  negamp=False, but amplitude < 0"
    if numpy.isnan(width_x) or numpy.isnan(height) or numpy.isnan(amplitude):
        raise ValueError("something is nan")
    if vheight:
        mylist = [height] + mylist
    return mylist
