import numpy as np

def moments(Xax, data, vheight=True, estimator=np.median, negamp=None,
        veryverbose=False, nsigcut=None, noise_guess=None,  **kwargs):
    """Returns (height, amplitude, x, width_x)
    the gaussian parameters of a 1D distribution by calculating its
    moments.  Depending on the input parameters, will only output 
    a subset of the above.
    "height" is the background level
    "amplitude" is the maximum (or minimum) of the data after background subtraction
    "x" is the first moment
    "width_x" is the second moment
    
    If using masked arrays, pass estimator=np.ma.median
    'estimator' is used to measure the background level (height)

    negamp can be used to force the peak negative (True), positive (False),
    or it will be "autodetected" (negamp=None)

    "nsigcut" - try to estimate the noise and only use data above/below nsigma
        above the noise.  Estimate the noise from the data unless passed as a
        keyword

    Theory:
    From first principles (in the absence of noise):
    integral(gaussian) = sqrt(2*pi*sigma^2) * amp
    sigma = integral / amp / sqrt(2*pi)
    
    in the presence of noise, this gets much more complicated

    """

    Xax = np.array(Xax)

    dx = np.abs(np.mean(np.diff(Xax))) # assume a regular grid
    integral = (data*dx).sum()
    height = estimator(data)

    if noise_guess is None:
        noise_guess = data.std() 
    height_cut_low = height-nsigcut*noise_guess if nsigcut is not None else height
    height_cut_high = height+nsigcut*noise_guess if nsigcut is not None else height
    
    # try to figure out whether pos or neg based on the minimum width of the pos/neg peaks
    Lpeakintegral = integral - height_cut_low*(data<height_cut_low).sum()*dx - (data[data>height_cut_low]*dx).sum()
    Lamplitude = data.min()-height
    Lwidth_x = Lpeakintegral / Lamplitude / np.sqrt(2*np.pi)

    Hpeakintegral = integral - height*(data>height_cut_high).sum()*dx - (data[data<height_cut_high]*dx).sum()
    Hamplitude = data.max()-height
    Hwidth_x = Hpeakintegral / Hamplitude / np.sqrt(2*np.pi)

    Lstddev = Xax[data<data.mean()].std()
    Hstddev = Xax[data>data.mean()].std()
    #print "Lstddev: %10.3g  Hstddev: %10.3g" % (Lstddev,Hstddev)
    #print "Lwidth_x: %10.3g  Hwidth_x: %10.3g" % (Lwidth_x,Hwidth_x)

    if negamp: # can force the guess to be negative
        xcen,amplitude,width_x = Xax[np.argmin(data)],Lamplitude,Lwidth_x
    elif negamp is None:
        if Hstddev < Lstddev: 
            xcen,amplitude,width_x, = Xax[np.argmax(data)],Hamplitude,Hwidth_x
        else:                                                                   
            xcen,amplitude,width_x, = Xax[np.argmin(data)],Lamplitude,Lwidth_x
    else:  # if negamp==False, make positive
        xcen,amplitude,width_x = Xax[np.argmax(data)],Hamplitude,Hwidth_x

    if veryverbose:
        print "negamp: %s  amp,width,cen Lower: %g, %g   Upper: %g, %g  Center: %g" %\
                (negamp,Lamplitude,Lwidth_x,Hamplitude,Hwidth_x,xcen)
    mylist = [amplitude,xcen,width_x]
    if negamp and amplitude > 0 and veryverbose: 
        print "WARNING: likely fit failure.  negamp=True, but amplitude > 0"
    if negamp is False and amplitude < 0 and veryverbose: 
        print "WARNING: likely fit failure.  negamp=False, but amplitude < 0"
    if np.isnan(width_x) or np.isnan(height) or np.isnan(amplitude):
        raise ValueError("something is nan")
    if vheight:
        mylist = [height] + mylist
    return mylist
