import numpy as np
import arithmetic

def smooth(data,smooth,smoothtype='gaussian',downsample=True,downsample_factor=None,
        convmode='same'):
    """
    Smooth and downsample the data array

    *smooth* [ float ] 
        Number of pixels to smooth by

    *smoothtype* [ 'gaussian','hanning', or 'boxcar' ]
        type of smoothing kernel to use

    *downsample* [ bool ]
        Downsample the data?

    *downsample_factor* [ int ]
        Downsample by the smoothing factor, or something else?

    *convmode* [ 'full','valid','same' ]
        see `numpy.convolve`.  'same' returns an array of the same length as
        'data' (assuming data is larger than the kernel)
    """
    
    roundsmooth = round(smooth) # can only downsample by integers

    if downsample_factor is None and downsample:
        downsample_factor = int(roundsmooth)
    elif downsample_factor is None:
        downsample_factor = 1

    if smooth > len(data) or downsample_factor > len(data):
        raise ValueError("Error: trying to smooth by more than the spectral length.")

    if smoothtype == 'hanning':
        kernel = np.hanning(2+roundsmooth)/np.hanning(2+roundsmooth).sum()
    elif smoothtype == 'gaussian':
        xkern  = np.linspace(-5*smooth,5*smooth,smooth*11)
        kernel = np.exp(-xkern**2/(2*(smooth/np.sqrt(8*np.log(2)))**2))
        kernel /= kernel.sum()
        if len(kernel) > len(data):
            lengthdiff = len(kernel)-len(data)
            if lengthdiff % 2 == 0: # make kernel same size as data
                kernel = kernel[lengthdiff/2:-lengthdiff/2]
            else: # make kernel 1 pixel smaller than data but still symmetric
                kernel = kernel[lengthdiff/2+1:-lengthdiff/2-1]
    elif smoothtype == 'boxcar':
        kernel = np.ones(roundsmooth)/float(roundsmooth)

    # deal with NANs or masked values
    if hasattr(data,'mask'):
        if type(data.mask) is np.ndarray:
            OK = True - data.mask
            data = arithmetic._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])
    if np.any(True - np.isfinite(data)):
        OK = np.isfinite(data)
        data = arithmetic._interp(np.arange(len(data)),np.arange(len(data))[OK],data[OK])

    if np.any(True - np.isfinite(data)):
        raise ValueError("NANs in data array after they have been forcibly removed.")

    smdata = np.convolve(data,kernel,convmode)[::downsample_factor]

    return smdata

def smooth_multispec(data,smoothfactor,**kwargs):
    """
    Smooth multiple spectra as from an ObsBlock (shape should be [speclen, nspec])
    """

    nobs = data.shape[1]

    newdata = np.array( [smooth(D,smoothfactor,**kwargs) for D in data.swapaxes(0,1)]).swapaxes(0,1)

    return newdata
