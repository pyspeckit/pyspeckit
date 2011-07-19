import numpy as np

def smooth(data,smooth,smoothtype='gaussian',downsample=True,downsample_factor=None,
        convmode='same'):
    """
    Smooth and downsample the data array
    smoothtype - 'gaussian','hanning', or 'boxcar'
    downsample - Downsample the data?
    downsample_factor - Downsample by the smoothing factor, or something else?
    convmode - see numpy.convolve.  'same' returns an array of the same length as
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

    smdata = np.convolve(data,kernel,convmode)[::downsample_factor]

    return smdata
