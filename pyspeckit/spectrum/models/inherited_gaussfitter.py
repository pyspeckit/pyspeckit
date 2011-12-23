import model
import numpy 

def gaussian(x,A,dx,w):
    """
    Returns a 1-dimensional gaussian of form
    H+A*numpy.exp(-(x-dx)**2/(2*w**2))
    
    [height,amplitude,center,width]
    
    """
    x = numpy.array(x) # make sure xarr is no longer a spectroscopic axis
    return A*numpy.exp(-(x-dx)**2/(2*w**2))

def gaussian_fitter(multisingle='multi'):

    return model.SpectralModel(gaussian, 3,
            parnames=['amplitude','shift','width'], 
            parlimited=[(False,False),(False,False),(True,False)], 
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=('A',r'\Delta x',r'\sigma'),
            multisingle=multisingle,
            )
