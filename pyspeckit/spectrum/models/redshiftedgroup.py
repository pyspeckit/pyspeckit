"""
========================================
Generalized redshifted line group fitter
========================================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import numpy as np
import model,fitter

ckms = 2.99792458e5

class redshiftedgroup(object):
    """
    Wrapper for the redshifted group model class.  Specify the offsets and rest
    frequencies when initializing.
    Widths will be matched.
    """

    def __init__(self, freq_dict):
        """
        Initialize the various parameters defining the redshiftedgroup transitions

        freq_dict - frequencies of the indvidual transitions indexed by name
        """
        self.freq_dict = freq_dict

        self.fitter = model.SpectralModel(self,3,
            parnames=['amp','center','width'], 
            parlimited=[(False,False), (False,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=("A","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )

        self.vheight_fitter = model.SpectralModel(fitter.vheightmodel(self),4,
            parnames=['height','amp','center','width'], 
            parlimited=[(False,False), (False,False), (True,False), (False,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0)],
            shortvarnames=("H","A","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )


    def __call__(self, *args, **kwargs):
        """
        Generate a model spectrum given an excitation temperature, optical depth, offset velocity, and velocity width.
        """
        return self.redshiftedgroup(*args,**kwargs)


    def redshiftedgroup(self, xarr, amp=1.0, redshift=0.0, width=1.0, 
            return_components=False ):
        """
        Generate a model spectrum given an amplitude, offset velocity, and velocity width.
        """

        # Convert X-units to frequency in Hz
        xarr = xarr.as_unit('Hz')

        # Generate an optical depth spectrum as a function of the X-axis 
        amp_cumul = np.zeros(len(xarr))
        # Error check: inputing NANs results in meaningless output - return without computing a model
        if np.any(np.isnan((amp,width,redshift))):
            if return_components:
                return [amp_cumul] * len(self.freq_dict)
            else:
                return amp_cumul

        components =[]
        for linename in self.freq_dict:
      
            redline = self.freq_dict[linename] / (1+redshift)
            if width == 0:
                amparr = xarr*0
            else:
                amparr = np.array(amp * np.exp(-(xarr-redline)**2/(2.0*width**2)))
                amparr[amparr!=amparr] = 0 # avoid nans
                components.append( amparr )
            amp_cumul += amparr

        # add a list of the individual 'component' spectra to the total components...

        if return_components:
            return np.array(components)
        else:
            return np.array(amp_cumul)
