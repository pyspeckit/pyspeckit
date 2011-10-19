"""
======================================
Generalized hyperfine component fitter
======================================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import numpy as np
import model,fitter

ckms = 2.99792458e5

class hyperfinemodel(object):
    """
    Wrapper for the hyperfine model class.  Specify the offsets and relative
    strengths when initializing, then you've got yourself a hyperfine modeler.
    """

    def __init__(self, line_names, voff_lines_dict, freq_dict, line_strength_dict):
        """
        Initialize the various parameters defining the hyperfine transitions

        line_names is a LIST of the line names to be used as indices for the dictionaries

        voff_lines_dict is a linename:v_off dictionary of velocity offsets for the hyperfine components.  Technically,
            this is redundant with freq_dict
        freq_dict - frequencies of the indvidual transitions

        line_strength_dict - Relative strengths of the hyperfine components, usually determined by their degeneracy and 
            Einstein A coefficients
        """
        self.line_names = line_names
        self.voff_lines_dict = voff_lines_dict
        self.freq_dict = freq_dict
        self.line_strength_dict = line_strength_dict

        self.fitter = model.SpectralModel(self,4,
            parnames=['Tex','tau','center','width'], 
            parlimited=[(False,False), (True,False), (False,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0)],
            shortvarnames=("T_{ex}","\\tau","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )

        self.vheight_fitter = model.SpectralModel(fitter.vheightmodel(self),5,
            parnames=['height','Tex','tau','center','width'], 
            parlimited=[(False,False), (False,False), (True,False), (False,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=("H","T_{ex}","\\tau","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )

        self.ampfitter = model.SpectralModel(self.hyperfine_amp,3,
            parnames=['amp','center','width'], 
            parlimited=[(False,False), (False,False), (True,False)], 
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=("amp","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunits='Hz' )

    def __call__(self, *args, **kwargs):
        """
        Generate a model spectrum given an excitation temperature, optical depth, offset velocity, and velocity width.
        """
        return self.hyperfine(*args,**kwargs)

    def hyperfine_amp(self, xarr, amp=None, xoff_v=0.0, width=1.0, 
            return_components=False, Tbackground=2.73, Tex=5.0, tau=0.1):
        """
        wrapper of self.hyperfine with order of arguments changed
        """ 
        return self.hyperfine(xarr, amp=amp, Tex=Tex, tau=tau, xoff_v=xoff_v,
                width=width, return_components=return_components,
                Tbackground=Tbackground)


    def hyperfine(self, xarr, Tex=5.0, tau=0.1, xoff_v=0.0, width=1.0, 
            return_components=False, Tbackground=2.73, amp=None ):
        """
        Generate a model spectrum given an excitation temperature, optical depth, offset velocity, and velocity width.
        """

        # Convert X-units to frequency in Hz
        xarr = xarr.as_unit('Hz')

        # Generate an optical depth spectrum as a function of the X-axis 
        tau_nu_cumul = np.zeros(len(xarr))
        # Error check: inputing NANs results in meaningless output - return without computing a model
        if np.any(np.isnan((tau,Tex,width,xoff_v))):
            if return_components:
                return [tau_nu_cumul] * len(self.line_names)
            else:
                return tau_nu_cumul

        components =[]
        for linename in self.line_names:
            voff_lines = np.array(self.voff_lines_dict[linename])
      
            lines = (1-voff_lines/ckms)*self.freq_dict[linename]
            if width == 0:
                tau_nu = xarr*0
            else:
                nuwidth = np.abs(width/ckms*lines)
                nuoff = xoff_v/ckms*lines
                # the total optical depth, which is being fitted, should be the sum of the components
                tau_line = (tau * self.line_strength_dict[linename]/np.sum(self.line_strength_dict.values())) 
          
                tau_nu = np.array(tau_line * np.exp(-(xarr+nuoff-self.freq_dict[linename])**2/(2.0*nuwidth**2)))
                tau_nu[tau_nu!=tau_nu] = 0 # avoid nans
                components.append( tau_nu )
            tau_nu_cumul += tau_nu

        # add a list of the individual 'component' spectra to the total components...

        if return_components:
            if amp is None:
                return (1.0-np.exp(-np.array(components)))*(Tex-Tbackground)
            else:
                comps = (1.0-np.exp(-np.array(components)))*(Tex-Tbackground)
                return comps/comps.max() * amp

        spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*(Tex-Tbackground)
      
        if amp is None:
            return spec
        else:
            return spec/spec.max() * amp



