"""
======================================
Generalized hyperfine component fitter
======================================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^
"""
import numpy as np
from astropy import units as u
from astropy import constants
import copy

from . import model
from . import fitter

# should be imported in the future
ckms = 2.99792458e5
hoverk = (constants.h.cgs/constants.k_B.cgs).value

class hyperfinemodel(object):
    """
    Wrapper for the hyperfine model class.  Specify the offsets and relative
    strengths when initializing, then you've got yourself a hyperfine modeler.

    There are a wide variety of different fitter attributes, each designed to
    free a different subset of the parameters.  Their purposes should be
    evident from their names.
    """

    def __init__(self, line_names, voff_lines_dict, freq_dict,
                 line_strength_dict, relative_strength_total_degeneracy):
        """
        Initialize the various parameters defining the hyperfine transitions

        Parameters
        ----------
        line_names: list
            list of the line names to be used as indices for the dictionaries
        voff_lines_dict: dict
            a linename:v_off dictionary of velocity offsets for the hyperfine
            components.  Technically, this is redundant with freq_dict
        freq_dict: dict
            frequencies of the indvidual transitions
        line_strength_dict: dict
            Relative strengths of the hyperfine components, usually determined
            by their degeneracy and Einstein A coefficients
        """
        self.line_names = tuple(line_names)
        self.voff_lines_dict = voff_lines_dict
        self.freq_dict = freq_dict
        self.line_strength_dict = line_strength_dict
        self.relative_strength_total_degeneracy = relative_strength_total_degeneracy

        self.fitter = model.SpectralModel(self,4,
            parnames=['Tex','tau','center','width'],
            # T_ex = 0 results in an infinity
            parlimited=[(True,False), (True,False), (False,False), (True,False)],
            parlimits=[(1e-5,0), (0,0), (0,0), (0,0)],
            # specify the parameter names (LaTeX is OK)
            shortvarnames=("T_{ex}","\\tau","v","\\sigma"),
            guess_types=['amplitude+2.73', 1.0, 'center', 'width'],
            fitunit='Hz')

        self.nlines = len(line_names)

        self.varyhf_fitter = model.SpectralModel(self.hyperfine_varyhf,3+self.nlines,
            parnames=['Tex','center','width']+['tau%s' % k for k in self.line_names],
            parlimited=[(True,False), (False,False), (True,False)]
                                                 + [(True,False),]*self.nlines,
            parlimits=[(1e-5,0), (0,0), (0,0)]+[(0,0),]*self.nlines,
            shortvarnames=("T_{ex}","v","\\sigma") +
                           tuple(("\\tau(\\mathrm{%s})" % k for k in self.line_names)),
            fitunit='Hz')

        self.varyhf_amp_fitter = model.SpectralModel(self.hyperfine_varyhf_amp, 2+self.nlines,
            parnames=['center','width']+['amp%s' % k for k in self.line_names],
            parlimited=[(False,False), (True,False)] + [(True,False),]*self.nlines,
            parlimits=[(0,0), (0,0)]+[(0,0),]*self.nlines,
            shortvarnames=("v","\\sigma") +
                           tuple(("amp(\\mathrm{%s})" % k for k in self.line_names)),
            fitunit='Hz')

        self.varyhf_amp_width_fitter = model.SpectralModel(self.hyperfine_varyhf_amp_width,1+self.nlines*2,
            parnames=['center']+['amp%s' % k for k in self.line_names]+['width%s' % k for k in self.line_names],
            parlimited=[(False,False)] + [(True,False),]*self.nlines + [(True,False)]*self.nlines,
            parlimits=[(0,0)]+[(0,0),]*self.nlines*2,
            shortvarnames=(("v",) +
                           tuple(("amp(\\mathrm{%s})" % k for k in self.line_names)) +
                           tuple(("\\sigma(\\mathrm{%s})" % k for k in self.line_names))),
            # specify the parameter names (TeX is OK)
            fitunit='Hz' )

        self.vheight_fitter = model.SpectralModel(fitter.vheightmodel(self),5,
            parnames=['height','Tex','tau','center','width'],
            parlimited=[(False,False), (True,False), (True,False), (False,False), (True,False)],
            parlimits=[(0,0), (1e-5,0), (0,0), (0,0), (0,0)],
            shortvarnames=("H","T_{ex}","\\tau","v","\\sigma"), # specify the parameter names (TeX is OK)
            guess_types=[0.0, 'amplitude+2.73', 1.0, 'center', 'width'],
            fitunit='Hz' )

        self.background_fitter = model.SpectralModel(self.hyperfine_addbackground,5,
            parnames=['Tbackground','Tex','tau','center','width'],
            parlimited=[(True,False), (True,False), (False,False), (True,False), (False,False), (True,False)],
            parlimits=[(0,0), (1e-5,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=('T_{BG}',"T_{ex}","\\tau","v","\\sigma"), # specify the parameter names (TeX is OK)
            guess_types=[2.73, 'amplitude+2.73', 1.0, 'center', 'width'],
            fitunit='Hz')

        self.background_contsub_fitter = model.SpectralModel(self.hyperfine_background,5,
            parnames=['Tbackground','Tex','tau','center','width'],
            parlimited=[(True,False), (True,False), (False,False), (True,False), (False,False), (True,False)],
            parlimits=[(0,0), (1e-5,0), (0,0), (0,0), (0,0), (0,0)],
            shortvarnames=('T_{BG}',"T_{ex}","\\tau","v","\\sigma"), # specify the parameter names (TeX is OK)
            guess_types=[0.0, 'amplitude+2.73', 1.0, 'center', 'width'],
            fitunit='Hz')

        self.ampfitter = model.SpectralModel(self.hyperfine_amp,3,
            parnames=['amp','center','width'],
            parlimited=[(False,False), (False,False), (True,False)],
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=("amp","v","\\sigma"), # specify the parameter names (TeX is OK)
            fitunit='Hz' )

        self.taufitter = model.SpectralModel(self.hyperfine_tau,3,
            parnames=['tau','center','width'],
            parlimited=[(True,False), (False,False), (True,False)],
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=(r'\tau',"v","\\sigma"), # specify the parameter names (TeX is OK)
            guess_types=[1.0, 'center', 'width'],
            fitunit='Hz')

        self.totaltaufitter = model.SpectralModel(self.hyperfine_tau_total,3,
            parnames=['tau','center','width'],
            parlimited=[(True,False), (False,False), (True,False)],
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=(r'\tau',"v","\\sigma"), # specify the parameter names (TeX is OK)
            guess_types=[1.0, 'center', 'width'],
            fitunit='Hz')

    def __copy__(self):
        # http://stackoverflow.com/questions/1500718/what-is-the-right-way-to-override-the-copy-deepcopy-operations-on-an-object-in-p
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        # A deep copy of the hyperfine model is OK to just do regular copies of
        # all attributes, since none of them are meant to be modified
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))
        return result

    def __call__(self, *args, **kwargs):
        """
        Generate a model spectrum given an excitation temperature, optical depth, offset velocity, and velocity width.
        """
        return self.hyperfine(*args,**kwargs)

    def hyperfine_amp(self, xarr, amp=None, xoff_v=0.0, width=1.0,
                      return_hyperfine_components=False, Tbackground=2.73,
                      Tex=5.0, tau=0.1):
        """
        wrapper of self.hyperfine with order of arguments changed
        """
        return self.hyperfine(xarr, amp=amp, Tex=Tex, tau=tau, xoff_v=xoff_v,
                width=width, return_hyperfine_components=return_hyperfine_components,
                Tbackground=Tbackground)

    def hyperfine_tau(self, xarr, tau, xoff_v, width, **kwargs):
        """ same as hyperfine, but with arguments in a different order, AND
        tau is returned instead of exp(-tau)"""
        return self.hyperfine(xarr, tau=tau, xoff_v=xoff_v, width=width,
                return_tau=True, **kwargs)

    def hyperfine_tau_total(self, xarr, tau_total, xoff_v, width, **kwargs):
        """ same as hyperfine, but with arguments in a different order, AND
        tau is returned instead of exp(-tau), AND the *peak* tau is used"""
        return self.hyperfine(xarr, tau_total=tau_total, xoff_v=xoff_v, width=width,
                return_tau=True, **kwargs)

    def hyperfine_varyhf(self, xarr, Tex, xoff_v, width, *args, **kwargs):
        """ Wrapper of hyperfine for using a variable number of peaks with specified
        tau """
        return self.hyperfine(xarr, Tex=Tex, xoff_v=xoff_v, width=width,
                              tau=dict(zip(self.line_names,args)),
                              vary_hyperfine_tau=True, **kwargs)

    def hyperfine_varyhf_amp(self, xarr, xoff_v, width, *args, **kwargs):
        """ Wrapper of hyperfine for using a variable number of peaks with specified
        amplitude (rather than tau).  Uses some opaque tricks: Tex is basically ignored,
        and return_tau means you're actually returning the amplitude,
        which is just passed in as tau"""
        return self.hyperfine(xarr, xoff_v=xoff_v, width=width,
                              tau=dict(zip(self.line_names,args)),
                              vary_hyperfine_tau=True,
                              return_tau=True, **kwargs)

    def hyperfine_varyhf_amp_width(self, xarr, xoff_v, *args, **kwargs):
        """ Wrapper of hyperfine for using a variable number of peaks with specified
        amplitude (rather than tau).  Uses some opaque tricks: Tex is basically ignored,
        and return_tau means you're actually returning the amplitude,
        which is just passed in as tau"""
        if len(args) % 2 != 0:
            raise ValueError("Incorrect number of arguments for varying amplitude"
                             " and width.  Need N amplitudes, N widths.")
        nargs = int(len(args)/2)
        return self.hyperfine(xarr, xoff_v=xoff_v,
                              tau=dict(zip(self.line_names,args[:nargs])),
                              width=dict(zip(self.line_names,args[nargs:])),
                              vary_hyperfine_tau=True,
                              vary_hyperfine_width=True,
                              return_tau=True, **kwargs)

    def hyperfine_addbackground(self, xarr, Tbackground=2.73, Tex=5.0, tau=0.1,
                                xoff_v=0.0, width=1.0, return_tau=False,
                                **kwargs):
        """
        Identical to hyperfine, but adds Tbackground as a constant continuum
        level
        """
        if return_tau:
            raise ValueError("Cannot return tau when adding a continuum background.")
        return (self.hyperfine(xarr, Tbackground=Tbackground, Tex=Tex, tau=tau,
                               xoff_v=xoff_v, width=width, return_tau=False,
                               **kwargs)
                + Tbackground)

    def hyperfine_background(self, xarr, Tbackground=2.73, Tex=5.0, tau=0.1,
                             xoff_v=0.0, width=1.0, return_tau=False,
                             **kwargs):
        """
        Identical to hyperfine, but with Tbackground free.  Assumes already
        background-subtracted
        """
        if return_tau:
            raise ValueError("Cannot return tau when adding a continuum background.")
        return self.hyperfine(xarr, Tbackground=Tbackground, Tex=Tex, tau=tau,
                              xoff_v=xoff_v, width=width, return_tau=False,
                              **kwargs)

    def hyperfine(self, xarr, Tex=5.0, tau=0.1, xoff_v=0.0, width=1.0,
                  return_hyperfine_components=False, Tbackground=2.73, amp=None,
                  return_tau=False, tau_total=None, vary_hyperfine_tau=False,
                  vary_hyperfine_width=False):
        """
        Generate a model spectrum given an excitation temperature, optical
        depth, offset velocity, and velocity width.

        Parameters
        ----------
        return_tau : bool
            If specified, return just the tau spectrum, ignoring Tex
        tau_total : bool
            If specified, use this *instead of tau*, and it tries to normalize
            to the *peak of the line*
        vary_hyperfine_tau : bool
            If set to true, allows the hyperfine transition amplitudes to vary and
            does not use the line_strength_dict.  If set, `tau` must be a dict
        """

        # Convert X-units to frequency in Hz
        try:
            xarr = xarr.as_unit('Hz').value
        except AttributeError:
            xarr = xarr.to('Hz').value

        # Ensure parameters are scalar / have no extra dims
        if not np.isscalar(Tex): Tex = Tex.squeeze()
        if not np.isscalar(xoff_v): xoff_v = xoff_v.squeeze()
        if vary_hyperfine_width:
            if not isinstance(width, dict):
                raise TypeError("If varying the amplitude of the hyperfine lines, must specify tau as a dict")
        else:
            if not np.isscalar(width): width = width.squeeze()
        if vary_hyperfine_tau:
            if not isinstance(tau, dict):
                raise TypeError("If varying the amplitude of the hyperfine lines, must specify tau as a dict")
        else:
            if not np.isscalar(tau): tau = tau.squeeze()

        # Generate an optical depth spectrum as a function of the X-axis
        tau_nu_cumul = np.zeros(len(xarr))
        # Error check: inputing NANs results in meaningless output - return without computing a model
        if (np.any(np.isnan((Tex,xoff_v))) or
           ((not vary_hyperfine_tau) and np.isnan(tau)) or
           ((not vary_hyperfine_width) and np.isnan(width))):
            if return_hyperfine_components:
                return [tau_nu_cumul] * len(self.line_names)
            else:
                return tau_nu_cumul

        if tau_total is not None:
            tau = 1

        components =[]
        for linename in self.line_names:
            voff_lines = np.array(self.voff_lines_dict[linename])

            lines = (1-voff_lines/ckms)*self.freq_dict[linename]
            if not vary_hyperfine_width and width == 0:
                tau_nu = xarr*0
            else:
                if vary_hyperfine_width:
                    nuwidth = np.abs(width[linename]/ckms*lines)
                else:
                    nuwidth = np.abs(width/ckms*lines)
                nuoff = xoff_v/ckms*lines
                if vary_hyperfine_tau:
                    tau_line = tau[linename]
                else:
                    # the total optical depth, which is being fitted, should be the sum of the components
                    tau_line = (tau * np.array(self.line_strength_dict[linename])/
                                np.array(self.relative_strength_total_degeneracy[linename]))

                tau_nu = np.array(tau_line *
                                  np.exp(-(xarr+nuoff-self.freq_dict[linename])**2 /
                                          (2.0*nuwidth**2)))
                tau_nu[tau_nu!=tau_nu] = 0 # avoid nans
            components.append(tau_nu)
            tau_nu_cumul += tau_nu

        # add a list of the individual 'component' spectra to the total components...

        if tau_total is not None:
            tau_max = tau_nu_cumul.max() # danger of undersampling...
            tau_nu_cumul *= tau_total/tau_max
            for c in components:
                c *= tau_total/tau_max

        if return_hyperfine_components:
            if return_tau:
                return components
            elif amp is None:
                return (1.0-np.exp(-np.array(components)))*(Tex-Tbackground)
            else:
                comps = (1.0-np.exp(-np.array(components)))*(Tex-Tbackground)
                return comps/comps.max() * amp

        if return_tau:
            return tau_nu_cumul
        else:

            # This is not the full equation of radiative transfer, but a
            # background-subtracted version
            # With "background" function B_nu = CMB, S_nu = absorber, and I_nu = received:
            # I_nu = B_nu * exp(-tau) + (1-exp(-tau)) * S_nu
            # This is a very good approximation for Rohlfs & Wilson eqn 15.29:
            #spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*(Tex-Tbackground)

            # this is the exact version of 15.29
            T0 = hoverk * xarr

            # division by zero should raise an exception, but zero-background
            # is allowed and just turns this term to zero=(exp(-inf))
            if Tbackground > 0:
                background_term = 1/(np.exp(T0/Tbackground)-1)
            else:
                background_term = 0

            with np.errstate(divide='raise'):
                # division by zero is disallowed and should raise an exception
                spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*T0*(1/(np.exp(T0/Tex)-1) - background_term)

            # This is the equation of radiative transfer using the RJ definitions
            # (eqn 1.37 in Rohlfs)
            # It is identical, except without T_background subtracted
            # spec = Tex+(np.exp(-np.array(tau_nu_cumul)))*(Tbackground-Tex)

            if amp is None:
                return spec
            else:
                return spec/spec.max() * amp
