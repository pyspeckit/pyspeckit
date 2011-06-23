import numpy as np
import itertools, cosmology

"""
To test:

import spectrum
spec = spectrum.Spectrum('sample_sdss.txt')
spec.plotter(xmin = 6400, xmax = 6800)
spec.specfit(guesses = [20, 6718.29, 5, 100, 6564.614, 20, 50, 6585.27, 20, 20, 6732.67, 5, 50, 6549.86, 5])
spec.specfit.measurements.identify()
spec.specfit.measurements.derive()
spec.specfit.measurements.lines

"""

cm_per_mpc = 3.08568e+24

class Measurements(object):
    def __init__(self, Spectrum, z = None, d = None, xunits = None):
        """
        This is called after a fit is run.  It will inherit the specfit object and derive as much as it can from modelpars.
        
        Notes: If z (redshift) or d (distance) are present, we can compute ingrated line luminosities rather than just fluxes.
            Provide distance in cm.
        
        Currently will only work with Gaussians. to generalize:
            1. make sure we manipulate modelpars correctly, i.e. read in entries corresponding to wavelength/frequency/whatever correctly.
            
        """
                    
        # Inherit specfit object    
        self.specfit = Spectrum.specfit
        self.speclines = Spectrum.speclines
        
        # If the header had units for flux values, use them!
        #if not Spectrum.units.strip(): self.fluxunits = 1
        #else: self.fluxunits = float(Spectrum.units)
                
        # This is where we'll keep our results                        
        self.lines = {}
        
        # Read in observed wavelengths
        tmp = np.reshape(self.specfit.modelpars, (len(self.specfit.modelpars) / 3, 3))
        self.obspos = np.sort(list(zip(*tmp)[1]))
        self.Nlines = len(self.obspos)
        
        # Read in modelpars, re-organize so it is a 2D array sorted by ascending wavelength
        self.modelpars = np.zeros_like(tmp)
        for i, element in enumerate(self.obspos):
            for j, arr in enumerate(tmp):
                if element == arr[1]: self.modelpars[i] = arr
                                   
        # Read in appropriate list of reference wavelengths/frequencies/whatever
        self.reflines = self.speclines.optical.optical_lines
        self.refpos = self.reflines['xarr']
        self.refname = self.reflines['name']
        
        # If distance or redshift has been provided, we can compute luminosities from fluxes
        if d is not None: self.d = d
        if z is not None: 
            self.cosmology = cosmology.Cosmology()
            self.d = self.cosmology.LuminosityDistance(z) * cm_per_mpc
        
    def identify(self):
        """
        Determine identify of lines in self.fitpars.  Fill entries of self.lines dictionary.
        
        Note: This method will be infinitely slow for more than 10 or so lines.
        """    
        
        self.IDresults= []
        self.odiff = np.abs(np.diff(self.obspos))
        self.rdiff = np.abs(np.diff(self.refpos))
        self.rdmin = 0.5 * min(self.rdiff)
        
        # If lines have multiple components...
        if np.any(self.odiff) < self.rdmin:
            where = np.argwhere(self.odiff < self.rdmin) 
            odiff = np.delete(self.odiff, where)
            multi = True
        else: 
            where = 0
            odiff = self.odiff
            multi = False        
        
        condition = (self.refpos >= 0.9 * min(self.obspos)) & (self.refpos <= 1.1 * max(self.obspos))   # Speeds things up
        refpos = self.refpos[condition]
                
        combos = itertools.combinations(refpos, self.Nlines - where)        
        for i, combo in enumerate(combos):
            rdiff = np.diff(combo)
            self.IDresults.append((np.sum(np.abs(odiff - rdiff)), combo))
            
        # Pick best solution
        MINloc = np.argmin(zip(*self.IDresults)[0])  # Location of best solution
        ALLloc = []                                  # x-values of best fit lines in reference dictionary
        
        # Fill lines dictionary        
        for element in self.IDresults[MINloc][1]: ALLloc.append(np.argmin(np.abs(self.refpos - element)))
        for i, element in enumerate(ALLloc): 
            line = self.refname[element]
            self.lines[line] = {}
            loc = np.argmin(np.abs(self.obspos - self.refpos[element]))                
            self.lines[line]['modelpars'] = list(self.modelpars[loc])
            
        # Track down odd lines
        if len(ALLloc) < self.Nlines:
            tmp = list(np.ravel(self.modelpars))
            for key in self.lines.keys():
                for element in self.lines[key]['modelpars']: tmp.pop(tmp.index(element))
                            
            try:  
                for i, x in enumerate(zip(*tmp)[1]):    
                    loc = np.argmin(np.abs(ALLloc - x))
                    line = self.refname[loc]
                    self.lines[line]['modelpars'].extend(tmp[i:i+3])
            except TypeError:
                loc = np.argmin(np.abs(tmp[1] - self.refpos))                       
                line = self.refname[loc]
                self.lines[line]['modelpars'].extend(tmp) 
                    
    def derive(self):
        """
        Calculate luminosity and FWHM for all spectral lines.
        """            
        
        for line in self.lines.keys():
            
            self.lines[line]['fwhm'] = self.compute_fwhm(self.lines[line]['modelpars'])
            self.lines[line]['flux'] = self.compute_flux(self.lines[line]['modelpars'])
            
            #if self.fluxunits != 1:
            #    try: self.lines[line]['lum'] = self.compute_luminosity(self.lines[line]['modelpars'])            
            #    except AttributeError: print 'Distance unknown.  Cannot compute line luminosities.'
    
    def compute_flux(self, pars):                                                                       
        """                                                                                                
        Calculate integrated flux of emission line.  Works for multi-component fits too.  Unnormalized.    
        """                                                                                                
                                                                                                           
        flux = 0                                                                                           
        for i in xrange(len(pars) / 3): flux += np.sqrt(2. * np.pi) * pars[3 * i] * pars[2 + 3 * i]                                    
                                                                                                           
        return flux                                                                                
                                                                                                           
    def compute_luminosity(self, pars):                                                                 
        """                                                                                                
        Determine luminosity of line (need distance).                                                      
        """                                                                                                
                                                                                                           
        lum = 0                                                                                            
        for i in xrange(len(pars) / 3): lum += self.compute_flux(pars) * 4. * np.pi * self.d**2                   
        return lum                                                                                         
        
    def compute_fwhm(self, pars):
        """
        Determine full-width at half maximum for multi-component fit numerically, or analytically if line
        has only a single component.  Uses bisection technique for the latter with absolute tolerance of 1e-4.
        """

        if len(pars) == 3: return 2. * np.sqrt(2. * np.log(2.)) * pars[2]
        else:
            
            atol = 1e-4
            pars2d = np.reshape(pars, (len(pars) / 3, 3))
            start = zip(*pars2d)[1][0]                    # start at central wavelength of first component
            
            # If the centroids are exactly the same for all components, we know the peak, and peak position
            if np.allclose(zip(*pars2d)[1], atol): 
                fmax = np.sum(zip(*pars2d)[0])            
                
            # Otherwise, we have to figure out where the multicomponent peak is
            else:    
                f = lambda x: self.specfit.fitter.multipeakgaussianslope(x, pars)
                xfmax = self.bisection(f, start)
                fmax = self.specfit.fitter.multipeakgaussian(xfmax, pars)
            
            hmax = 0.5 * fmax    
                
            # current height relative to half max - we want to minimize this function.  Could be asymmetric.
            f = lambda x: self.specfit.fitter.multipeakgaussian(x, pars) - hmax                   
            xhmax1 = self.bisection(f, start)
            xhmax2 = self.bisection(f, start + (start - xhmax1))
                                        
            return abs(xhmax2 - xhmax1)
            
    def bisection(self, f, x_guess):
        """
        Find root of function using bisection method.  Absolute tolerance of 1e-4 is being used.
        """

        x1, x2 = self.bracket_root(f, x_guess)
        
        # Narrow bracketed range with bisection until tolerance is met
        i = 0
        while abs(x2 - x1) > 1e-4:
            midpt = np.mean([x1, x2])
            fmid = f(midpt)
    
            if np.sign(fmid) < 0: x1 = midpt
            else: x2 = midpt
            
            if fmid == 0.0: break
            
        return x2    
        
    def bracket_root(self, f, x_guess, atol = 1e-4):
        """
        Bracket root by finding points where function goes from positive to negative.
        """
        
        f1 = f(x_guess)
        f2 = f(x_guess + 1)
        df = f2 - f1
                
        # Determine whether increasing or decreasing x_guess will lead us to zero
        if (f1 > 0 and df < 0) or (f1 < 0 and df > 0): sign = 1
        else: sign = -1
        
        # Find root bracketing points
        xpre = x_guess
        xnow = x_guess + sign 
        fpre = f1
        fnow = f(xnow)
        while (np.sign(fnow) == np.sign(fpre)):
            xpre = xnow
            xnow += sign * 0.1
            fpre = f(xpre)
            fnow = f(xnow)
                    
        x1 = min(xnow, xpre)
        x2 = max(xnow, xpre)
        
        if not np.all([np.sign(fpre), np.sign(fnow)]): 
            x1 -= 1e-4
            x2 += 1e-4
                                
        return x1, x2    
        
            
        

