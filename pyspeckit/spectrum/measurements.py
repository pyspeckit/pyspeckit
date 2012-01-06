import numpy as np
import itertools, cosmology

"""
To test:

import spectrum
spec = spectrum.Spectrum('sample_sdss.txt')
spec.plotter(xmin = 6400, xmax = 6800)
spec.specfit(guesses = [20, 6718.29, 5, 100, 6564.614, 20, 50, 6585.27, 20, 20, 6732.67, 5, 50, 6549.86, 5])
spec.measure()
"""

cm_per_mpc = 3.08568e+24

class Measurements(object):
    def __init__(self, Spectrum, z=None, d=None, xunits=None, fluxnorm=None,
            miscline=None, misctol=10, ignore=None, derive=True, debug=False):
        """
        This can be called after a fit is run.  It will inherit the specfit
        object and derive as much as it can from modelpars.  Just do:
        spec.measure(z, xunits, fluxnorm)
        
        Notes: If z (redshift) or d (distance) are present, we can compute
        integrated line luminosities rather than just fluxes.  Provide distance
        in cm.
        
        Currently will only work with Gaussians.  To generalize:
            1. make sure we manipulate modelpars correctly, i.e. read in
            entries corresponding to wavelength/frequency/whatever correctly.
            
        misclines = dictionary
            miscline = {{'name': H_alpha', 'wavelength': 6565, 'etc': 0}, {}}
            
        """
        self.debug = debug
                    
        # Inherit specfit object    
        self.specfit = Spectrum.specfit
        self.speclines = Spectrum.speclines
                        
        # Bit of a hack - help identifying unmatched lines
        self.miscline = miscline        
        self.misctol = misctol
                
        # Flux units in case we are interested in line luminosities or just having real flux units
        if fluxnorm is not None: self.fluxnorm = fluxnorm
        else: self.fluxnorm= 1
                
        # This is where we'll keep our results                        
        self.lines = {}
        
        # Read in observed wavelengths
        tmp1 = np.reshape(self.specfit.modelpars, (len(self.specfit.modelpars) / 3, 3))
        tmp2 = np.reshape(self.specfit.modelerrs, (len(self.specfit.modelerrs) / 3, 3))
        
        if ignore is not None:
            tmp1 = np.delete(tmp1, ignore, 0)
            tmp2 = np.delete(tmp2, ignore, 0)
                
        order = np.argsort(zip(*tmp1)[1])
        self.obspos = np.sort(list(zip(*tmp1)[1]))
        self.Nlines = len(self.obspos)
                
        # Read in modelpars and modelerrs, re-organize so they are 2D arrays sorted by ascending wavelength
        self.modelpars = np.zeros_like(tmp1)
        self.modelerrs = np.zeros_like(tmp2)
        for i, element in enumerate(order): 
            self.modelpars[i] = tmp1[element]
            self.modelerrs[i] = tmp2[element]
                                                                      
        # Read in appropriate list of reference wavelengths/frequencies/whatever
        self.reflines = self.speclines.optical.optical_lines
        self.refpos = self.reflines['xarr']
        self.refname = self.reflines['name']
                        
        # If distance or redshift has been provided, we can compute luminosities from fluxes
        if d is not None: self.d = d
        else: self.d = None
        if z is not None: 
            self.cosmology = cosmology.Cosmology()
            self.d = self.cosmology.LuminosityDistance(z) * cm_per_mpc
            self.redshift = z
        else:
            self.redshift = 0.0
            
        self.identify()
        if derive: self.derive()
        
    def identify(self):
        """
        Determine identity of lines in self.modelpars.  Fill entries of self.lines dictionary.
        
        Note: This method will be infinitely slow for more than 10 or so lines.
        """    
        
        # List to keep track of line identification.  Each entry is (cost, (line1, line2, lin3,...))
        self.IDresults = []
        
        # Spacing between observed lines (odiff) and reference lines (rdiff)
        self.odiff = np.abs(np.diff(self.obspos))
        self.rdiff = np.abs(np.diff(self.refpos))
        self.rdmin = 0.5 * min(self.rdiff)
        
        # If lines have multiple components...
        if np.any(self.odiff) < self.rdmin:
            where = np.ravel(np.argwhere(self.odiff < self.rdmin))
            odiff = np.delete(self.odiff, where)
            multi = True
        else: 
            where = 0
            odiff = self.odiff
            multi = False

        # need to account for redshift if self.redshift is set
        # WRONG! Assume REST frame, do shifting elsewhere...
        refpos = self.refpos #* (1.0+self.redshift)
        refname = self.refname     
                            
        # Don't include elements of reference array that are far away from the observed lines (speeds things up)
        condition = (refpos >= 0.9 * min(self.obspos)) & (refpos <= 1.1 * max(self.obspos)) 
        refpos = refpos[condition]
        refname = refname[condition]
                        
        combos = itertools.combinations(refpos, self.Nlines - len(where))        
        for i, combo in enumerate(combos):
            rdiff = np.diff(combo)
            self.IDresults.append((np.sum(np.abs(odiff - rdiff)), combo))
                        
        # Pick best solution
        MINloc = np.argmin(zip(*self.IDresults)[0])  # Location of best solution
        ALLloc = []                                  # x-values of best fit lines in reference dictionary
                                
        # Determine indices of matched reference lines       
        for element in self.IDresults[MINloc][1]: 
            ALLloc.append(np.argmin(np.abs(refpos - element)))                    
        
        # Fill lines dictionary            
        for i, element in enumerate(ALLloc): 
            line = refname[element]
            self.lines[line] = {}
            loc = np.argmin(np.abs(self.obspos - refpos[element]))                
            self.lines[line]['modelpars'] = list(self.modelpars[loc])            
            self.lines[line]['modelerrs'] = list(self.modelerrs[loc])            
                                                               
        # Track down odd lines (i.e. broad components of lines already identified)
        # This won't yet work for lines that are truly unidentified
        if len(ALLloc) < self.Nlines:
                                                
            # Eliminate all modelpars/errs that belong to lines that were identified
            tmp1 = list(np.ravel(self.modelpars))
            tmp2 = list(np.ravel(self.modelerrs))
            for key in self.lines.keys():
                for element in self.lines[key]['modelpars']: 
                    loc = np.argmin(np.abs(element - tmp1))
                    tmp1 = np.delete(tmp1, loc)
                    tmp2 = np.delete(tmp2, loc)
             
            # Loop over unmatched modelpars/errs, find name of unmatched line, extend corresponding dict entry
            if self.miscline is None:                                          
                try:  
                    for i, x in enumerate(zip(*tmp1)[1]):    
                        loc = np.argmin(np.abs(ALLloc - x))
                        line = refname[loc]
                        self.lines[line]['modelpars'].extend(tmp1[i:i+3])
                        self.lines[line]['modelerrs'].extend(tmp2[i:i+3])
                except TypeError:
                    loc = np.argmin(np.abs(tmp1[1] - refpos))                       
                    line = refname[loc]
                    self.lines[line]['modelpars'].extend(tmp1)
                    self.lines[line]['modelerrs'].extend(tmp2)
            
            # If we've know a-priori which lines the unmatched lines are likely to be, use that information        
            else:
                
                for i, miscline in enumerate(self.miscline):
                    try:  
                        for j, x in enumerate(zip(*tmp1)[1]):    
                            if abs(x - self.lines[miscline['name']]['modelpars'][1]) < self.misctol[i]:
                                self.lines[line]['modelpars'].extend(tmp1[j:j+3])
                                self.lines[line]['modelerrs'].extend(tmp2[j:j+3])
                    except TypeError:
                        if abs(tmp1[1] - self.lines[self.miscline[0]]['modelpars'][1]) < self.misctol:
                            self.lines[self.miscline[0]]['modelpars'].extend(tmp1)
                            self.lines[self.miscline[0]]['modelerrs'].extend(tmp2)
                            break #?
        
        
                                              
        self.separate() 
                            
    def derive(self):
        """
        Calculate luminosity and FWHM for all spectral lines.
        """            
        
        for line in self.lines.keys():
            if self.debug:
                print "Computing parameters for line %s" % line
            
            self.lines[line]['fwhm'] = self.compute_fwhm(self.lines[line]['modelpars'])
            self.lines[line]['flux'] = self.compute_flux(self.lines[line]['modelpars'])
            self.lines[line]['amp'] = self.compute_amplitude(self.lines[line]['modelpars'])
            self.lines[line]['pos'] = self.lines[line]['modelpars'][1]
            
            if self.d is not None:
                self.lines[line]['lum'] = self.compute_luminosity(self.lines[line]['modelpars'])            
                
    def separate(self):
        """
        For multicomponent lines, separate into broad and narrow components (assume only one of components is narrow).
        """
        
        for key in self.lines.keys():
            modpars = self.lines[key]['modelpars']
            moderrs = self.lines[key]['modelerrs']
            if len(modpars) > 3:
                modpars2d = np.reshape(modpars, (len(modpars) / 3, 3))
                moderrs2d = np.reshape(moderrs, (len(moderrs) / 3, 3))
                sigma = zip(*modpars2d)[2]
                minsigma = min(np.abs(sigma))
                i_narrow = list(np.abs(sigma)).index(minsigma)
            else: continue
                        
            self.lines["{0}_N".format(key)] = {}         
            self.lines["{0}_N".format(key)]['modelpars'] = []   
            self.lines["{0}_N".format(key)]['modelerrs'] = []   
            self.lines["{0}_B".format(key)] = {}            
            self.lines["{0}_B".format(key)]['modelpars'] = [] 
            self.lines["{0}_B".format(key)]['modelerrs'] = [] 
                        
            for i, arr in enumerate(modpars2d):
                if i == i_narrow: 
                    self.lines["{0}_N".format(key)]['modelpars'] = arr
                    self.lines["{0}_N".format(key)]['modelerrs'] = moderrs2d[i]
                else: 
                    self.lines["{0}_B".format(key)]['modelpars'].extend(arr)
                    self.lines["{0}_B".format(key)]['modelerrs'].extend(moderrs2d[i])
                    
    def compute_flux(self, pars):                                                                       
        """                                                                                                
        Calculate integrated flux of emission line.  Works for multi-component fits too.  Unnormalized.    
        """                                                                                                
                                                                                                           
        flux = 0                                                                                           
        for i in xrange(len(pars) / 3): flux += np.sqrt(2. * np.pi) * pars[3 * i] * abs(pars[2 + 3 * i])
                                                                                                           
        return flux * self.fluxnorm
        
    def compute_amplitude(self, pars):
        """
        Calculate amplitude of emission line.  Should be easy - add multiple components if they exist.
        Currently assumes multiple components have the same centroid.
        """
        
        amp = 0
        for i in xrange(len(pars) / 3): amp += pars[3 * i]
        return amp * self.fluxnorm
                                                                                                           
    def compute_luminosity(self, pars):                                                                 
        """                                                                                                
        Determine luminosity of line (need distance and flux units).                                                      
        """                                                                                                
                                                                                                           
        lum = 0                                                                                            
        for i in xrange(len(pars) / 3): lum += self.compute_flux(pars) * 4. * np.pi * self.d**2                   
        return lum                                                                                         
        
    def compute_fwhm(self, pars):
        """
        Determine full-width at half maximum for multi-component fit numerically, or analytically if line
        has only a single component.  Uses bisection technique for the former with absolute tolerance of 1e-4.
        """

        if len(pars) == 3: 
            return 2. * np.sqrt(2. * np.log(2.)) * abs(pars[2])
        else:
            atol = 1e-4
            pars2d = np.reshape(pars, (len(pars) / 3, 3))
            start = zip(*pars2d)[1][0]                    # start at central wavelength of first component
            
            # If the centroids are exactly the same for all components, we know the peak, and peak position
            if np.allclose(zip(*pars2d)[1], atol): 
                fmax = np.sum(zip(*pars2d)[0])            
                
            # Otherwise, we have to figure out where the multicomponent peak is
            else:    
                f = lambda x: self.specfit.fitter.slope(x)
                xfmax = self.bisection(f, start)
                fmax = self.specfit.fitter.n_modelfunc(pars)(np.array([xfmax, xfmax]))[0]
                            
            hmax = 0.5 * fmax    
                                
            # current height relative to half max - we want to minimize this function.  Could be asymmetric.
            f = lambda x: self.specfit.fitter.n_modelfunc(pars)(np.array([x])) - hmax                   
            xhmax1 = self.bisection(f, start)
            xhmax2 = self.bisection(f, start + (start - xhmax1))
                                                    
            return abs(xhmax2 - xhmax1)      
            
    def bisection(self, f, x_guess):
        """
        Find root of function using bisection method.  Absolute tolerance of 1e-4 is being used.
        """

        x1, x2 = self.bracket_root(f, x_guess)
        
        # Narrow bracketed range with bisection until tolerance is met
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
        
    def to_tex(self):
        """
        Write out fit results to tex format.
        """    
        
        pass
            
        

