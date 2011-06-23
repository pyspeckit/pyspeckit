import numpy as np
import itertools
import speclines

opt_waves = np.array(
             [1033.30, 1215.67, 1239.42, 1305.53, 1335.52, 1399.8, 1545.86, 1640.4, 
              1665.85, 1857.4, 1908.27, 2326.0, 2439.5, 2800.32, 3346.79, 3426.85,
              3728.30, 3798.976, 3836.47, 3889.0, 3934.777, 3969.588, 4072.3, 
              4102.89, 4305.61, 4341.68, 4364.436, 4862.68, 4960.295, 5008.240,
              5176.7, 5895.6, 6302.046, 6365.536, 6549.86, 6564.61, 6585.27,
              6707.89, 6718.29, 6732.67])
opt_lines = np.array(
             ['OVI', 'Ly_alpha', 'NV', 'OI', 'CII', 'SiIV+OIV', 'CIV', 'HeII', 
              'OIII', 'AlIII', 'CIII', 'CII', 'NeIV', 'MgII', 'NeV', 'NeV', 'OII', 
              'H_theta', 'H_eta', 'HeI', 'K', 'H', 'SII', 'H_delta', 'G', 'H_gamma',
              'OIII', 'H_beta', 'OIII', 'OIII', 'Mg', 'Na', 'OI', 'OI', 'NII', 'H_alpha',
              'NII', 'Li', 'SII', 'SII'])
              
optical_lines = {'line': opt_lines, 'wavelength': opt_waves, 'units': 'Angstrom', 'vac': True}

"""
To test:

import spectrum
spec = spectrum.Spectrum('sample_sdss.txt')
spec.plotter(xmin = 6400, xmax = 6800)
spec.specfit(guesses = [20, 6718.29, 5, 100, 6564.614, 20, 50, 6585.27, 20, 20, 6732.67, 5, 50, 6549.86, 5])
spec.specfit.measurements.identify()


"""

class Measurements(object):
    def __init__(self, specfit, z = None, d = None, xunits = None):
        """
        This is called after a fit is run.  It will inherit the specfit object and derive as much as it can from modelpars.
        
        Notes: If z (redshift) or d (distance) are present, we can compute ingrated line luminosities rather than just fluxes.
        
        Currently will only work with Gaussians. to generalize:
            1. make sure we manipulate modelpars correctly, i.e. read in entries corresponding to wavelength/frequency/whatever correctly.
            
        """
        
        print optical_lines.name
            
        # Inherit specfit object    
        self.specfit = specfit  
        
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
        self.refpos = optical_lines['wavelength']
        
        
    def identify(self):
        """
        Determine identify of lines in self.fitpars.  Fill entries of self.lines dictionary.
        
        probs if two lines have very similar wavelengths in fit - account for this somehow
        """    
        
        self.IDresults= []
        odiff = np.diff(self.obspos)
        condition = (self.refpos >= 0.9 * min(self.obspos)) & (self.refpos <= 1.1 * max(self.obspos))   # Speeds things up
        refpos = self.refpos[condition]
        
        combos = itertools.combinations(refpos, self.Nlines)        
        for i, combo in enumerate(combos):
            rdiff = np.diff(combo)
            self.IDresults.append((np.sum(np.abs(odiff - rdiff)), combo))
            
        # Pick best solution
        MINloc = np.argmin(zip(*self.IDresults)[0])  # Location of best solution
        ALLloc = []                                  # Locations of best fit lines in reference dictionary
                
        for element in self.IDresults[MINloc][1]: ALLloc.append(np.argmin(np.abs(optical_lines['wavelength'] - element)))
        for i, element in enumerate(ALLloc): 
            line = optical_lines['line'][element]
            self.lines[line] = {}
            self.lines[line]['modelpars'] = self.modelpars[i]
    
    def derive(self):
        """
        Calculate luminosity and FWHM for all spectral lines.
        """            
        
        for line in self.lines.keys():
            
            # If only single component fit
            if len(self.lines[line].shape) == 1:
                self.fwhm[line] = 2. * np.sqrt(2. * np.log(2.)) * self.lines[line][2]
            
            # If multi-component fit
            else:
                self.fwhm[line] = self.find_fwhm(self.lines[line])
                
            try: self.lum[line] = self.find_lum(self.lines[line], self.d_L)
            except TypeError: pass
    
    def integrated_flux(self, pars):                                                                       
        """                                                                                                
        Calculate integrated flux of emission line.  Works for multi-component fits too.  Unnormalized.    
        """                                                                                                
                                                                                                           
        flux = 0                                                                                           
        for i in xrange(len(pars) / 3):                                                                    
            flux += np.sqrt(2. * np.pi) * pars[0 + 3*i] * pars[2 + 3*i]                                    
                                                                                                           
        return flux                                                                                        
                                                                                                           
    def integrated_luminosity(self, pars):                                                                 
        """                                                                                                
        Determine luminosity of line (need distance).                                                      
        """                                                                                                
                                                                                                           
        lum = 0                                                                                            
        for i in xrange(len(pars) / 3):                                                                    
            lum += self.fluxnorm * self.integrated_flux(pars) * 4. * np.pi * self.d_L**2                   
                                                                                                           
        return lum                                                                                         
        
    def integrated_fwhm(self, pars):
        """
        Determine full-width at half maximum for multi-component fit numerically.
        """
        
        fmax = np.sum(zip(*pars)[0])   # full max
        hmax = 0.5 * fmax              # half max
                
        start = np.mean(zip(*pars)[1])
        maxnow = self.multigauss(start, pars)
        
        wave = start + 1
        while maxnow > hmax:
            maxnow = self.multigauss(wave, pars)
            wave += 1
                        
        return 2. * (np.interp(hmax, [self.multigauss(wave - 1, pars), self.multigauss(wave - 2, pars)], [wave - 1, wave - 2]) - start)
            
        
        
            
        

