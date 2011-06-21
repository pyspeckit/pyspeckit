import numpy as np
from speclines import *

class Measurements:
    def __init__(self, specfit, z = None, d = None):
        """
        This is called after a fit is run.  It will inherit the specfit object and derive as much as it can from modelpars.
        
        Notes: If z (redshift) or d (distance) are present, we can compute ingrated line luminosities rather than just fluxes.
        
        Currently will only work with Gaussians. to generalize:
            1. make sure we manipulate modelpars correctly, i.e. read in entries corresponding to wavelength/frequency/whatever correctly.
            
        """
            
        self.specfit = specfit
        self.modelpars = np.reshape(self.specfit.modelpars, (len(self.specfit.modelpars), self.specfit.npars))
        self.obspos = zip(*self.modelpars)[1]
        
        # Depending on units of xarr, read in set of reference spectral line positions
        #self.refpos = 
        
        self.identify()
        
        #self.derive()
        
    def identify(self):
        """
        Determine identify of lines in self.fitpars.  Fill entries of self.lines dictionary.
        """    
        
        
    
    
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
            
        
        
            
        

