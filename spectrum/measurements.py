"""
measurements.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-04-03.

Description: Read in list of Gaussian fit parameters for spectral lines, construct spectrum object.

Example:
from speclines import *
import h5py
f = h5py.File('J004719.39_144212.6_pca_out.hdf5')
fits = f['modelpars'].value
spec = spectrum(fits)
 
"""

import misc, scipy
import numpy as np
from cosmology import *
from constants import *

class measurements(object):
    def __init__(self, SpecFit):
    
    #fitpars, speclinefile = None, 
        #fluxnorm = 1e-17, z = None, d_L = None, linesep = True):
        """
        Initialize a spectrum object, i.e. identify spectral lines that have been fit, and calculate
        their FWHMs and luminosities.
            
            fitpars: list
            ex: fitpars = [A1, x1, sigma1, A2, x2, sigma2]
            
            fluxnorm: float
            ex: fluxnorm = 1e-17 # erg / s for SDSS spectra (default)
            
            linesep: bool
            Whether or not to separate multi-component fits into their broad and narrow components.
            
            d_L = luminosity distance in Mpc
            
        """
        self.SpecFit = SpecFit
        self.cosmology = cosmology()
                
        #self.fitpars = fitpars
        #self.z = z
        #self.d_L = d_L
        #self.linesep = linesep
        self.Nlines = len(self.SpecFit.modelpars) / 3
        self.lines = {}
        self.lum = {}       # erg / s
        self.fwhm = {}      # wavelength units currently
        
        #if self.d_L is not None: self.d_L *= cm_per_mpc
        #else:
        #    if self.z is None: print "Redshift unknown.  Cannot derive line luminosities. "
        #    else: self.d_L = self.cosmology.LuminosityDistance(self.z) * cm_per_mpc
        #
        ## Read in common SDSS spectral lines
        #if speclinefile == None:
        #    speclinefile = "{0}/sdss_speclines.dat".format(os.environ.get("PYSPECKIT"))
        #
        #try: 
        #    self.wavelengths, self.names = misc.readtab(speclinefile) 
        #    self.wavelengths = np.array(self.wavelengths, float)
        #except IOError: print "Failed to read {0}.".format(speclinefile)
        #
        #self.identify()
        #
        ## Calculate derived parameters
        #self.derive()
        
    def identify(self):
        """
        Determine identify of lines in self.fitpars.  Fill entries of self.lines dictionary.
        """    
        
        maxdiff = max(np.diff(self.wavelengths))
        for i in np.arange(self.Nlines):
            diff = list(abs(self.wavelengths - self.fitpars[3 * i + 1]))
            minloc = diff.index(min(diff))
            
            if np.allclose(self.fitpars[3 * i + 1], self.wavelengths[minloc], atol = maxdiff):
                
                lines = np.array([self.fitpars[3 * i], self.fitpars[3 * i + 1], self.fitpars[3 * i + 2]])
                
                try: self.lines[self.names[minloc]].append(lines)
                except KeyError: self.lines[self.names[minloc]] = [lines]
    
        self.cleanup()
    
    def cleanup(self):
        """
        Go through spec.lines and split up lines like NII, SII, and OIII into separate entries, naming them according
        to wavelength.  Also, split broad and narrow lines apart (if linesep = True), and turn lists into numpy arrays.
        """
        
        # Split doublets, rename
        for line in self.lines.keys():
            
            if line == 'SII':
                tmp = np.ravel(self.lines[line])
                if tmp[1] < tmp[4]: i = 0
                else: i = 1
                    
                self.lines['SII_6716'] = self.lines[line][i]
                self.lines['SII_6731'] = self.lines[line][i - 1]
                del(self.lines['SII'])
                
            if line == 'NII':
                tmp = np.ravel(self.lines[line])
                if tmp[1] < tmp[4]: i = 0
                else: i = 1
                    
                self.lines['NII_6548'] = self.lines[line][i]    
                self.lines['NII_6583'] = self.lines[line][i - 1]
                del(self.lines['NII'])
                
            if line == 'OIII':
                tmp = np.ravel(self.lines[line])
                if tmp[1] < tmp[4]: i = 0
                else: i = 1
                    
                self.lines['OIII_4959'] = self.lines[line][i]    
                self.lines['OIII_5007'] = self.lines[line][i - 1]
                del(self.lines['OIII'])
                        
        # Split broad and narrow components
        if self.linesep:
                    
            for line in self.lines.keys(): 
                if len(np.array(self.lines[line]).shape) == 1: continue
                
                minsigma = np.inf
                for fit in self.lines[line]: minsigma = min(minsigma, fit[2])
                
                for fit in self.lines[line]:
                    if fit[2] == minsigma: self.lines["{0}_N".format(line)] = fit
                    else: 
                        try: self.lines["{0}_B".format(line)].append(fit)
                        except KeyError: self.lines["{0}_B".format(line)] = [fit]
                        
                del(self.lines[line])
                   
        # Convert lists to numpy arrays        
        for line in self.lines.keys():    
            if type(self.lines[line]) is list: self.lines[line] = np.array(self.lines[line])
    
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
                                       
    def gaussian(self, x, pars):
        """
        Return 1D Gaussian (already continuum subtracted).
        
            pars = [A, x0, sigma]
            
        """
        return pars[0] * np.exp(-(x - pars[1])**2 / 2. / pars[2]**2)                
     
    def multigauss(self, x, pars):
        """
        Multi-component Gaussian.
        """ 
        result = 0.
        for fit in pars: 
            result += self.gaussian(x, fit)
            
        return result
                    
    def find_lum(self, pars, d):
        """
        Determine integrated flux of line.  Unnormalized currently.
        """
                
        if len(pars.shape) == 1: 
            func = lambda x: self.gaussian(x, pars)
            return scipy.integrate.quad(func, pars[1] - 5. * pars[2], pars[1] + 5. * pars[2])[0] * 4. * np.pi * d**2
        else: 
            func = lambda x: self.multigauss(x, pars)
            return scipy.integrate.quad(func, pars[0][1] - 5. * pars[0][2], pars[0][1] + 5. * pars[0][2])[0] * 4. * np.pi * d**2
                
    def find_fwhm(self, pars):
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
            
        
        
            
        

