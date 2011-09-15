""" 
cosmology.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2010-03-01.

Description: Simple cosmology calculator.

Notes: 
    -Everything here uses cgs.
    -I have assumed a flat universe for all calculations, i.e. OmegaCurvatureNow = 0.0.
    -WMAP VII cosmological parameters by default.

"""

import numpy as np
try:
    from scipy.integrate import romberg
    scipyOK=True
except:
    scipyOK=False

c = 29979245800.0
G = 6.673*10**-8
km_per_mpc = 3.08568 * 10**13 * 10**6
cm_per_mpc = 3.08568 * 10**13 * 10**5 * 10**6
sqdeg_per_std = (180.0**2) / (np.pi**2)

class Cosmology:
    def __init__(self, OmegaMatterNow = 0.272, OmegaLambdaNow = 0.728, OmegaBaryonNow = 0.044,
        h_70 = 0.702, sigma_8 = 0.807):
        
        self.OmegaMatterNow = OmegaMatterNow
        self.OmegaLambdaNow = OmegaLambdaNow
        self.OmegaBaryonNow = OmegaBaryonNow
        self.h_70 = h_70
        self.sigma_8 = sigma_8
        self.OmegaCDMNow = self.OmegaMatterNow - self.OmegaBaryonNow
        
        self.HubbleParameterNow = h_70 * 100 / km_per_mpc
        self.CriticalDensityNow = 3 * self.HubbleParameterNow**2 / 8 / np.pi / G

    def ScaleFactor(self, z):
        return 1.0 / (1.0 + z)
        
    def EvolutionFunction(self, z):
        return np.sqrt(self.OmegaMatterNow * (1.0 + z)**3  + self.OmegaLambdaNow)
        
    def HubbleParameter(self, z):	
        return self.HubbleParameterNow * np.sqrt(self.OmegaMatterNow * (1.0 + z)**3 + 
            self.OmegaLambdaNow) 
    
    def OmegaMatter(self, z):
        return self.OmegaMatterNow * (1.0 + z)**3 / self.EvolutionFunction(z)**2
    
    def OmegaLambda(self, z):
	    return self.OmegaLambdaNow / self.EvolutionFunction(z)**2
    
    def MeanMatterDensity(self, z):
        return self.OmegaMatter(z) * self.CriticalDensity(z)
        
    def MeanBaryonDensity(self, z):
        return (self.OmegaBaryonNow / self.OmegaMatterNow) * self.MeanMatterDensity(z)
    
    def CriticalDensity(self, z):
        return (3.0 * self.HubbleParameter(z)**2) / (8.0 * np.pi * G)
        
    def LookbackTime(self, z_i, z_f):
        """
        Returns lookback time from z_i to z_f in seconds, where z_i < z_f.
        """
        Integrand = lambda z: (1.0 / (1.0 + z) / self.EvolutionFunction(z))
        
        if scipyOK:
            return (romberg(Integrand, z_i, z_f) / self.HubbleParameterNow)    
        
    def TimeToRedshiftConverter(self, z_i, dt):
        """
        Given a redshift and elapsed time dt (in seconds), returns final redshift after dt has passed.
        Uses the high redshift approximation (OmegaLambda -> 0).
        """
        return ((1. + z_i)**(-3. / 2.) + (3. * self.HubbleParameterNow * np.sqrt(self.OmegaMatterNow) * \
            dt / 2.))**(-2. / 3.) - 1.      
    
    def ComovingRadialDistance(self, z_i, z_f):
        """
        Returns comoving radial distance (comoving transverse distance for Omega_K -> 0).
        """
        
        integrand = lambda z: 1. / self.EvolutionFunction(z)
        return romberg(integrand, z_i, z_f) * c / self.HubbleParameterNow / cm_per_mpc
        
    def LuminosityDistance(self, z_f):
        """
        Returns luminosity distance in Mpc.  Assumes we mean distance from us (z = 0).
        """
        
        return (1. + z_f) * self.ComovingRadialDistance(0., z_f)
        
    
