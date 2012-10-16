import numpy as np
import pyspeckit.spectrum.models.inherited_voigtfitter as IV

xarr = np.linspace(-100,100,1000)
dx = np.diff(xarr).mean()

V1 = IV.voigt(xarr,1,0,1,1,normalized=False)
V2 = IV.voigt(xarr,1,0,1,1,normalized=True)

assert np.sqrt(2*np.pi) - 0.05 < V1.sum()*dx < np.sqrt(2*np.pi) + 0.05 
assert 0.99 < V2.sum()*dx < 1.01
