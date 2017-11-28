from __future__ import print_function
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from ..spectrum import models
from . import fith2co

path_to_data = "/Users/adam/work/h2co/radex/troscompt_grid_March2012"

texgrid1 = pyfits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex1.fits')
taugrid1 = pyfits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau1.fits')
texgrid2 = pyfits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tex2.fits')
taugrid2 = pyfits.getdata(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau2.fits')
hdr    = pyfits.getheader(path_to_data+'/1-1_2-2_T=5to55_lvg_troscompt_100square_opgrid_tau2.fits')
# # this deserves a lot of explanation:
# # models.formaldehyde.formaldehyde_radex is the MODEL that we are going to fit
# # models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
# # and annotations
# # all of the parameters after the first are passed to the model function
formaldehyde_radex_fitter = models.model.SpectralModel(
        models.formaldehyde.formaldehyde_radex_orthopara_temp, 6,
        parnames=['density','column','orthopara','temperature','center','width'],
        parvalues=[4,12,1.0,15.0,0,1],
        parlimited=[(True,True), (True,True), (True,True), (True,True), (False,False), (True,False)],
        parlimits=[(1,8), (11,16), (-3,np.log10(3.0)), (5,55), (0,0), (0,0)],
        parsteps=[0.01,0.01,0,0,0,0],
        fitunit='Hz',
        texgrid=((4,5,texgrid1),(14,15,texgrid2)),
        taugrid=((4,5,taugrid1),(14,15,taugrid2)),
        hdr=hdr,
        shortvarnames=("n","N",'OP','T',"v","\\sigma"),
        )
