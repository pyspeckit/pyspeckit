"""
Silly example to showcase dsypher's reddening tools
"""
import pyspeckit

sp = pyspeckit.OpticalSpectrum('sn2009ip_halpha.fits')

sp.plotter()

sp.deredden(a_v=1)
sp.plotter(clear=False,color='b')

sp.redden(a_v=2)
sp.plotter(clear=False,color='r')
