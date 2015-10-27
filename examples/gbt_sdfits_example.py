from __future__ import print_function
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
from pylab import *;import numpy,scipy,matplotlib;
import pyspeckit
from pyspeckit.spectrum.readers import gbt

sdfitsfile = '/Users/adam/observations/gbt/h2co_pilot/AGBT09C_049_02.raw.acs.fits'
gbtdata = pyfits.open(sdfitsfile)
objectname = 'G32.80+0.19'
bintable=gbtdata[1]
A049 = gbt.GBTSession(bintable)
whobject = bintable.data['OBJECT'] == objectname
blocks = A049.load_target(objectname).blocks
onMoff1 = blocks['A9OFF1'] - blocks['A9OFF2']
onMoff1on = blocks['A9ON1'] - blocks['A9ON2']

onoffs = dict((name[:-1],blocks[name[:-1]+'1'] - blocks[name[:-1]+'2']) for name in blocks if name[-1]=='1')

av1 = onMoff1.average()
av1.plotter()
av2 = onMoff1on.average()
av2.plotter(axis=av1.plotter.axis,clear=False,color='b')

calonA9 = blocks['A9ON1'].average()
caloffA9 = blocks['A9OFF1'].average()
print("tsys: ",gbt.dcmeantsys(calonA9, caloffA9, calonA9.header['TCAL']))

G32 = reduced_nods = A049.reduce_target(objectname, verbose=True)
print(G32)

G32.average_IFs(debug=True)
G32.average_pols()
G32.spectra['if0'].plotter()

G37 = A049.reduce_target('G37.87-0.40')
print(G37)
G37.average_IFs(debug=True)
G37.spectra['if0'].plotter()

name = 'A9'
num = '1'
av2 = blocks[name+'OFF'+num].average(debug=True)
#while av2.data.mean() > 0.01:
#    av2 = blocks[name+'OFF'+num].average(debug=True)
#    print(av2)

colors = ['k','b','r','g']
for fd,col in zip(['A9','A13','C25','C29'],colors):
    G32.reduced_scans[fd].plotter(figure=figure(4), clear=False,color=col)

G32['if1'].plotter(figure=figure(4), color='magenta', clear=False)
G32['if1fd1'].plotter(figure=figure(4), color='orange', clear=False)
G32['if1fd2'].plotter(figure=figure(4), color='cyan', clear=False)

G32b = pyspeckit.Spectrum('/Users/adam/work/h2co/data/pilot/G32.80+0.19_h2co_Jy.fits')
G32b.xarr.convert_to_unit('Hz')
((G32b+G32b.header['CONTINUU'])*1.91).plotter(figure=figure(4),clear=False,color='purple')
