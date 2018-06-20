from __future__ import print_function
# Test measurements class on an SDSS AGN spectrum

import pyspeckit
import astropy.units as u
import numpy as np

# Some lines
Hbeta = 4862.721
OIIIa = 4960.295
OIIIb = 5008.239
NIIa = 6549.86
Halpha = 6564.614
NIIb = 6585.27
SIIa = 6718.29
SIIb = 6732.68

# Offsets between SII lines in each complex
SIIb_off = SIIb - SIIa
OIIIb_off = OIIIb - OIIIa
NIIa_off_Ha = NIIa - Halpha
NIIb_off_Ha = NIIb - Halpha
SIIa_off_Ha = SIIa - Halpha 
SIIb_off_Ha = SIIb - Halpha 

# Guesses
narrow = 5.
broad = 30.

# Initialize spectrum object
spec = pyspeckit.Spectrum('sample_sdss.txt')
spec.unit = 'erg s^{-1} cm^{-2} \\AA^{-1}'
spec.xarr.set_unit = u.dimensionless_unscaled

# H-alpha
spec.specfit.selectregion(xmin = Halpha - 5, xmax = Halpha + 5)
ampHa = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])

# SII
spec.specfit.selectregion(xmin = SIIa - 20, xmax = SIIb + 20)
smallamp = np.max(spec.data[spec.specfit.xmin:spec.specfit.xmax])    

spec.specfit(guesses = [smallamp, SIIa, narrow, smallamp, SIIb, narrow],
        tied=['', '', 'p[-1]', '', 'p[1] + {0}'.format(SIIb_off), ''],
        negamp=False, quiet=True, multifit=None, show_components=True)
    
ampSIIa, lamSIIa, sigmaSII, ampSIIb, lamSIIb, sigmaSII = spec.specfit.modelpars

opz = lamSIIa / SIIa    # Calculate (1 + z)

# Now, use SII/OIII fit information to help fit NII-Halpha complex
spec.specfit.selectregion(xmin = NIIa - 100, xmax = NIIb + 100)

guesses = [2 * smallamp, NIIa * opz, sigmaSII, ampHa, Halpha * opz, sigmaSII,
          6 * smallamp, NIIb * opz, sigmaSII,
          ampSIIa, lamSIIa, sigmaSII, ampSIIb, lamSIIb, sigmaSII]
tied = ['', 'p[4] + {0}'.format(NIIa_off_Ha), 'p[{0}]'.format(len(guesses) - 1), 
        '', '', 'p[{0}]'.format(len(guesses) - 1), 
        '3 * p[0]', 'p[4] + {0}'.format(NIIb_off_Ha), 'p[{0}]'.format(len(guesses) - 1),
        '', 'p[4] + {0}'.format(SIIa_off_Ha), 'p[{0}]'.format(len(guesses) - 1), 
        '', 'p[4] + {0}'.format(SIIb_off_Ha), '']                                                
fixed = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

minp = []
maxp = []
lmin = []
lmax = []
for i, element in enumerate(guesses):
    minp.append(element - 0.2 * element)
    maxp.append(element + 0.2 * element)

    # Don't limit amplitude
    if i % 3 == 0:
        lmin.append(False)
        lmax.append(False)
    else:
        lmin.append(True)
        lmax.append(True)
        
# Plot, and do final fit
spec.plotter(xmin = NIIa - 100, xmax = SIIb + 30)

spec.specfit(guesses=guesses, tied=tied, fixed=fixed, negamp=False,
        limitedmin=lmin, limitedmax=lmax, minpars=minp, maxpars=maxp,
        show_components=True)
    
guesses.extend([50, Halpha * opz, 50])
tied.extend(['', 'p[4]', ''])
fixed.extend([0, 0, 0]) 

lmin.extend([True, False, False])
lmax.extend([False, False, False])
minp.extend([0, 0, 0]) 
maxp.extend([0, 0, 0]) 

spec.specfit(guesses = guesses, tied = tied, fixed = fixed, negamp = False,
    limitedmin = lmin, limitedmax = lmax, minpars = minp, maxpars = maxp)

guesses.extend([10, Halpha * opz, 25])
tied.extend(['', 'p[4]', ''])
fixed.extend([0, 0, 0])   

lmin.extend([True, False, False])
lmax.extend([False, False, False])
minp.extend([0, 0, 0]) 
maxp.extend([0, 0, 0])                                 

spec.plotter(xmin=NIIa - 100, xmax=SIIb + 30)

spec.specfit(guesses=guesses, tied=tied, fixed=fixed, negamp=False,
    limitedmin=lmin, limitedmax=lmax, minpars=minp, maxpars=maxp, annotate=False,
    show_components=True)
    
spec.plotter.refresh()

spec.measure(z = 0.05, fluxnorm = 1e-17)

# Overplot positions of lines and annotate
y = spec.plotter.ymax * 0.85
for i, line in enumerate(spec.measurements.lines.keys()):
    x = spec.measurements.lines[line]['modelpars'][1]
    spec.plotter.axis.plot([x]*2, [spec.plotter.ymin.value, spec.plotter.ymax.value], ls = '--', color = 'k')
    try: spec.plotter.axis.annotate(spec.speclines.optical.lines[line][-1], 
        (x, y), rotation = 90, ha = 'right', va = 'center')
    except KeyError: pass

spec.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
spec.plotter.axis.set_ylabel(r'Flux $(10^{-17} \mathrm{erg/s/cm^2/\AA})$')

spec.plotter.refresh()

print("Line   Pos   Flux (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)   Amplitude")
for line in spec.measurements.lines.keys():
    print(line, spec.measurements.lines[line]['pos'],
          spec.measurements.lines[line]['flux'],
          spec.measurements.lines[line]['fwhm'],
          spec.measurements.lines[line]['lum'],
          spec.measurements.lines[line]['amp'])

"""
Correct result:
Line   Pos   Flux (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)   Amplitude
NIIa 6551.28531669 1.45663617387e-15 5.53353450775 8.58090267112e+39 2.47295525368e-16
NIIb 6586.69531669 4.36990852161e-15 5.53353450775 2.57427080134e+40 7.41886576104e-16
SIIb 6734.10531669 1.25342148139e-15 5.53353450775 7.38378459262e+39 2.12795431905e-16
SIIa 6719.71531669 1.33401863855e-15 5.53353450775 7.85857464219e+39 2.26478544189e-16
H_alpha_B 6566.03931669 4.83056497428e-14 68.3479980471 5.6912781152e+41 5.97961753041e-16
H_alpha_N 6566.03931669 3.92425150256e-15 5.53353450775 2.31173856619e+40 6.66226649048e-16
H_alpha 6566.03931669 5.22299012454e-14 10.9497680664 9.23043874266e+41 1.26418840209e-15
    
"""
