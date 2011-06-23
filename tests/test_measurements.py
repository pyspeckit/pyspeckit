# Test measurements class on an SDSS AGN spectrum

import spectrum
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
spec = spectrum.Spectrum('sample_sdss.txt')

# H-alpha
spec.specfit.selectregion(xmin = Halpha - 5, xmax = Halpha + 5)
ampHa = np.max(spec.data[spec.specfit.gx1:spec.specfit.gx2])

# SII
spec.specfit.selectregion(xmin = SIIa - 20, xmax = SIIb + 20)
smallamp = np.max(spec.data[spec.specfit.gx1:spec.specfit.gx2])    

spec.specfit(guesses = [smallamp, SIIa, narrow, smallamp, SIIb, narrow], 
    tied = ['', '', 'p[-1]', '', 'p[1] + {0}'.format(SIIb_off), ''], negamp = False, quiet = True, multifit = True)
    
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

spec.specfit(guesses = guesses, tied = tied, fixed = fixed, negamp = False,
    limitedmin = lmin, limitedmax = lmax, minpars = minp, maxpars = maxp)
    
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

spec.plotter(xmin = NIIa - 100, xmax = SIIb + 30)

spec.specfit(guesses = guesses, tied = tied, fixed = fixed, negamp = False,
    limitedmin = lmin, limitedmax = lmax, minpars = minp, maxpars = maxp)
    
spec.plotter.refresh()

spec.measure(z = 0.05, fluxnorm = 1e-17)

# Overplot positions of lines and annotate
y = spec.plotter.ymax * 0.85
for i, line in enumerate(spec.measurements.lines.keys()):
    x = spec.measurements.lines[line]['modelpars'][1]
    spec.plotter.axis.plot([x]*2, [spec.plotter.ymin, spec.plotter.ymax], ls = '--', color = 'k')
    try: spec.plotter.axis.annotate(spec.speclines.optical.lines[line][-1], 
        (x, y), rotation = 90, ha = 'right', va = 'center')
    except KeyError: pass

spec.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
spec.plotter.axis.set_ylabel(r'Flux $(10^{-17} \mathrm{erg/s/cm^2/\AA})$')

print "Line   Flux (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)"
for line in spec.measurements.lines.keys():
    print line, spec.measurements.lines[line]['flux'], spec.measurements.lines[line]['fwhm'], spec.measurements.lines[line]['lum']

raw_input('Done.')
    