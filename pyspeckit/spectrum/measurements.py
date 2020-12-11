from __future__ import print_function
import numpy as np
from six.moves import xrange
import itertools
from . import cosmology
from collections import OrderedDict

"""
To test:

import spectrum
spec = spectrum.Spectrum('sample_sdss.txt')
spec.plotter(xmin = 6400, xmax = 6800)
spec.specfit(guesses = [20, 6718.29, 5, 100, 6564.614, 20, 50, 6585.27, 20, 20, 6732.67, 5, 50, 6549.86, 5])
spec.measure()
"""

cm_per_mpc = 3.08568e+24

class Measurements(object):
    def __init__(self, Spectrum, z=None, d=None, fluxnorm=None,
                 miscline=None, misctol=10., ignore=None, derive=True, debug=False,
                 restframe=False, ptol=2, sort=False):
        """
        This can be called after a fit is run.  It will inherit the specfit
        object and derive as much as it can from modelpars.  Just do:
        spec.measure(z, xunits, fluxnorm)

        Notes: If z (redshift) or d (distance) are present, we can compute
        integrated line luminosities rather than just fluxes.  Provide distance
        in cm.

        Only works with Gaussians.  To generalize:
            1. make sure we manipulate modelpars correctly, i.e. read in
            entries corresponding to wavelength/frequency/whatever correctly.

        Parameters
        ----------
        z: float or None
            redshift
        d: float or None
            distance in cm (used for luminosities)
        fluxnorm: bool
            Normalize the fluxes?
        miscline: dictionary
            miscline = [{'name': H_alpha', 'wavelength': 6565}]
        misctol: tolerance (in Angstroms) for identifying an unmatched line
            to the line(s) we specify in miscline dictionary.
        sort: bool
            Sort the entries in order of observed wavelength (or velocity or
            frequency)

        """
        self.debug = debug

        self.restframe = restframe

        # Inherit specfit object
        self.specfit = Spectrum.specfit
        self.speclines = Spectrum.speclines

        # Bit of a hack - help identifying unmatched lines
        self.miscline = miscline
        self.misctol = misctol

        # Flux units in case we are interested in line luminosities or just having real flux units
        if fluxnorm is not None:
            self.fluxnorm = fluxnorm
        else:
            self.fluxnorm = 1

        # This is where we'll keep our results
        self.lines = OrderedDict()

        # Read in observed wavelengths
        tmp1 = np.reshape(self.specfit.modelpars, (int(len(self.specfit.modelpars) / 3), 3))
        tmp2 = np.reshape(self.specfit.modelerrs, (int(len(self.specfit.modelerrs) / 3), 3))

        if ignore is not None:
            tmp1 = np.delete(tmp1, ignore, 0)
            tmp2 = np.delete(tmp2, ignore, 0)

        # each tmp1 contains amplitude,wavelength,width
        # (Assumes gaussians)
        wavelengths = tmp1[:,1]
        
        # sort by wavelength
        if sort:
            order = np.argsort(wavelengths)
            self.obspos = wavelengths[order]
        else:
            order = np.arange(wavelengths.size)
            self.obspos = wavelengths

        self.Nlines = wavelengths.size

        # Read in modelpars and modelerrs, re-organize so they are 2D arrays sorted by ascending wavelength
        self.modelpars = np.zeros_like(tmp1)
        self.modelerrs = np.zeros_like(tmp2)
        for i, element in enumerate(order):
            self.modelpars[i] = tmp1[element]
            self.modelerrs[i] = tmp2[element]

        # Read in appropriate list of reference wavelengths/frequencies/whatever
        self.reflines = self.speclines.optical.get_optical_lines()
        self.refpos = self.reflines['xarr']
        self.refname = self.reflines['name']

        # Redshift reference lines if restframe = True
        if self.restframe and z is not None:
            self.refpos *= (1.0 + z)

        # If distance or redshift has been provided, we can compute luminosities from fluxes
        if d is not None:
            self.d = d
        else:
            self.d = None
        if z is not None:
            self.cosmology = cosmology.Cosmology()
            self.d = self.cosmology.LuminosityDistance(z) * cm_per_mpc

        self.unmatched = self.identify_by_position(ptol=ptol)

        #if np.sum(unmatched) >= 2:
        #    self.identify_by_spacing(unmatched)
        if derive:
            self.derive()

    def identify_by_position(self, ptol):
        """
        Match observed lines to nearest reference line.  Don't use spacing at all.

        ptol = tolerance (in angstroms) to accept positional match
        """

        if not hasattr(self, 'lines'):
            self.lines = OrderedDict()

        # Fill lines dictionary
        unmatched = np.zeros_like(self.obspos)
        for i, pos in enumerate(self.obspos):

            # Check miscline directory for match
            matched = False
            if self.miscline is not None:

                for line in self.miscline:
                    if abs(pos - line['wavelength']) > ptol:
                        continue

                    matched = True
                    name = line['name']
                    break

            if not matched:
                diff = np.abs(pos - self.refpos)
                loc = np.argmin(diff)

                if diff[loc] <= ptol:
                    matched = True

                    name = self.refname[loc]
                    if name in self.lines.keys():
                        name += '_1'

                        num = int(name[-1])
                        while name in self.lines.keys():
                            num += 1
                            name = '%s_%i' % (self.refname[loc], num)

            if matched:
                self.lines[name] = {}
                self.lines[name]['modelpars'] = list(self.modelpars[i])
                self.lines[name]['modelerrs'] = list(self.modelerrs[i])
            else:
                name = 'unknown_1'
                num = 1
                while name in self.lines.keys():
                    num += 1
                    name = 'unknown_%i' % num

                self.lines[name] = {}
                self.lines[name]['modelpars'] = list(self.modelpars[i])
                self.lines[name]['modelerrs'] = list(self.modelerrs[i])
                unmatched[i] = 1

        return unmatched

    def identify_by_spacing(self):
        """
        Determine identity of lines in self.modelpars.  Fill entries of self.lines dictionary.

        Note: This method will be infinitely slow for more than 10 or so lines.
        """

        if self.unmatched is None:
            self.unmatched = np.ones_like(self.obspos)

        # Remove lines that were already identified
        obspos = self.obspos[self.unmatched == 1]

        # Spacing between observed lines (odiff) and reference lines (rdiff)
        self.odiff = np.abs(np.diff(obspos))
        self.rdiff = np.abs(np.diff(self.refpos))

        # Don't try to identify lines with separations smaller than the smallest
        # separation in our reference library
        self.rdmin = 0.99 * min(self.rdiff)

        # If lines have multiple components (i.e. spacing much closer than ref lines),
        # delete them from ID list.
        if np.any(self.odiff) < self.rdmin:
            where = np.ravel(np.argwhere(self.odiff < self.rdmin))
            odiff = np.delete(self.odiff, where)
            multi = True
        else:
            where = 0
            odiff = self.odiff
            multi = False

        refpos = self.refpos
        refname = self.refname

        # Don't include elements of reference array that are far away from the observed lines (speeds things up)
        condition = (refpos >= 0.99 * min(self.obspos)) & (refpos <= 1.01 * max(self.obspos))
        refpos = refpos[condition]
        refname = refname[condition]

        if len(refpos) == 0:
            print('WARNING: No reference lines in this wavelength regime.')
        elif len(refpos) < self.Nlines:
            print('WARNING: More observed lines than reference lines in this band.')

        # Construct all possible (N-element) combos of reference lines
        combos = itertools.combinations(refpos, min(self.Nlines, len(refpos)))

        # List to keep track of line identification.  Each entry is (cost, (line1, line2, line3,...))
        self.IDresults = []
        for i, combo in enumerate(combos):
            rdiff = np.diff(combo)

            if len(odiff) == len(rdiff):
                result = (np.sum(np.abs(odiff - rdiff)), combo)
                self.IDresults.append(result)
            else: # If more/less observed lines than reference lines, try excluding observed lines one at a time
                if len(odiff) > len(rdiff):
                    subcombos = itertools.combinations(odiff, len(rdiff))
                    for subcombo in subcombos:
                        result = (np.sum(np.abs(subcombo - rdiff)), combo)
                        self.IDresults.append(result)
                else:
                    subcombos = itertools.combinations(rdiff, len(odiff))
                    for subcombo in subcombos:
                        result = (np.sum(np.abs(odiff - subcombo)), combo)
                        self.IDresults.append(result)

        # Pick best solution
        best = np.argmin(zip(*self.IDresults)[0])  # Location of best solution
        ALLloc = []                                # x-values of best fit lines in reference dictionary

        # Determine indices of matched reference lines
        for element in self.IDresults[best][1]:
            ALLloc.append(np.argmin(np.abs(refpos - element)))

        # Fill lines dictionary
        for i, element in enumerate(ALLloc):
            line = refname[element]
            self.lines[line] = {}
            loc = np.argmin(np.abs(self.obspos - refpos[element]))
            self.lines[line]['modelpars'] = list(self.modelpars[loc])
            self.lines[line]['modelerrs'] = list(self.modelerrs[loc])

        # Track down odd lines (i.e. broad components of lines already identified)
        # This won't yet work for lines that are truly unidentified
        if len(ALLloc) < self.Nlines:

            # Figure out which modelpars/errs that belong to lines that were already identified
            mpars = self.modelpars.copy()
            merrs = self.modelerrs.copy()
            for line in self.lines:
                wavelengths = zip(*mpars)[1]
                i = np.argmin(np.abs(zip(*mpars)[1] - self.lines[line]['modelpars'][1]))
                mpars = np.delete(mpars, i, 0)
                merrs = np.delete(merrs, i, 0)

            # Loop over unmatched modelpars/errs, find name of unmatched line, extend corresponding dict entry
            if self.miscline is None:
                for i, x in enumerate(zip(*mpars)[1]):
                    self.lines['unknown%i' % i] = {}
                    self.lines['unknown%i' % i]['modelpars'] = mpars[i]
                    self.lines['unknown%i' % i]['modelerrs'] = merrs[i]

            # If we've know a-priori which lines the unmatched lines are likely to be, use that information
            else:
                print(self.miscline)
                for i, miscline in enumerate(self.miscline):
                    for j, x in enumerate(zip(*mpars)[1]):
                        if abs(x - miscline['wavelength']) < self.misctol:
                            name = miscline['name']
                        else:
                            name = 'unknown%i' % j

                        self.lines[name] = {}
                        self.lines[name]['modelpars'] = mpars[j]
                        self.lines[name]['modelerrs'] = merrs[j]

        self.separate()

    def derive(self):
        """
        Calculate luminosity and FWHM for all spectral lines.
        """

        for line in self.lines.keys():
            if self.debug:
                print("Computing parameters for line %s" % line)

            self.lines[line]['fwhm'] = self.compute_fwhm(self.lines[line]['modelpars'])
            self.lines[line]['flux'] = self.compute_flux(self.lines[line]['modelpars'])
            self.lines[line]['amp'] = self.compute_amplitude(self.lines[line]['modelpars'])
            self.lines[line]['pos'] = self.lines[line]['modelpars'][1]

            if self.d is not None:
                self.lines[line]['lum'] = self.compute_luminosity(self.lines[line]['modelpars'])

    def separate(self):
        """
        For multicomponent lines, separate into broad and narrow components (assume only one of components is narrow).
        """

        for key in self.lines.keys():
            modpars = self.lines[key]['modelpars']
            moderrs = self.lines[key]['modelerrs']
            if len(modpars) > 3:
                modpars2d = np.reshape(modpars, (len(modpars) / 3, 3))
                moderrs2d = np.reshape(moderrs, (len(moderrs) / 3, 3))
                sigma = zip(*modpars2d)[2]
                minsigma = min(np.abs(sigma))
                i_narrow = list(np.abs(sigma)).index(minsigma)
            else:
                continue

            self.lines["{0}_N".format(key)] = {}
            self.lines["{0}_N".format(key)]['modelpars'] = []
            self.lines["{0}_N".format(key)]['modelerrs'] = []
            self.lines["{0}_B".format(key)] = {}
            self.lines["{0}_B".format(key)]['modelpars'] = []
            self.lines["{0}_B".format(key)]['modelerrs'] = []

            for i, arr in enumerate(modpars2d):
                if i == i_narrow:
                    self.lines["{0}_N".format(key)]['modelpars'] = arr
                    self.lines["{0}_N".format(key)]['modelerrs'] = moderrs2d[i]
                else:
                    self.lines["{0}_B".format(key)]['modelpars'].extend(arr)
                    self.lines["{0}_B".format(key)]['modelerrs'].extend(moderrs2d[i])

    def compute_flux(self, pars):
        """
        Calculate integrated flux of emission line.  Works for multi-component fits too.  Unnormalized.
        """

        flux = 0
        niter = (len(pars) / 3)
        assert niter == int(niter)
        for i in xrange(int(niter)):
            flux += np.sqrt(2. * np.pi) * pars[3 * i] * abs(pars[2 + 3 * i])

        return flux * self.fluxnorm

    def compute_amplitude(self, pars):
        """
        Calculate amplitude of emission line.  Should be easy - add multiple components if they exist.
        Currently assumes multiple components have the same centroid.
        """

        amp = 0
        niter = (len(pars) / 3)
        for i in xrange(int(niter)):
            amp += pars[3 * i]
        return amp * self.fluxnorm

    def compute_luminosity(self, pars):
        """
        Determine luminosity of line (need distance and flux units).
        """

        lum = 0
        niter = (len(pars) / 3)
        for i in xrange(int(niter)):
            lum += self.compute_flux(pars) * 4. * np.pi * self.d**2
        return lum

    def compute_fwhm(self, pars):
        """
        Determine full-width at half maximum for multi-component fit numerically, or analytically if line
        has only a single component.  Uses bisection technique for the former with absolute tolerance of 1e-4.
        """

        if len(pars) == 3:
            return 2. * np.sqrt(2. * np.log(2.)) * abs(pars[2])
        else:
            atol = 1e-4
            niter = (len(pars) / 3)
            pars2d = np.reshape(pars, (int(niter), 3))
            start = list(zip(*pars2d))[1][0]                    # start at central wavelength of first component

            # If the centroids are exactly the same for all components, we know the peak, and peak position
            if np.allclose(zip(*pars2d)[1], atol):
                fmax = np.sum(zip(*pars2d)[0])

            # Otherwise, we have to figure out where the multicomponent peak is
            else:
                f = lambda x: self.specfit.fitter.slope(x)
                xfmax = self.bisection(f, start)
                fmax = self.specfit.fitter.n_modelfunc(pars)(np.array([xfmax, xfmax]))[0]

            hmax = 0.5 * fmax

            # current height relative to half max - we want to minimize this function.  Could be asymmetric.
            f = lambda x: self.specfit.fitter.n_modelfunc(pars)(np.array([x])) - hmax
            xhmax1 = self.bisection(f, start)
            xhmax2 = self.bisection(f, start + (start - xhmax1))

            return abs(xhmax2 - xhmax1)

    def bisection(self, f, x_guess):
        """
        Find root of function using bisection method.  Absolute tolerance of 1e-4 is being used.
        """

        x1, x2 = self.bracket_root(f, x_guess)

        # Narrow bracketed range with bisection until tolerance is met
        while abs(x2 - x1) > 1e-4:
            midpt = np.mean([x1, x2])
            fmid = f(midpt)

            if np.sign(fmid) < 0: x1 = midpt
            else: x2 = midpt

            if fmid == 0.0: break

        return x2

    def bracket_root(self, f, x_guess, atol = 1e-4):
        """
        Bracket root by finding points where function goes from positive to negative.
        """

        f1 = f(x_guess)
        f2 = f(x_guess + 1)
        df = f2 - f1

        # Determine whether increasing or decreasing x_guess will lead us to zero
        if (f1 > 0 and df < 0) or (f1 < 0 and df > 0): sign = 1
        else: sign = -1

        # Find root bracketing points
        xpre = x_guess
        xnow = x_guess + sign
        fpre = f1
        fnow = f(xnow)
        while (np.sign(fnow) == np.sign(fpre)):
            xpre = xnow
            xnow += sign * 0.1
            fpre = f(xpre)
            fnow = f(xnow)

        x1 = min(xnow, xpre)
        x2 = max(xnow, xpre)

        if not np.all([np.sign(fpre), np.sign(fnow)]):
            x1 -= 1e-4
            x2 += 1e-4

        return x1, x2

    def to_tex(self):
        """
        Write out fit results to tex format.
        """

        pass



