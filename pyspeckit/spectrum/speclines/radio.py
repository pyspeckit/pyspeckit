"""
Storage for radio spectral line information.
"""
from __future__ import print_function

import numpy as np
import copy
import re
import warnings
from astropy.table import Table
from .. import units

greekletterlist = ["Alpha","Nu","Beta","Xi","Gamma","Omicron","Delta","Pi","Epsilon","Rho","Zeta","Sigma","Eta","Tau","Theta","Upsilon","Iota","Phi","Kappa","Chi","Lambda","Psi","Mu","Omega"]
greekletterlist += ["alpha","nu","beta","xi","gamma","omicron","delta","pi","epsilon","rho","zeta","sigma","eta","tau","theta","upsilon","iota","phi","kappa","chi","lambda","psi","mu","omega"]

greekletter = re.compile("&([A-Za-z][a-z]*);?")
sub = re.compile("<sub>(.*)</sub>")
sup = re.compile("<sup>(.*)</sup>")
junk = re.compile("&([^a-zA-Z]*);?")
R = re.compile('\\W')  # find and remove all non-alphanumeric characters

import os
selfpath = os.path.split(os.path.abspath(__file__))[0]

def LineName(species,qn):
    if greekletter.search(species) is not None:
        #spn = greekletter.sub(greekletter.search(species).groups()[0][0],species)
        name = greekletter.sub(greekletter.search(species).groups()[0][0],qn.replace("(","").replace(")",""))
    else:
        name = species+qn
    if sub.search(name) is not None:
        name = sub.sub("("+sub.search(name).groups()[0]+")",name)
    if sup.search(name) is not None:
        name = sup.sub("("+sup.search(name).groups()[0]+")",name)
    return name.replace(" ","")

def LatexName(species,qn):
    if greekletter.search(species) is not None:
        spn = greekletter.sub(r"$\\%s$" % greekletter.search(species).groups()[0],species)
        name = qn.replace("(","$_{").replace(")","}$").replace("&","$\\").replace(";","$")
        if greekletter.search(name) is not None:
            name = greekletter.sub("\\%s" % greekletter.search(name).groups()[0],name) 
    else:
        spn = species
        name = spn+"$(%s)$" % (qn.replace("(","_{").replace(")","}"))
        if greekletter.search(name) is not None:
            name = greekletter.sub(r"\\%s" % greekletter.search(name).groups()[0],name) 
    if junk.search(name) is not None:
        name = junk.sub("",name)
    if sub.search(name) is not None:
        name = sub.sub("_{"+sub.search(name).groups()[0]+"}",name)
    if sup.search(name) is not None:
        name = sup.sub("^{"+sup.search(name).groups()[0]+"}",name)
    if "&" in name:
        raise ValueError("But... that's... not... possible!")
    return name.replace(" ","")

#def radio_lines(minwav=0.0,maxwav=1.0):
#    lines = dict([(LineName(species,qn),[freq,'GHz',True, LatexName(species,qn)]) for
#        species,freq,qn in zip(splat.Species,splat.FreqGHz,splat.ResolvedQNs)])
#
#    #splat = query_splatalogue(minwav,maxwav)


class radio_lines(object):
    def __init__(self, spectrum, voff=0.0, webquery=False, **kwargs):
        """
        Initialize the radio lines class
        Requires a spectrum object

        (webquery is deprecated and ignored; use astroquery.splatalogue to
        query splatalogue directly)
        """
        self.Spectrum = spectrum

        self.minfreq_GHz = self.Spectrum.xarr.umin(unit='GHz')*(1+voff/units.speedoflight_kms)
        self.maxfreq_GHz = self.Spectrum.xarr.umax(unit='GHz')*(1+voff/units.speedoflight_kms)

        self.table = get_splat_table(webquery=webquery, 
                minwav=units.speedoflight_ms/self.maxfreq_GHz/1e9,
                maxwav=units.speedoflight_ms/self.minfreq_GHz/1e9,
                waveunits='m')

        self.voff = voff
        self._lines = []
        self._linenames = []

    def show(self, voff=None, ymax_scale=0.8, userecommended=True,
            maxupperstateenergy=None, minupperstateenergy=None, color='r',
            verbose=False, force=False, regexp=None, **kwargs):
        """
        Display vertical lines (using 'vlines') at the position of each
        discovered line
        """
        ymin = self.Spectrum.plotter.ymin
        ymax = self.Spectrum.plotter.ymax
        if ymax < 0:
            ymax_scale = 1.0/ymax_scale

        self.hide()
        voff = self.voff if voff is None else voff

        frequency = np.asarray(self.table['frequency'])

        mask = np.ones(len(self.table),dtype='bool')
        mask *= ((frequency > self.Spectrum.xarr.umin(unit='GHz')*(1+voff/units.speedoflight_kms)) *
                (frequency < self.Spectrum.xarr.umax(unit='GHz')*(1+voff/units.speedoflight_kms)))

        if userecommended:
            if 'frequencyrecommended' in self.table.colnames:
                mask *= np.asarray(self.table['frequencyrecommended'], dtype='bool')
            else:
                warnings.warn("Table has no 'frequencyrecommended' column; "
                              "ignoring userecommended.")
        if maxupperstateenergy is not None:
            mask *= (np.asarray(self.table['upperstateenergyK']) < maxupperstateenergy)
        if minupperstateenergy is not None:
            mask *= (np.asarray(self.table['upperstateenergyK']) > minupperstateenergy)
        if regexp is not None:
            import re
            reg = re.compile(regexp)
            mask *= np.array([reg.search(sn) is not None for sn in self.table['Species']])

        freqoff = voff * 1e3 / units.speedoflight_ms * frequency[mask]

        if mask.sum() > 0:
            if mask.sum() > 30 and not force:
                print("WARNING: show() will plot %i lines!  Use force=True if you want this to happen anyway." % (mask.sum()))
                return
            if verbose: print("Labeled %i lines." % mask.sum())
            xcoords = self.Spectrum.xarr.x_to_coord(frequency[mask]-freqoff, 'GHz')
            # x_to_coord may return a Quantity; matplotlib needs plain floats
            xcoords = np.atleast_1d(getattr(xcoords, 'value', xcoords))
            self._lines += [self.Spectrum.plotter.axis.vlines(
                    xcoords, ymin, ymax, color=color, **kwargs)]
            self._linenames += [self.Spectrum.plotter.axis.text(
                    XCOORD,
                    ymax*ymax_scale,
                    NAME,
                    rotation='vertical', color=color)
                for XCOORD,NAME in zip(xcoords,self.table['LatexName'][mask])]
        else:
            if verbose: print("No lines found in range [%g, %g] GHz" % (self.minfreq_GHz, self.maxfreq_GHz))

        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

    def hide(self):
        """
        Remove all annotations and lines
        """
        for text in self._linenames:
            if text in self.Spectrum.plotter.axis.texts:
                text.remove()
        self._linenames = []
        for line in self._lines:
            if line in self.Spectrum.plotter.axis.collections:
                line.remove()
        self._lines = []

        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()
 
def get_splat_table(webquery=False, savename=None, filename=None, **kwargs):
    """
    Load a splatalogue line table (colon-delimited splatalogue.csv export)
    into an `astropy.table.Table`.

    Parameters
    ----------
    webquery : bool
        Deprecated & ignored.  The legacy SLAP web-query interface has been
        removed; use `astroquery.splatalogue` to query splatalogue directly.
    savename : str or None
        If specified, write the resulting table to this filename.
    filename : str or None
        Path to a colon-delimited splatalogue CSV export.  Defaults to
        ``splatalogue.csv`` in the ``pyspeckit/spectrum/speclines`` directory.
    """
    if webquery:
        warnings.warn("pyspeckit's SLAP-based splatalogue web query has been "
                      "removed; use astroquery.splatalogue instead.  Falling "
                      "back to the local splatalogue.csv table.",
                      DeprecationWarning)

    if filename is None:
        filename = os.path.join(selfpath, "splatalogue.csv")
    if not os.path.exists(filename):
        raise IOError("Splatalogue table file {0} not found.  Download a "
                      "colon-delimited CSV export from splatalogue.online, or "
                      "use astroquery.splatalogue to build a line table."
                      .format(filename))

    splat = Table.read(filename, format='ascii', delimiter=':', guess=False)

    for cn in splat.colnames:
        if cn != R.sub('',cn):
            splat.rename_column(cn,R.sub('',cn))

    # entries with no computed (or measured) frequency are read as masked (or
    # blank-string) values; coerce them to float columns filled with -999
    for freqcol in ('FreqGHz', 'MeasFreqGHz'):
        if freqcol in splat.colnames:
            col = splat[freqcol]
            if col.dtype.kind in 'SU':
                data = np.array([x if x not in ('', None) else '-999'
                                 for x in col], dtype='float')
                splat[freqcol] = data
            elif hasattr(col, 'filled'):
                splat[freqcol] = col.filled(-999)

    latex_names = [LatexName(species,qn) for species,qn in zip(splat['Species'],splat['ResolvedQNs'])]
    line_names  = [LineName(species,qn) for species,qn in zip(splat['Species'],splat['ResolvedQNs'])]
    splat['LineName'] = line_names
    splat['LatexName'] = latex_names
    splat['frequency'] = splat['FreqGHz']

    if savename is not None:
        splat.write(savename)

    return splat
