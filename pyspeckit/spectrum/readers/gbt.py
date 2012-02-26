"""
GBTIDL SDFITS file
"""
import pyfits
import pyspeckit
import numpy as np

def read_gbt_scan(sdfitsfile, obsnumber=0):
    """
    Read a single scan from a GBTIDL SDFITS file
    """

    bintable = _get_bintable(sdfitsfile)

    data = bintable.data[obsnumber]['DATA']

    header = pyfits.Header()

    for par in bintable.data.dtype.names:
        if par not in ('DATA',):
            header.update(par[:8], bintable.data[obsnumber][par])
    header.update('CUNIT1','Hz')

    HDU = pyfits.PrimaryHDU(data=data,header=header)

    sp = pyspeckit.Spectrum(HDU,filetype='pyfits')

    return sp

def read_gbt_target(sdfitsfile, objectname, verbose=True):
    """
    Give an object name, get all observations of that object as an 'obsblock'
    """

    bintable = _get_bintable(sdfitsfile)

    whobject = bintable.data['OBJECT'] == objectname
    if verbose:
        print "Number of individual scans for Object %s: %i" % (objectname,whobject.sum())

    calON = bintable.data['CAL'] == 'T'

    n_nods = np.unique(bintable.data['PROCSIZE'])

    blocks = {}
    for sampler in np.unique(bintable.data[whobject]['SAMPLER']):
        whsampler = bintable.data['SAMPLER'] == sampler
        nods = np.unique(bintable.data['PROCSEQN'][whsampler*whobject])
        for nod in nods:
            whnod = bintable.data['PROCSEQN'] == nod
            for onoff in ('ON','OFF'):
                calOK = (calON - (onoff=='ON'))
                whOK = (whobject*whsampler*calOK*whnod)
                if verbose:
                    print "Number of spectra for sampler %s, nod %i, cal%s: %i" % (sampler,nod,onoff,whOK.sum())
                crvals = bintable.data[whOK]['CRVAL1']
                maxdiff = np.diff(crvals).max()
                freqres = np.max(bintable.data[whOK]['FREQRES'])
                if maxdiff < freqres:
                    splist = [read_gbt_scan(bintable,ii) for ii in np.where(whOK)[0]]
                    blocks[sampler+onoff+str(nod)] = pyspeckit.ObsBlock(splist,force=True)
                    blocks[sampler+onoff+str(nod)]._arithmetic_threshold = np.diff(blocks[sampler+onoff+str(nod)].xarr).min() / 5.
                else:
                    print "Maximum frequency difference > frequency resolution: %f > %f" % (maxdiff, freqres)

    return blocks

def reduce_gbt_target(sdfitsfile, objectname, verbose=True):
    """
    Do a nodded on/off observation...
    """
    pass

def _get_bintable(sdfitsfile):
    """
    Private function: given a filename, HDUlist, or bintable, return a bintable
    """

    if isinstance(sdfitsfile, pyfits.HDUList):
        bintable = sdfitsfile[1]
    elif isinstance(sdfitsfile, pyfits.BinTableHDU):
        bintable = sdfitsfile
    else:
        bintable = pyfits.open(sdfitsfile)[1]

    return bintable
