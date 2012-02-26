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
                calOK = (calON - (onoff=='OFF'))
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

def list_targets(sdfitsfile):
    bintable = _get_bintable(sdfitsfile)
    print "\n".join(np.unique(bintable.data['OBJECT']))

def reduce_gbt_target(sdfitsfile, objectname, verbose=True):
    """
    Wrapper - read an SDFITS file, get an object, reduce it (assuming nodded) and return it
    """
    # for efficiency, this should be stored - how?
    blocks = read_gbt_target(sdfitsfile, objectname, verbose=verbose)

    reduced_nods = reduce_blocks(blocks)

    return reduced_nods

def reduce_blocks(blocks, verbose=False):
    """
    Do a nodded on/off observation given a dict of observation blocks as
    produced by read_gbt_target
    """

    # strip off trailing digit to define pairs
    nodpairs = set([s[:-1].replace("ON","").replace("OFF","") for s in blocks])

    reduced_nods = {}

    # for each pair...
    for sampname in nodpairs:
        on1 = sampname+"ON1"
        off1 = sampname+"OFF1"
        on2 = sampname+"ON2"
        off2 = sampname+"OFF2"

        feednumber = blocks[on1].header.get('FEED')
        
        on1avg = blocks[on1].average()
        off1avg = blocks[off1].average()
        on2avg = blocks[on2].average()
        off2avg = blocks[off2].average()

        # first find TSYS
        tsys1 = dcmeantsys(on1avg,off1avg,on1avg.header.get('TCAL'))
        tsys2 = dcmeantsys(on2avg,off2avg,on2avg.header.get('TCAL'))
        if verbose:
            print "Nod Pair %s (feed %i) has tsys1=%f tsys2=%f" % (sampname, feednumber, tsys1, tsys2)

        # then get the total power
        tp1 = totalpower(on1avg,off1avg)
        tp2 = totalpower(on2avg,off2avg)

        # then do the signal-reference bit
        if feednumber == 1:
            nod = sigref(tp1,tp2,tsys2)
        elif feednumber == 2:
            nod = sigref(tp2,tp1,tsys1)

        reduced_nods[sampname] = nod

    return reduced_nods

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

def dcmeantsys(calon,caloff,tcal,debug=False):
    """
    from GBTIDL's dcmeantsys.py
    ;  mean_tsys = tcal * mean(nocal) / (mean(withcal-nocal)) + tcal/2.0
    """

    # Use the inner 80% of data to calculate mean Tsys
    nchans = calon.data.shape[0]
    pct10 = nchans/10
    pct90 = nchans - pct10

    meanoff = np.mean(caloff.slice(pct10,pct90,units='pixels').data)
    meandiff = np.mean(calon.slice(pct10,pct90,units='pixels').data - 
                        caloff.slice(pct10,pct90,units='pixels').data)

    meanTsys = ( meanoff / meandiff * tcal + tcal/2.0 )
    if debug:
        print "pct10: %i  pct90: %i mean1: %f mean2: %f tcal: %f tsys: %f" % (pct10,pct90,meanoff,meandiff,tcal,meanTsys)

    return meanTsys

def totalpower(calon, caloff):
    """
    Do a total-power calibration of an on/off data set
    (see dototalpower.pro in GBTIDL)
    """

    if hasattr(calon,'average'):
        ON = calon.average()
    else:
        ON = calon

    if hasattr(caloff,'average'):
        OFF = caloff.average()
    else:
        OFF = caloff

    TP = (ON+OFF)/2.

    return TP

def sigref(nod1, nod2, tsys_nod2):
    """
    Signal-Reference ('nod') calibration
    ; ((dcsig-dcref)/dcref) * dcref.tsys 
    see GBTIDL's dosigref
    """

    return (nod1-nod2)/nod2*tsys_nod2

"""
TEST CODE
for name in reduced_nods:
    for num in '1','2':
        av1 = blocks[name+'ON'+num].average()
        av2 = blocks[name+'OFF'+num].average()
        tsys = dcmeantsys(av1, av2, blocks[name+'OFF'+num].header['TCAL'],debug=True)
        if tsys < 5:
            print "%s %s: %f" % (name,num,tsys),
            print av1,np.mean(av1.slice(409,3687,units='pixels').data)
            print av2,np.mean(av2.slice(409,3687,units='pixels').data)
            print av1-av2,np.mean(av1.slice(409,3687,units='pixels').data-av2.slice(409,3687,units='pixels').data)
            print blocks[name+'OFF'+num].header['TCAL']


"""
