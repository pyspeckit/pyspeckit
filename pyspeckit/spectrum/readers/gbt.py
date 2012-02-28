"""
GBTIDL SDFITS file
"""
import pyfits
import pyspeckit
import numpy as np
try:
    import coords
except ImportError:
    coords = None

def unique_targets(sdfitsfile):
    bintable = _get_bintable(sdfitsfile)
    return uniq(bintable.data['OBJECT'])

def count_integrations(sdfitsfile, target):
    """
    Return the number of integrations for a given target
    (uses one sampler; assumes same number for all samplers)
    """
    bintable = _get_bintable(sdfitsfile)

    whobject = bintable.data['OBJECT'] == target
    any_sampler = bintable.data['SAMPLER'][whobject][0]
    whsampler = bintable.data['SAMPLER'][whobject] == any_sampler

    return (whsampler).sum()


def list_targets(sdfitsfile, doprint=True):
    """
    List the targets, their location on the sky...
    """
    bintable = _get_bintable(sdfitsfile)

    # Things to include:
    #   Scan           Source      Vel    Proc Seqn   RestF nIF nInt nFd     Az    El
    #        1      G33.13-0.09     87.4     Nod    1  14.488   4   15   2  209.7  47.9
    strings = [ "\n%18s  %10s %10s %26s%8s %9s %9s %9s" % ("Object Name","RA","DEC","%12s%14s"%("RA","DEC"),"N(ptgs)","Exp.Time","requested","n(ints)") ]
    for objectname in uniq(bintable.data['OBJECT']):
        whobject = bintable.data['OBJECT'] == objectname
        RA,DEC = bintable.data['TRGTLONG'][whobject],bintable.data['TRGTLAT'][whobject]
        RADEC = zip(RA,DEC)
        midRA,midDEC = np.median(RA),np.median(DEC)
        npointings = len(set(RADEC))
        if coords is not None:
            sexagesimal = coords.Position((midRA,midDEC)).hmsdms()
        else:
            sexagesimal = ""
        firstsampler = bintable.data['SAMPLER'][whobject][0]
        whfirstsampler = bintable.data['SAMPLER'] == firstsampler
        exptime = bintable.data['EXPOSURE'][whobject*whfirstsampler].sum()
        duration = bintable.data['DURATION'][whobject*whfirstsampler].sum()
        n_ints = count_integrations(bintable, objectname)
        strings.append( "%18s  %10f %10f %26s%8i %9g %9g %9i" % (objectname,midRA,midDEC,sexagesimal, npointings, exptime, duration, n_ints) )

    if doprint:
        print "\n".join(strings)

    return strings

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
    sp.xarr.frame = 'topo'

    # HACK - temporary!
    # Convert xarr to LSR units
    sp.xarr.convert_to_unit('m/s')
    sp.xarr -= header['VFRAME']
    sp.xarr.convert_to_unit('Hz')

    return sp

def read_gbt_target(sdfitsfile, objectname, verbose=False):
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

def reduce_gbt_target(sdfitsfile, objectname, verbose=False):
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
    nodpairs = uniq([s[:-1].replace("ON","").replace("OFF","") for s in blocks])

    reduced_nods = {}

    # for each pair...
    for sampname in nodpairs:
        on1 = sampname+"ON1"
        off1 = sampname+"OFF1"
        on2 = sampname+"ON2"
        off2 = sampname+"OFF2"

        feednumber = blocks[on1].header.get('FEED')
        # don't need this - if feednumber is 1, ref feed should be 2, and vice-versa
        reference_feednumber = blocks[on1].header.get('SRFEED')
        
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
    elif isinstance (sdfitsfile, pyfits.FITS_rec):
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
        print caloff
        print caloff.slice(pct10,pct90,units='pixels')
        print calon
        print calon.slice(pct10,pct90,units='pixels')
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

def find_matched_freqs(reduced_blocks, debug=False):
    """
    Use frequency-matching to find which samplers observed the same parts of the spectrum

    *WARNING* These IF numbers don't match GBTIDL's!  I don't know how to get those to match up!
    """

    # IF order 
    sampler_numbers = [int(name[1:]) for name in reduced_blocks.keys()]
    sorted_pairs = sorted(zip(sampler_numbers,reduced_blocks.keys()))
    sorted_names = zip(*sorted_pairs)[1]

    # how many IFs?
    frequencies = dict((name,reduced_blocks[name].header.get('OBSFREQ')) for name in sorted_names)
    if debug: print frequencies
    round_frequencies = dict((name,
        round_to_resolution(reduced_blocks[name].header.get('OBSFREQ'),
                            reduced_blocks[name].header.get('FREQRES')))
                 for name in sorted_names)
    if debug: print round_frequencies
    unique_frequencies = uniq(round_frequencies.values()) # uniq is an order-preserving function
    if debug: print unique_frequencies
    nIFs = len(unique_frequencies)
    if debug:
        print frequencies.keys()
        print round_frequencies.keys()
        print unique_frequencies

    IF_dict = {}
    # most annoying part of this whole process... (besides the obvious 1+1=23452342345 error...)
    # ifnumbers are approximately defined by the order things are written.... but only very approximately
    for ifnum,freq in enumerate(unique_frequencies):
        IF_dict[ifnum] = [sampler for (sampler,rf) in round_frequencies.iteritems() if rf == freq]
        if debug: print ifnum,freq,IF_dict[ifnum]

    return IF_dict

def uniq(seq):
    """ from http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


def round_to_resolution(frequency, resolution):
    """
    kind of a hack, but round the frequency to the nearest integer multiple of the resolution,
    then multiply it back into frequency space
    """
    # doesn't really work... divisor = 10**np.floor(np.log10(resolution))
    divisor = resolution
    return (np.round(frequency/divisor)) * divisor

def find_pols(block):
    """
    Get a dictionary of the polarization for each sampler
    """
    pols = dict((name,polnum_to_pol[int(block[name].header.get('CRVAL4'))]) for name in block)
    return pols

def find_feeds(block):
    """
    Get a dictionary of the feed numbers for each sampler
    """
    feeds = dict((name,block[name].header.get('FEED')) for name in block)
    return feeds

def identify_samplers(block):
    """
    Identify each sampler with an IF number, a feed number, and a polarization
    """
    feeds = find_feeds(block)
    pols = find_pols(block)
    ifs = find_matched_freqs(block)

    IDs = dict(
            (name,
                {'pol': pols[name],
                 'feed': feeds[name],
                 'IF': [k for k,v in ifs.iteritems() if name in v][0]})
                for name in block.keys()
                )

    return IDs
    
def average_pols(block):
    """
    Average the polarizations for each feed in each IF
    """

    # will have names like "(IFnum)(feed)"
    averaged_pol_dict = {}

    ifdict = find_matched_freqs(block)
    feeddict = find_feeds(block)
    IDs = identify_samplers(block)

    for ifnum,ifsampler in ifdict.iteritems():
        for sampler,feednum in feeddict.iteritems():
            if sampler not in ifsampler:
                continue
            newname = "if%ifd%i" % (ifnum, feednum)
            matched_samplers = [sampler_name for (sampler_name,ID) in IDs.iteritems()
                    if ID['feed'] == feednum and ID['IF'] == ifnum]
            if len(matched_samplers) != 2:
                raise ValueError("Too few/many matches: %s" % matched_samplers)

            if newname not in averaged_pol_dict:
                average = (block[matched_samplers[0]] + block[matched_samplers[1]]) / 2.
                averaged_pol_dict[newname] = average

    return averaged_pol_dict

def average_IF(block, debug=False):
    """
    Average the polarizations for each feed in each IF
    """

    # will have names like "(IFnum)"
    averaged_dict = {}

    ifdict = find_matched_freqs(block, debug=debug)

    import operator
    for ifnum,ifsamplers in ifdict.iteritems():
        if debug: print "if%i: freq %g" % (ifnum, block[ifsamplers[0]].header['OBSFREQ'])
        averaged_dict["if%i" % ifnum] = reduce(operator.add,[block[name] for name in ifsamplers]) / len(ifsamplers)

    return averaged_dict


polnum_to_pol = {
    1: 'I',
    2: 'Q',
    3: 'U',
    4: 'V',
   -1: 'RR',
   -2: 'LL',
   -3: 'RL',
   -4: 'LR',
   -5: 'XX',
   -6: 'YY',
   -7: 'XY',
   -8: 'YX'}


class GBTSession(object):
    """
    A class wrapping all of the above features
    """
    def __init__(self, sdfitsfile):
        """
        Load an SDFITS file or a pre-loaded FITS file
        """
        self.bintable = _get_bintable(sdfitsfile)
        self.targets = dict((target,None) for target in unique_targets(self.bintable))
        self.print_header = "\n".join( ("Observer: " + self.bintable.data[0]['OBSERVER'],
            "Project: %s" % self.bintable.header['PROJID'],
            "Backend: %s" % self.bintable.header['BACKEND'],
            "Telescope: %s" % self.bintable.header['TELESCOP'],
            "Bandwidth: %s" % self.bintable.data[0]['BANDWID'],
            "Date: %s" % self.bintable.data[0]['DATE-OBS']) )

    def __repr__(self):
        self.instance_info = super(GBTSession,self).__repr__()
        if not hasattr(self,'StringDescription'):
            self.StringDescription = list_targets(self.bintable, doprint=False)
        return self.print_header+"\n"+"\n".join(self.StringDescription)
    
    def __str__(self):
        if not hasattr(self,'StringDescription'):
            self.StringDescription = list_targets(self.bintable, doprint=False)
        return self.print_header+"\n"+"\n".join(self.StringDescription)

    def load_target(self, target, **kwargs):
        """
        Load a Target...
        """
        self.targets[target] = GBTTarget(self, target)
        return self.targets[target]

    def reduce_target(self, target, **kwargs):
        """
        Reduce the data for a given object name
        """
        if self.targets[target] is None:
            self.load_target(target)

        self.targets[target].reduce(**kwargs)

        return self.targets[target]

    def reduce_all(self):
        """

        """
        for target in self.targets:
            self.reduce_target(target)

class GBTTarget(object):
    """
    A collection of ObsBlocks or Spectra
    """
    def __init__(self, Session, target, **kwargs):
        """
        Container for the individual scans of a target from a GBT session
        """
        self.name = target
        self.Session = Session
        self.blocks = read_gbt_target(Session.bintable, target, **kwargs)
        self.spectra = {}

    def __getitem__(self, ind):
        return self.spectra[ind]

    def __repr__(self):
        self.instance_info = super(GBTTarget,self).__repr__()
        self.StringDescription = [("Object %s with %i scan blocks and %i 'reduced' spectra" %
                (self.name,len(self.blocks),len(self.spectra)))]
        if hasattr(self,'reduced_scans'):
            self.StringDescription += ["%s" % ID for ID in identify_samplers(self.reduced_scans)]
        return "\n".join(self.StringDescription)

    def reduce(self, obstype='nod', **kwargs):
        """
        Reduce nodded observations (they should have been read in __init__)
        """
        self.reduced_scans = reduce_blocks(self.blocks, **kwargs)
        self.spectra.update(self.reduced_scans)

    def average_pols(self):
        if hasattr(self,'reduced_scans'):
            self.averaged_pols = average_pols(self.reduced_scans)
            self.spectra.update(self.averaged_pols)
        else:
            raise ValueError("Must reduce the scans before averaging pols")

    def average_IFs(self, **kwargs):
        if hasattr(self,'reduced_scans'):
            self.averaged_IFs = average_IF(self.reduced_scans, **kwargs)
            self.spectra.update(self.averaged_IFs)
        else:
            raise ValueError("Must reduce the scans before averaging IFs")

"""
TEST CODE
for name in reduced_nods:
    for num in '1','2':

name='A10'
num='1'
tsys = 10
while (tsys>5):
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
