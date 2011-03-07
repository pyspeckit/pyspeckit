

class Spectrum(object):
    """
    The core class for the spectroscopic toolkit.  Contains the data and error
    arrays along with wavelength / frequency / velocity information in various
    formats.

    Will contain functions / classes to:
        -read and write FITS-compliant spectroscopic data sets
            -read fits binaries?
        -plot a spectrum
        -fit a spectrum both interactively and non-interactively
            -with gaussian, voigt, lorentzian, and multiple (gvl) profiles
        -fit a continuum "baseline" to a selected region with an n-order
          polynomial
        -perform fourier transforms and operations in fourier space on a spectrum
    """

    def __init__(self,filename):
        """
        Initialize the Spectrum.  Accepts files in the following formats:
            - .fits
            - .txt (not functional yet)

        Ideally, we'd like to do something where each independent type of data
        file has its own subclass for file reading, but we haven't figured out
        how to automatically determine which subclass to use yet.
        """

        import readers,fitters,plotters,baseline

        if ".fit" in filename: # allow .fit or .fits
            try: 
                self.data,self.error,self.xarr,self.header = readers.open_1d_fits(filename)
                self.parse_header(self.header)
            except TypeError as inst:
                print "Failed to read fits file."
                print inst
        elif ".txt" in filename:
            try:
                self.xarr,self.data,self.error,self.Table = readers.open_1d_txt(filename)
                self.parse_text_header(self.Table)
            except Exception as inst:
                print "Reading txt failed.",inst.args
                print inst
                # try to use readcol to parse it as a wavelength / data / error file?
                # I think this is a typical format?
                pass

        self.plotter = plotters.Plotter(self)
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)

    def parse_text_header(self,Table):
        """
        Grab relevant parameters from a table header (xaxis type, etc)
        """
        self.xtype = Table.data.dtype.names[0]
        self.xunits = Table.columns[self.xtype].unit
        self.ytype = Table.data.dtype.names[1]
        self.units = Table.columns[self.ytype].unit

    def parse_header(self,hdr,specname=None, wcstype=''):
        """
        Parse parameters from a .fits header into required spectrum structure
        parameters

        This should be moved to the FITSSpectrum subclass when that is available
        """

        self.xunits = hdr.get('CUNIT1'+wcstype)

        if hdr.get('BUNIT'):
            self.units = hdr.get('BUNIT').strip()
        else:
            self.units = 'undefined'

        if hdr.get('REFFREQ'+wcstype):
            self.reffreq = hdr.get('REFFREQ'+wcstype)
        else:
            self.reffreq = None

        if hdr.get('CTYPE1'+wcstype):
            self.xtype = hdr.get('CTYPE1'+wcstype)
        else:
            self.xtype = 'VLSR'

        if specname is not None:
            self.specname = specname
        elif hdr.get('OBJECT'):
            self.specname = hdr.get('OBJECT')

