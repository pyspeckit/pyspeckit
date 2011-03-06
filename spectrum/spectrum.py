

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

        import readers,fitters,plotters,baseline

        if ".fits" in filename:
            try: 
                self.data,self.error,self.xarr,self.header = readers.open_1d_fits(filename)
                self.parse_header(self.header)
            except TypeError:
                # do something else?
                pass
        elif ".txt" in filename:
            try:
                self.xarr,self.data,self.error = readers.open_1d_txt(filename)
            except:
                # try to use readcol to parse it as a wavelength / data / error file?
                # I think this is a typical format?
                pass

        self.plotter = plotters.Plotter(self)
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)

    def parse_header(self,hdr,specname=None, wcstype=''):
        """
        Parse parameters from a .fits header into required spectrum structure
        parameters
        """
        if hdr.get('CUNIT1'+wcstype) in ['m/s','M/S']:
            self.xunits = 'km/s' # change to km/s because you're converting units
        else:
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

