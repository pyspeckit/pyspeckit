

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

    import readers,fitters,plotters,baseline

    def __init__(self,filename):
        if ".fits" in filename:
            try: 
                self.data,self.error,self.xarr,self.header = self.readers.open_1d_fits(filename)
                self.readers.parse_header(self.header)
            except TypeError:
                # do something else?
                pass
        elif ".txt" in filename:
            try:
                self.xarr,self.data,self.error = self.readers.open_1d_txt(filename)
            except:
                # try to use readcol to parse it as a wavelength / data / error file?
                # I think this is a typical format?
                pass


        self.specfit = fitters.Specfit(self)
        self.plotter = plotters.Plotter(self)
        self.baseline = baseline.Baseline(self)


