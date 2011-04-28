"""
Author: Adam Ginsburg
Created: 3/17/2011
"""
# import parent package
import spectrum
# import local things
import mapplot
import readers

class Cube(spectrum.Spectrum):

    def __init__(self,filename, x0=0, y0=0, **kwargs):
        """
        Initialize the Cube.  Accepts files in the following formats:
            - .fits

        x0,y0 - initial spectrum to use (defaults to lower-left corner)
        """

        if ".fit" in filename: # allow .fit or .fits
            try: 
                self.cube,self.xarr,self.header,self.fitsfile = readers.open_3d_fits(filename, **kwargs)
                self.data = self.cube[:,y0,x0]
                self.error = None
                self.parse_header(self.header)
            except TypeError as inst:
                print "Failed to read fits file."
                print inst
        else:
            raise TypeError('Not a .fits cube - what type of file did you give it?')

        self.fileprefix = filename.rsplit('.', 1)[0]    # Everything prior to .fits or .txt
        self.plotter = spectrum.plotters.Plotter(self)
        self.specfit = spectrum.fitters.Specfit(self)
        self.baseline = spectrum.baseline.Baseline(self)
        # Initialize writers
        self.writer = {}
        for writer in spectrum.writers.writers: self.writer[writer] = spectrum.writers.writers[writer](self)
        self.mapplot = mapplot.MapPlotter(self)

    def plot_spectrum(self, x, y, **kwargs):
        """
        Fill the .data array with a real spectrum and plot it
        """

        self.data = self.cube[:,y,x]

        self.plotter(**kwargs)
