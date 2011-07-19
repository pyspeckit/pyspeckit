try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False

# rewrite this garbage

class write_txt(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum

    def write_data(self, clobber = True):
        """
        Write all fit information to an ASCII file.
        """
        
        fn = "{0}_fit.dat".format(self.Spectrum.fileprefix)
        if not clobber:
            i = 1
            while os.path.exists(fn):
                fn = "{0}_fit({1}).dat".format(self.Spectrum.fileprefix, i)
                i += 1
                
        f = open(fn, 'w')
        
        # Print header
        print >> f, "# Column 1: {0}".format("x-values")
        print >> f, "# Column 2: {0}".format("model spectrum")
        for i, element in enumerate(self.Spectrum.specfit.modelcomponents):
            print >> f, "# Column {0}: model spectrum component {1}".format(i + 3, i + 1)        
        print >> f, "# Column {0}: residuals".format(i + 4)
        print >> f, ""    
                
        components = zip(*self.Spectrum.specfit.modelcomponents)        
        for i, element in enumerate(self.Spectrum.specfit.model):
            line = "{0:10}{1:10}".format(self.Spectrum.xarr[self.Spectrum.specfit.gx1:self.Spectrum.specfit.gx2][i], 
                round(self.Spectrum.specfit.model[i], 5))
            for j, component in enumerate(components[i]): line += "{0:10}".format(round(component, 5))    
            line += "{0:10}".format(round(self.Spectrum.specfit.residuals[i], 5))       
                
            print >> f, line
            
        print >> f, ""
            
        f.close()

