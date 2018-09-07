from __future__ import print_function
import os
try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False

# rewrite this garbage

class write_txt(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum

    def write_data(self, overwrite = True):
        """
        Write all fit information to an ASCII file.
        """

        fn = "{0}_fit.dat".format(self.Spectrum.fileprefix)
        if not overwrite:
            i = 1
            while os.path.exists(fn):
                fn = "{0}_fit({1}).dat".format(self.Spectrum.fileprefix, i)
                i += 1

        with open(fn, 'w') as f:

            # Print header
            print("# Column 1: {0}".format("x-values"), file=f)
            print("# Column 2: {0}".format("model spectrum"), file=f)
            for i, element in enumerate(self.Spectrum.specfit.modelcomponents):
                print("# Column {0}: model spectrum component {1}".format(i + 3, i + 1), file=f)
            print("# Column {0}: residuals".format(i + 4), file=f)
            print("", file=f)

            components = zip(*self.Spectrum.specfit.modelcomponents)
            for i, element in enumerate(self.Spectrum.specfit.model):
                line = "{0:10}{1:10}".format(self.Spectrum.xarr[self.Spectrum.specfit.gx1:self.Spectrum.specfit.gx2][i],
                    round(self.Spectrum.specfit.model[i], 5))
                for j, component in enumerate(components[i]): line += "{0:10}".format(round(component, 5))
                line += "{0:10}".format(round(self.Spectrum.specfit.residuals[i], 5))

                print(line, file=f)

            print("", file=f)
