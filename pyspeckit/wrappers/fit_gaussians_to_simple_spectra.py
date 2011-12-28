"""
======================
Simple Gaussian Fitter
======================
"""
import pyspeckit

def fit_gaussians_to_simple_spectra(filename, units='km/s', doplot=True,
        baseline=True, plotresiduals=False, figuresavename=None,
        croprange=None, savename=None, **kwargs):
    """
    As stated in the name title, will fit Gaussians to simple spectra!

    kwargs will be passed to specfit

    *figuresavename* [ None | string ] 
        After fitting, save the figure to this filename if specified

    *croprange* [ list of 2 floats ]
        Crop the spectrum to (min,max) in the specified units

    *savename* [ None | string ]
        After fitting, save the spectrum to this filename

    Note that this wrapper can be used from the command line:

        python fit_gaussians_to_simple_spectra.py spectrum.fits
    
    """

    # load a FITS-compliant spectrum
    spec = pyspeckit.Spectrum(filename)
    if spec.xarr.units != units:
        spec.xarr.frequency_to_velocity()
        spec.xarr.convert_to_unit(units)

    if croprange is not None and len(croprange) == 2:
        spec.crop(*croprange)

    if doplot:
        spec.plotter()
    if baseline:
        spec.baseline()

    spec.specfit(**kwargs)
    spec.specfit(guesses=spec.specfit.modelpars)

    if plotresiduals:
        spec.specfit.plotresiduals()

    if figuresavename is not None:
        spec.plotter.figure.savefig(figuresavename)

    if savename is not None:
        spec.write(savename)

    return spec

if __name__ == "__main__":

    import sys
    import optparse

    parser=optparse.OptionParser()
    parser.add_option("--units","--unit","-u",help="Units.  Default is km/s",default="km/s")
    parser.add_option("--plot","--doplot","-p",help="Plot spectrum.  Default is "+\
            "True, but will not save unless a savename is specified.  You can "+\
            "set the savename with this option.",default=True)
    parser.add_option("--baseline","-b",help="Baseline the spectrum?  Default is True",default=True)
    parser.add_option("--plotresiduals",help="Plot the residual spectrum?  Default is False",default=False)
    parser.add_option("--savename","-s",help="Savename for the .fits or .txt (or other format) file.",default=None)
    parser.add_option("--crop","-c",help="Crop the spectrum?  Syntax is [x1,x2]",default=None)
    parser.set_description("""Fit a gaussian to a spectrum, optionally after baselining.  Save a 
    figure and/or a binary image.  Extra arguments will be passed via kwargs to spectrum.specfit""")

    options,args = parser.parse_args()
    croprange = [int(i) for i in options.crop.strip('[]').split(',')]
    if len(args) > 1: kwargs = args[1:]

    fit_gaussians_to_simple_spectra(args[0],units=options.units,
            doplot=bool(options.plot), baseline=options.baseline,
            plotresiduals=options.plotresiduals, figuresavename=options.plot,
            croprange=croprange, savename=options.savename, **kwargs)
