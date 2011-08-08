"""
Wrapper to fit ammonia spectra.  Generates a reasonable guess at the position and velocity using a gaussian fit
"""

def fitnh3(spectrum, vrange=[-100,100], vrangeunits='km/s', quiet=False,
        Tex=20,Tkin=15,column=1e15,fortho=1.0): 

    if vrange:
        spectrum.xarr.convert_to_unit(vrangeunits)
        spectrum.crop(*vrange)

    spectrum.specfit(fittype='gaussian',negamp=False)
    ampguess,vguess,widthguess = spectrum.specfit.modelpars

    spectrum.specfit(fittype='ammonia',quiet=quiet,multifit=True,guesses=[Tex,Tkin,column,widthguess,vguess,fortho])

    return spectrum
