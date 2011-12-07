"""
This wrapper is intended to replicate the inputs of agpy's showspec.splat_1d code.
It is uglier than the spectrum class.
"""

import pyspeckit

def splat_1d(filename=None,vmin=None,vmax=None,button=None,dobaseline=False,
        exclude=None,smooth=None,order=1,savepre=None,vcrop=True,
        vconv=None,vpars=None,spec=None,xtora=None,ytodec=None,
        specname=None,quiet=True,specnum=0,errspecnum=None,wcstype='',
        offset=0.0, continuum=0.0, annotatebaseline=False, plotspectrum=True,
        smoothto=None, xunits=None, units=None, conversion_factor=None,
        smoothtype='gaussian',convmode='same',maskspecnum=None,
        fignum=1, axis=None, autorefresh=False, title=None, color=None,
        label=None,clear=False, negamp=None, plotscale=1.0, voff=0.0, **kwargs):

    sp = pyspeckit.Spectrum(filename, specnum=specnum, errspecnum=errspecnum, wcstype=wcstype, **kwargs)
    sp.plotter.xmin = vmin
    sp.plotter.xmax = vmax
    sp.plotter.offset = offset
    if type(continuum) is str:
        sp.plotter.offset += sp.header.get(continuum)
    else:
        sp.plotter.offset += continuum
    if vcrop and vmin and vmax:
        sp.crop(vmin,vmax)

    if button != None:
        print "button keyword doesn't do anything"
    if maskspecnum:
        print "maskspecnum doesn't do anything [not implemented]"

    if specname:
        sp.plotter.title = specname
    if title:
        sp.plotter.title = title

    if smoothto:
        smooth = abs(smoothto/sp.xarr.cdelt())
    if smooth:
        sp.smooth(smooth, smoothtype=smoothtype, convmode=convmode)

    if xunits:
        sp.xarr.convert_to_unit(xunits)
    if units:
        sp.units = units

    if dobaseline:
        sp.baseline(order=order, exclude=exclude, annotate=annotatebaseline)

    if plotspectrum:
        if color is None:
            sp.plotter(figure=fignum, axis=axis, autorefresh=autorefresh, clear=clear, silent=quiet, reset_ylimits=True, title=title, plotscale=plotscale, xoffset=voff, label=label)
        else:
            sp.plotter(figure=fignum, axis=axis, autorefresh=autorefresh, color=color, clear=clear, silent=quiet, reset_ylimits=True, title=title, plotscale=plotscale, xoffset=voff, label=label)

    if savepre is not None:
        glon,glat = sp.header.get("GLON"),sp.header.get("GLAT")
        if glon is None or glat is None:
            raise ValueError("Header doesn't contain glon/glat.  Default savepre will not work.")
        if glat < 0: pm="" 
        else: pm = "+"
        savename = savepre + "G%07.3f%0s%07.3f_" % (glon,pm,glat) + sp.header.get('MOLECULE').replace(' ','') + sp.header.get('TRANSITI').replace(' ','')
        sp.plotter.savefig(savename+'.png')

    sp.fitspec = sp.specfit
    sp.spectrum = sp.data
    sp.axis = sp.plotter.axis
    sp.dv = sp.xarr.cdelt()
    sp.vind = sp.xarr
    sp.cube = sp.data
    sp.save = sp.write

    return sp
