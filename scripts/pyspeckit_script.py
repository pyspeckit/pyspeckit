#/bin/env ipython -i --matplotlib
"""
pyspeckit command line startup script
"""
from __future__ import print_function
import sys
# remove script file's parent directory from path
# (otherwise, can't import pyspeckit)
#sys.path.pop(0)
from pyspeckit.spectrum.classes import Spectrum, Spectra
from pyspeckit.cubes.SpectralCube import Cube, CubeStack
from pyspeckit import wrappers as pw
import optparse
import os
import re

if __name__ == "__main__":

    import matplotlib
    import itertools
    import pylab

    parser=optparse.OptionParser()
    parser.add_option("--verbose","-v",help="Be loud? Default True",default=False,action='store_true')
    parser.add_option("--debug","-d",help="Debug mode.  Default False",default=False,action='store_true')
    parser.add_option("--doplot","-p",help="Plot? Default True",default=True)
    parser.add_option("--fitgaussian",help="Fit a gaussian?",default=False,action='store_true')
    parser.add_option("--fitnh3",help="Fit NH3?",default=False,action='store_true')
    parser.add_option("--threed",'--3d','--cube',help="Data cube?",default=False,action='store_true')
    parser.add_option("--filetype",help="File type to use.", default=None)
    parser.add_option("--smooth",help="Smooth the spectrum (by how much)?",default=False)
    parser.add_option("--wcstype",help="What wcstype to use?  Can be a list: A,B,C,T,V where elements correspond to input spectra",default=None)
    parser.add_option("--specnum",help="What specnum?",default=0)
    parser.add_option("--hdu",help="What HDU number?",default=None)
    parser.add_option("--unmerged",help="[tspec only] Is the tspec file NOT merged?",default=False,action='store_true')

    options,args = parser.parse_args()

    verbose = options.verbose

    if verbose:
        print("Args: ",args)
        print("Options: ",options)

    if options.debug:
        print("DEBUG MODE.  Using a different excepthook.")
        def info(type, value, tb):
            if hasattr(sys, 'ps1') or not sys.stderr.isatty():
                # we are in interactive mode or we don't have a tty-like
                # device, so we call the default hook
                sys.__excepthook__(type, value, tb)
            else:
                import traceback, pdb
                # we are NOT in interactive mode, print the exception...
                traceback.print_exception(type, value, tb)
                print()
                # ...then start the debugger in post-mortem mode.
                pdb.pm()
        
        sys.excepthook = info

    specnum = int(options.specnum)

    if options.wcstype:
        if "," in options.wcstype:
            wcstype = options.wcstype.split(",")
        else:
            wcstype = options.wcstype
    else:
        wcstype = ''

    # specify kwargs before passing both for brevity and to allow
    # for some kwargs to not be specified at all
    kwargs = {'specnum':specnum,
            'wcstype':wcstype,
            'verbose':verbose,
            'filetype':options.filetype}
    if options.hdu is not None:
        kwargs['hdu'] = int(options.hdu)
    
    if len(args) > 1:
        if options.threed:
            cubelist = [Cube(fname) for fname in args]
            splist = cubelist
            cube = CubeStack(cubelist)
            options.doplot = False
        elif len(wcstype) == len(args):
            splist = [Spectrum(a, **kwargs) for a, w in
                zip(args, wcstype)]
            sp = Spectra(splist)
        else:
            splist = [Spectrum(a,**kwargs) for a in args]
            sp = Spectra(splist)
        linestyles = itertools.cycle(["steps-mid","steps-mid--"])
        colors = itertools.cycle(matplotlib.cm.spectral(pylab.linspace(0,1,len(splist))))
    else:
        if len(wcstype)==1:
            sp = Spectrum(*args,**kwargs)
        else:
            if options.threed:
                cube = Cube(*args)
                options.doplot = False
            else:
                sp = Spectrum(*args,**kwargs)
    if options.smooth > 0: sp.smooth(float(options.smooth))
    if options.doplot: sp.plotter()

    if options.fitgaussian:
        sp.specfit()

    if options.fitnh3:
        pw.fitnh3.fitnh3(sp)

    import IPython
    IPython.embed()
