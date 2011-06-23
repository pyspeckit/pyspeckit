"""
Storage for radio spectral line information.
"""

import numpy as np
import re

try:
    import atpy
    atpyOK = True
except ImportError:
    import readcol
    atpyOK = False

try:
    import query_splatalogue
    webOK = True
except ImportError:
    webOK = False

greekletter = re.compile("&([a-z]*);")
R = re.compile('\W')  # find and remove all non-alphanumeric characters

import os
selfpath = os.path.split(os.path.abspath(__file__))[0]

def LineName(species,qn):
    if greekletter.search(species) is not None:
        #spn = greekletter.sub(greekletter.search(species).groups()[0][0],species)
        name = greekletter.sub(greekletter.search(species).groups()[0][0],qn.replace("(","").replace(")",""))
    else:
        name = species+qn
    return name

def LatexName(species,qn):
    if greekletter.search(species) is not None:
        spn = greekletter.sub(r"$\\%s$" % greekletter.search(species).groups()[0],species)
        name = qn.replace("(","_{").replace(")","}").replace("&","$\\").replace(";","")
    else:
        spn = species
        name = spn+"$(%s)$" % (qn.replace("(","_{").replace(")","}"))
    return name

#def radio_lines(minwav=0.0,maxwav=1.0):
#    lines = dict([(LineName(species,qn),[freq,'GHz',True, LatexName(species,qn)]) for
#        species,freq,qn in zip(splat.Species,splat.FreqGHz,splat.ResolvedQNs)])
#
#    #splat = query_splatalogue(minwav,maxwav)
 
def radio_lines(webquery=False,**kwargs):
    if webquery and webOK:
        splat = query_splatalogue.query_splatalogue(**kwargs)
    elif atpyOK:
        splat = atpy.Table(selfpath+"/splatalogue.csv",type='ascii',delimiter=':')
    else:
        splat = readcol.readcol(selfpath+"/splatalogue.csv",fsep=":",asStruct=True)

    for cn in splat.columns:
        if cn != R.sub('',cn):
            splat.rename_column(cn,R.sub('',cn))

    latex_names = [LatexName(species,qn) for species,qn in zip(splat.Species,splat.ResolvedQNs)]
    line_names  = [LineName(species,qn) for species,qn in zip(splat.Species,splat.ResolvedQNs)]
    if hasattr(splat,'add_column'):
        splat.add_column("LineName",line_names)
        splat.add_column("LatexName",latex_names)

    return splat
