"""
Storage for radio spectral line information.
"""

import numpy as np
import readcol
import re
greekletter = re.compile("&([a-z]*);")

splat = readcol.readcol("splatalogue.csv",fsep=":",asStruct=True)

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

lines = dict([(LineName(species,qn),[freq,'GHz',True, LatexName(species,qn)]) for
    species,freq,qn in zip(splat.Species,splat.FreqGHz,splat.ResolvedQNs)])

