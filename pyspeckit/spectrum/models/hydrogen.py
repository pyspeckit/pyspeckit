# -*- coding: latin-1 -*-
"""
Hydrogen Models
---------------

Hydrogen in HII regions is typically assumed to follow Case B recombination
theory.

The values for the Case B recombination coefficients are given by `Hummer &
Storey (1987)
<http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1987MNRAS.224..801H&db_key=AST>`_.
They are also computed in `Hummer (1994)
<http://adsabs.harvard.edu/abs/1994MNRAS.268..109H>`_ and tabulated at a `wiki
<http://wiki.hmet.net/index.php/HII_Case_B_Recombination_Coefficients>`_.  I
had to OCR and pull out by hand some of the coefficients.

.. Hummer & Storey 1987: Recombination-line intensities for hydrogenic ions.  I - Case B calculations for H I and He II


However, Storey & Hummer 1995 might be a better resource now?
https://cdsarc.unistra.fr/viz-bin/cat/VI/64


Module API
^^^^^^^^^^

"""
import numpy as np
from six import iteritems
from . import model
from .. import units

# log(temperature), alphaBalphaB = [[1.0,9.283],
alphaB=[[1.2,8.823],
    [1.4,8.361],
    [1.6,7.898],
    [1.8,7.435],
    [2.0,6.973],
    [2.2,6.512],
    [2.4,6.054],
    [2.6,5.599],
    [2.8,5.147],
    [3.0,4.700],
    [3.2,4.258],
    [3.4,3.823],
    [3.6,3.397],
    [3.8,2.983],
    [4.0,2.584],
    [4.2,2.204],
    [4.4,1.847],
    [4.6,1.520],
    [4.8,1.226],
    [5.0,0.9696],
    [5.2,0.7514],
    [5.4,0.5710],
    [5.6,0.4257],
    [5.8,0.3117],
    [6.0,0.2244],
    [6.2,0.1590],
    [6.4,0.1110],
    [6.6,0.07642],
    [6.8,0.05199],
    [7.0,0.03498]]

# beta is 946 / H 03b2 / oct 1662


# the command to make this:
# %s/\v.*\n?([0-9]\.[0-9]+).*jH(.+\n?|[^β]\n?[0-9]\.[0-9]+) ?([0-9]\.[0-9]+) ([0-9]\.[0-9]+) ([0-9]\.[0-9]+)/[\1, \3, \4, \5]\r/
# from this:
# http://www.scribd.com/doc/70881169/Physics-of-the-Interstellar-and-Intergalactic-Medium
# with some by-hand editing too
# vacuum wavelengths computed by:
#  [(ii,1.0/10973731.6 * 1./(1/2.**2 - 1/float(ii)**2) * 1e6 * (1+1/1836.15266) ) for ii in range(10,20)]

# Table 14.2 in Draine 2001
temperatures = {5000:0,10000:1,20000:2}
table14dot2 = {
#Balmer-line intensities relative to Hβ 0.48627 µm
'balmer':[["a" ,0.65646, 3.03, 2.86, 2.74],
          ["b" ,0.48627, 1., 1., 1.,],
          ["g" ,0.43418, 0.459, 0.469, 0.475],
          ["d" ,0.41030, 0.252, 0.259, 0.264],
          ["e" ,0.39713, 0.154, 0.159, 0.163],
          ["8" ,0.38902, 0.102, 0.105, 0.106],
          ["9" ,0.38365, 0.0711, 0.0732, 0.0746],
          ["10",0.37990, 0.0517, 0.0531, 0.0540],
          ["11",0.377174, 0.0000, 0.0398, 0.0000],
          ["12",0.375125, 0.0000, 0.0306, 0.0000],
          ["13",0.373547, 0.0000, 0.0240, 0.0000],
          ["14",0.372303, 0.0000, 0.0193, 0.0000],
          ["15",0.371306, 0.0000, 0.0157, 0.0000],
          ["16",0.370494, 0.0000, 0.0130, 0.0000],
          ["17",0.369824, 0.0000, 0.0108, 0.0000],
          ["18",0.369264, 0.0000, 0.00918, 0.0000],
          ["19",0.368792, 0.0000, 0.00785, 0.0000],
          ["20",0.368389,0.0000,0.006790,0.0000] ,
          ["21",0.368044,0.0000,0.005920,0.0000] ,
          ["22",0.367745,0.0000,0.005210,0.0000] ,
          ["23",0.367484,0.0000,0.004610,0.0000] ,
          ["24",0.367256,0.0000,0.004120,0.0000] ,
          ["25",0.367054,0.0000,0.003700,0.0000] ,
          ["26",0.366876,0.0000,0.003350,0.0000] ,
          ["27",0.366718,0.0000,0.003040,0.0000] ,
          ["28",0.366576,0.0000,0.002780,0.0000] ,
          ["29",0.366448,0.0000,0.002550,0.0000] ,
          ["30",0.366333,0.0000,0.002350,0.0000] ,
          ["31",0.366230,0.0000,0.001650,0.0000] ,
          ["32",0.366136,0.0000,0.001220,0.0000] ,
          ["33",0.366050,0.0000,0.000919,0.0000] ,
          ["34",0.365972,0.0000,0.000708,0.0000]],

#Paschen (n→3) line intensities relative to corresponding Balmer lines
'paschen':[["a", 1.8756, 0.405, 0.336, 0.283],
           ["b", 1.2821, 0.399, 0.347, 0.305],
           ["g", 1.0941, 0.391, 0.348, 0.311],
           ["d", 1.0052, 0.386, 0.348, 0.314],
           ["e", 0.95487, 0.382, 0.348, 0.316],
           ["9", 0.92317, 0.380, 0.347, 0.317],
           ["10",0.90175, 0.380, 0.347, 0.317]],


#Brackett (n→4) line intensities relative to corresponding Balmer lines
'brackett':[["a", 4.0523, 0.223, 0.169, 0.131],
            ["b", 2.6259, 0.219, 0.174, 0.141],
            ["g", 2.1661, 0.212, 0.174, 0.144],
            ["d", 1.9451, 0.208, 0.173, 0.145],
            ["e", 1.8179, 0.204, 0.173, 0.146],
            ["10",1.7367, 0.202, 0.172, 0.146]],

#pfund (n→5) line intensities relative to corresponding Balmer lines
'pfund':[["a", 7.4599, 0.134, 0.0969, 0.0719],
         ["b", 4.6538, 0.134, 0.101, 0.0774] ,
         ["g", 3.7406, 0.130, 0.101, 0.0790] ,
         ["d", 3.2970, 0.127, 0.100, 0.0797] ,
         ["e",3.0392, 0.125, 0.0997, 0.0801],
         ["11",2.873004,0.00, 0.084171, 0.0000],
         ["12",2.758276,0.00, 0.098693, 0.0000],
         ["13",2.675139,0.00, 0.098333, 0.0000],
         ["14",2.612655,0.00, 0.097927, 0.0000],
         ["15",2.564334,0.00, 0.097452, 0.0000],
         ["16",2.526098,0.00, 0.096923, 0.0000],
         ["17",2.495261,0.00, 0.097222, 0.0000],
         ["18",2.469994,0.00, 0.096187, 0.0000],
         ["19",2.449007,0.00, 0.095541, 0.0000],
         ["20",2.431369,0.00, 0.094551, 0.0000],
         ["21",2.416392,0.00, 0.093750, 0.0000],
         ["22",2.403559,0.00, 0.092514, 0.0000],
         ["23",2.392474,0.00, 0.091540, 0.0000],
         ["24",2.382830,0.00, 0.075728, 0.0000],
         ["25",2.374384,0.00, 0.083784, 0.0000],
         ["26",2.366943,0.00, 0.087761, 0.0000],
         ["27",2.360353,0.00, 0.086513, 0.0000],
         ["28",2.354488,0.00, 0.083094, 0.0000],
         ["29",2.349243,0.00, 0.083922, 0.0000],
         ["30",2.344534,0.00, 0.082979, 0.0000],
         ["31",2.340290,0.00, 0.077576, 0.0000],
         ["32",2.336451,0.00, 0.073197, 0.0000],
         ["33",2.332966,0.00, 0.071055, 0.0000],
         ["34",2.329793,0.00, 0.069774, 0.0000]],


#Humphreys (n→6) line intensities relative to corresponding Balmer lines
'humphreys':[["a", 6.12372, 0.0855, 0.0601, 0.0435],
             ["b", 7.5026, 0.0867, 0.0632, 0.0471],
             ["g", 5.9083, 0.0850, 0.0634, 0.0481],
             ["d", 5.1287, 0.0833, 0.0632, 0.0486],
             ["e", 4.673, 0.0833, 0.0632, 0.0486]] # intensities are copied from delta
}

r_to_hbeta = {}
wavelength = {}

for line in table14dot2['balmer']:
    r_to_hbeta["balmer"+line[0]] = np.array(line[2:])
    wavelength["balmer"+line[0]] = line[1]

for series,values in iteritems(table14dot2):
    if series == 'balmer':
        continue
    for line in values:
        r_to_hbeta[series+line[0]] = np.array(line[2:]) * r_to_hbeta['balmer'+line[0]]
        wavelength[series+line[0]] = line[1]
        
name_num_map = {1: 'Lyman',
                2: 'Balmer',
                3: 'Paschen',
                4: 'Brackett',
                5: 'Pfund',
                6: 'Humphreys',
               }
greek_num_map = {1: 'alpha',
                 2: 'beta',
                 3: 'gamma',
                 4: 'delta',
                 5: 'epsilon',
                }

def retrieve_storey1995(temperature=10000, case='b'):
    """
    Retrieve the Hummer and Storey tables from Vizier:
    https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/64
    
    based on the temperature and recombination case.  
    
    
    The data are parsed into a dictionary using the function
    parse_storey1995 below
    """
    
    import requests
    import gzip
    import io
    
    # https://cdsarc.cds.unistra.fr/ftp/VI/64/r1b0010.d.gz
    temval = int(temperature / 100)
    # r: prefix. ?
    # 1: hydrogen
    # b: case b
    assert case in ('a', 'b')
    filename = f'r1{case}{temval:04d}.d.gz'
    resp = requests.get(f'https://cdsarc.cds.unistra.fr/ftp/VI/64/{filename}')
    zz = gzip.open(io.BytesIO(resp.content))
    data = zz.read().decode()
    
    return parse_storey1995(data.split("\n"))
    
    
def parse_storey1995(data, e_or_r='E'):
    """
    Returns a dictionary with keys like:"
    (100.0, 1, 1000.0, 'B', 2, 124)
    
    where the keys are
    
    dens, Z, temp, case, nmin, nc
    
    dens is the density in cm^-3, Z is the number of electrons,
    temp is the temperature in K, case is case A or B recombination
    (Case B, no transitions into n=1, no Lyman transitions, are included)
    nmin and nc - I don't know what these are
    
    Each entry contains a mapping from [NU,NL] -> emissivity, where
    the emissivity is in erg/s/cm^3 and NU, NL are the upper and lower
    electronic levels
    
    """
    entries = {}
    next_defines_entry = False
    coefficient_reading = False
    r_or_e = 'E' if e_or_r == 'R' else 'R'
    for row in data:
        if next_defines_entry:
            #print('entry definition:', row)
            dens, Z, temp, case, nmin, nc = row.split()
            dens, temp = map(float, (dens, temp))
            Z, nmin, nc = map(int, (Z, nmin, nc))
            entries[(dens, Z, temp, case, nmin, nc)] = {}
            next_defines_entry = False
            continue

        if row.startswith(f' {e_or_r}_NU'):
            #print('row definition:', row)
            NE = float(row.split("NE=")[-1].strip().split()[0])
            TE = float(row.split("TE=")[-1].strip().split()[0])
            E_NU = int(row.split(f"{e_or_r}_NU=")[-1].strip().split()[0])
            coefficient_reading = True
            continue

        if row.startswith(f" {r_or_e}_NU"):
            coefficient_reading = False
            continue

        if coefficient_reading:
            sp = row.strip().split()
            for nl, rate in zip(sp[:-1:2], sp[1::2]):
                entries[(dens, Z, temp, case, nmin, nc)][(E_NU, int(nl))] = float(rate)

        if row.startswith(' DENS'):
            next_defines_entry = True
            coefficient_reading = False
            continue

    return entries    


# not used right now, but it could be so I'm listing it here
def rrl(n,dn=1,amu=1.007825,Z=1):
    """
    compute Radio Recomb Line freqs in GHz
    from Brown, Lockman & Knapp ARAA 1978 16 445

    UPDATED:
        Gordon & Sorochenko 2009, eqn A6

    Parameters
    ----------
    n : int
        The number of the lower level of the recombination line (H1a is Lyman
        alpha, for example)
    dn : int
        The delta-N of the transition.  alpha=1, beta=2, etc.
    amu : float
        The mass of the central atom
    Z : int
        The ionization parameter for the atom.  Z=1 is neutral, Z=2 is singly
        ionized, etc.  For hydrogen, only z=1 makes sense, since ionized
        hydrogen has no electrons and therefore cannot have recombination
        lines.

    Returns
    -------
    frequency in GHz
    """
    nu = 3.289842e6*(1-5.48593e-4/amu)*(1/float(n)**2 - 1/float(n+dn)**2)
    nu = 3.28984196e6 * Z**2 * (amu-5.48579911e-4*(Z+1))/(amu-5.48579911e-4*Z) * (1/float(n)**2 - 1/float(n+dn)**2)
    return nu

def find_lines(xarr):
    """
    Given a :class:`pyspeckit.units.SpectrosopicAxis` instance, finds all the
    lines that are in bounds.  Returns a list of line names.
    """
    return [linename for (linename, lam) in iteritems(wavelength)
            if xarr.as_unit('micron').in_range(lam)]

def hydrogen_fitter(sp, temperature=10000, tiedwidth=False):
    """
    Generate a set of parameters identifying the hydrogen lines in your
    spectrum.  These come in groups of 3 assuming you're fitting a gaussian to
    each.  You can tie the widths or choose not to.

    *temperature* [ 5000, 10000, 20000 ]
        The case B coefficients are computed for 3 temperatures

    *tiedwidth* [ bool ]
        Should the widths be tied?

    Returns a list of `tied` and `guesses` in the xarr's units
    """

    tnum = temperatures[temperature]

    lines = find_lines(sp.xarr)
    reference_line = lines[0]
    subtied = [
        ['p[0]*%f' % (r_to_hbeta[line][tnum]/r_to_hbeta[reference_line][tnum]),
         'p[1]+%f' % (wavelength[line]-wavelength[reference_line]),
         'p[2]' if tiedwidth else '']
        for line in lines[1:]]
    tied = ['','',''] + [t for sublist in subtied for t in sublist]

    dx = np.median(sp.xarr.dxarr)

    subguesses = [[sp.data[sp.xarr.as_unit('micron').x_to_pix(wavelength[line])],
        sp.xarr.x_to_coord(wavelength[line], 'micron'),
        dx*2.0] for line in lines]
    guesses = [g for sublist in subguesses for g in sublist]

    return tied, guesses

def hydrogen_model(xarr, amplitude=1.0, width=0.0, velocity=0.0, a_k=0.0,
                   temperature=10000):
    """
    Generate a set of parameters identifying the hydrogen lines in your
    spectrum.  These come in groups of 3 assuming you're fitting a gaussian to
    each.  You can tie the widths or choose not to.

    Parameters
    ----------
    sp : pyspeckit.Spectrum
        The spectrum to fit
    temperature : [ 5000, 10000, 20000 ]
        The case B coefficients are computed for 3 temperatures
    a_k : float
        The K-band extinction normalized to 2.2 microns.  Simple
        exponential.
    width : float
        Line width in km/s
    velocity : float
        Line center in km/s
    amplitude : float
        arbitrary amplitude of the first line (all other lines
        will be scaled to this value)

    Returns
    -------
    np.ndarray with same shape as sp.xarr

    """

    tnum = temperatures[temperature]

    lines = find_lines(xarr)
    reference_line = lines[0]

    model = np.array(xarr * 0)
    xarr = xarr.as_unit('micron')


    lw = width / units.speedoflight_kms * wavelength[reference_line]
    center = wavelength[reference_line]
    model += amplitude * np.exp(-(xarr-center)**2 / (2.0*lw)**2)

    for line in lines[1:]:
        relamp = (r_to_hbeta[line][tnum]/r_to_hbeta[reference_line][tnum])
        lw = width / units.speedoflight_kms * wavelength[line]
        center = wavelength[line]

        model += amplitude * relamp * np.exp(-(xarr-center)**2 / (2.0*lw)**2)

    if a_k > 0:
        # alpha=1.8 comes from Martin & Whittet 1990.  alpha=1.75 from Rieke and Lebofsky 1985
        Al = a_k * (xarr/2.2)**(-1.75)
        model *= np.exp(-Al)

    return model

def add_to_registry(sp):
    """
    Add the Hydrogen model to the Spectrum's fitter registry
    """
    # can't have absorption in recombination case
    extincted_hydrogen_emission = model.SpectralModel(hydrogen_model, 4,
                                                      shortvarnames=('A','\\sigma','\\Delta x','A_K'),
                                                      parnames=['amplitude','width','velocity','extinction'],
                                                      parlimited=[(True,False),(True,False),(False,False),
                                                                  (True,False)],
                                                      parlimits=[(0,0), (0,0),
                                                                 (0,0), (0,0)],
            fitunit='micron')

    sp.Registry.add_fitter('hydrogen', extincted_hydrogen_emission,
                           extincted_hydrogen_emission.npars)
