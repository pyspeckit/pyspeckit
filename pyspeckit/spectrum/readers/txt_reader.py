"""
==========================
PySpecKit ASCII Reader
==========================

Routines for reading in ASCII format spectra.  If astropy.io.ascii is not
installed, will use a very simple routine for reading in the data.

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
.. moduleauthor:: Jordan Mirocha <mirochaj@gmail.com>
"""
from __future__ import print_function
from astropy.io import ascii
import numpy as np

from ...specwarnings import warn
from . import readcol
from .. import units

def open_1d_txt(filename, xaxcol=0, datacol=1, errorcol=2,
                text_reader='simple', format=None, **kwargs):
    """
    Attempt to read a 1D spectrum from a text file assuming wavelength as the
    first column, data as the second, and (optionally) error as the third.

    Reading can be done either with astropy.io.ascii or a 'simple' reader.  If
    you have an IPAC, CDS, or formally formatted table, you'll want to use
    astropy.io.ascii and specify a format.

    If you have a simply formatted file of the form, e.g.
    # name name
    # unit unit
    data data
    data data

    kwargs are passed to astropy.io.ascii.read
    """
    if text_reader in ('simple','readcol'):
        if text_reader == 'simple':
            data, error, XAxis, T = simple_txt(filename, xaxcol=xaxcol,
                                               datacol=datacol,
                                               errorcol=errorcol, **kwargs)
        elif text_reader == 'readcol':
            Tlist = readcol.readcol(filename, twod=False, **kwargs)
            XAxis = units.SpectroscopicAxis(Tlist[xaxcol])
            data = Tlist[datacol]
            error = Tlist[errorcol]

            T = dummy_class()

            T.data = dummy_class()
            T.data.dtype = dummy_class()
            T.columns = {}
            T.columns[T.data.dtype.names[xaxcol]] = dummy_class()
            T.columns[T.data.dtype.names[datacol]] = dummy_class()


    elif text_reader in ('ascii', 'astropy', 'asciitable'):
        T = ascii.read(filename, format=format, **kwargs)

        xarr = T.data[T.data.dtype.names[xaxcol]]
        data = T.data[T.data.dtype.names[datacol]]
        if len(T.columns) > errorcol:
            error = T.data[T.data.dtype.names[errorcol]]
        else:
            # assume uniform, zero error
            error = data*0

        if 'xunits' in T.keywords:
            xunits = T.keywords['xunits']
        else:
            xunits = 'unknown'

        XAxis = units.SpectroscopicAxis(xarr,xunits)

    # Need this in Spectrum class to correctly parse header
    T.xaxcol = xaxcol
    T.datacol = datacol

    return data, error, XAxis, T

def simple_txt(filename, xaxcol=0, datacol=1, errorcol=2, skiplines=0, **kwargs):
    """
    Very simple method for reading columns from ASCII file.
    """

    with open(filename, 'r') as f:

        hdr = None
        colunits = []
        coldata = []
        for ii, line in enumerate(f):

            # Ignore blank lines
            if not line.strip():
                continue

            # Possibly read in header
            if line.split()[0][0] == '#':
                if (ii) == (0+skiplines):
                    hdr = line[1:].split()
                if (ii) == (1+skiplines):
                    colunits = line[1:].split()

                continue

            if ii < skiplines:
                continue

            coldata.append(line.split())

            for j, element in enumerate(coldata[-1]):
                try:
                    coldata[-1][j] = float(element)
                except ValueError:
                    coldata[-1][j] = str(element)

    coldata = list(zip(*coldata))

    if not colunits:
        colunits = ['unknown'] * len(coldata)
    if not hdr:
        hdr = ['unknown'] * len(coldata)

    # Prepare to return data
    data = coldata[datacol]
    xarr = coldata[xaxcol]
    if errorcol > len(coldata) - 1:
        error = np.array(data)*0
    else:
        error = coldata[errorcol]

    if len(error) != len(data):
        raise ValueError("Data and Error lengths do not match.")

    XAxis = units.SpectroscopicAxis(xarr, colunits[xaxcol])

    # Create atPy style Table instance
    T = dummy_class()
    T.data = dummy_class()
    T.data.dtype = dummy_class()
    T.data.dtype.names = hdr
    T.columns = {}
    T.columns[T.data.dtype.names[xaxcol]] = dummy_class()
    T.columns[T.data.dtype.names[xaxcol]].unit = colunits[xaxcol]
    T.columns[T.data.dtype.names[datacol]] = dummy_class()
    T.columns[T.data.dtype.names[datacol]].unit = colunits[datacol]

    return np.array(data), np.array(error), XAxis, T

class dummy_class:
    def __init__(self):
        pass


