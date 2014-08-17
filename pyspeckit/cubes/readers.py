"""
Author: Adam Ginsburg
Created: 3/17/2011
"""
import numpy as np
import spectral_cube
import spectral_cube.io.fits
from .. import spectrum
from ..spectrum import units
import operator
from astropy import log

def open_3d_fits(filename, scale_keyword=None, scale_action=operator.div,
                 **kwargs):
    """
    Grabs all the relevant pieces of a simple FITS-compliant 3d data cube

    Parameters
    ----------
    scale_keyword : str
        A string to use to scale the data using the action
        ``scale_action``

    """
    cube = spectral_cube.io.fits.load_fits_cube(filename)

    if scale_keyword is not None:
        scale_value = cube.header[scale_keyword]
        data = scale_action(cube.filled_data[:], scale_value)
    else:
        data = cube.filled_data

    xunit = cube.spectral_axis.unit.to_string().replace(" ","")

    vconv = spectral_cube.spectral_axis.determine_vconv_from_ctype(cube.wcs.wcs.ctype[2]).__name__[8:]
    XAxis = units.SpectroscopicAxis(cube.spectral_axis.value,
                                    unit=xunit,
                                    xtype=units.unit_type_dict[xunit],
                                    refX=cube.wcs.wcs.restfrq,
                                    refX_units='Hz',
                                    velocity_convention=vconv,
                                    )

    return data,XAxis,cube.header,cube
