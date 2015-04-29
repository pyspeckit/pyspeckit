import pyspeckit
from pyspeckit.spectrum import units
import numpy as np
import pytest
from astropy import units as u

convention = ['optical','radio','relativistic']

unit_type_dict = {
    'Hz' :'frequency', 'kHz':'frequency', 'MHz':'frequency', 'GHz':'frequency',
    'hz' :'frequency', 'khz':'frequency', 'mhz':'frequency', 'ghz':'frequency',
    'THz':'frequency', 
    'meters/second':'velocity', 'm/s':'velocity', 'kilometers/s':'velocity',
    'kilometers/second':'velocity',
    'centimeters/second':'velocity',
    'megameters/s':'velocity',
    'megameters/second':'velocity','Mm/s':'velocity',
    'km/s':'velocity', 'kms':'velocity', 'centimeters/s':'velocity',
    'kms-1':'velocity', 'km s-1':'velocity', 'm s-1':'velocity',  'ms-1':'velocity',
    'cm/s':'velocity', 'cms':'velocity', 
    'meters':'wavelength','m':'wavelength',
    'centimeters':'wavelength','cm':'wavelength',
    'millimeters':'wavelength','mm':'wavelength',
    'nanometers':'wavelength','nm':'wavelength',
    'micrometers':'wavelength','micron':'wavelength','microns':'wavelength','um':'wavelength',
    'kilometers':'wavelength','km':'wavelength',
    'megameters':'wavelength','Mm':'wavelength',
    'angstrom':'wavelength','angstroms':'wavelength','A':'wavelength',
    'unknown':'pixels',
    None: 'pixels',
}

params = [(a,b,c,d) for a in unit_type_dict if a not in ('unknown',None)
                  for b in unit_type_dict if b not in ('unknown',None)
                  for c in convention
                  for d in ['GHz','cm']]
params = params[:5]

@pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
def test_convert_to(unit_from, unit_to, convention, ref_unit):
    x = units.SpectroscopicAxis(np.arange(5), unit=u.angstrom, velocity_convention='optical', equivalencies=u.doppler_optical(3*u.AA))
    x.convert_to_unit(unit_to,convention=convention)


def test_equivalencies_1():
    x = units.SpectroscopicAxis(np.arange(5), unit=u.angstrom, velocity_convention='optical', equivalencies=u.doppler_optical(3*u.AA))
    assert x.equivalencies is not None # == u.doppler_optical(3*u.AA)

@pytest.mark.parametrize(('velocity_convention'), convention)
def test_equivalencies_2(velocity_convention):
    x = SpectroscopicAxis(np.arange(5), unit=u.angstrom, refX=3*u.AA, velocity_convention=velocity_convention)
    assert x.equivalencies == u.doppler_optical(3*u.AA)

@pytest.mark.parametrize(('velocity_convention'), convention)
def test_equivalencies_3(velocity_convention):
    x = SpectroscopicAxis(np.arange(5), unit=u.angstrom, refX=3, refX_unit='angstrom', velocity_convention=velocity_convention)
    assert x.equivalencies == u.doppler_optical(3*u.AA)

def test_initialize_units():
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=u.dimensionless_unscaled)
    assert xarr.unit == u.dimensionless_unscaled
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=u.m)
    assert xarr.unit == u.m

@pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
def test_convert_units(unit_from,unit_to,convention,ref_unit):
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_unit=ref_unit)
    xarr.convert_to_unit(unit_to,convention=convention)
    assert xarr.unit == unit_to

def test_convert_units2():
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit='Hz')

    velocity_arr = xarr.as_unit('km/s', velocity_convention='optical', refX=5*u.GHz)

    xarr2 = units.SpectroscopicAxis(np.linspace(1,10,10),unit='Hz',
                                    velocity_convention='optical',
                                    refX=5*u.GHz)
    velocity_arr2 = xarr2.to(u.km/u.s)

    assert np.all(velocity_arr.value==velocity_arr2.value)



@pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
def test_convert_back(unit_from, unit_to,convention,ref_unit):
    if unit_from in ('cms','cm/s') or unit_to in ('cms','cm/s'):
        xvals = np.linspace(1000,10000,10)
        threshold = 1e-6
    else:
        xvals = np.linspace(1,10,10)
        threshold = 1e-15
    # all conversions include a * or / by speedoflight_ms
    threshold = np.spacing(units.speedoflight_ms) * 100
    if 'megameter' in unit_from or 'Mm' in unit_from:
        threshold *= 10
    xarr = units.SpectroscopicAxis(xvals,unit=unit_from,refX=5,refX_unit=ref_unit,velocity_convention=convention)
    xarr.convert_to_unit(unit_to,convention=convention)
    xarr.convert_to_unit(unit_from,convention=convention)
    assert all(np.abs((xarr.value - xvals)/xvals) < threshold)
    assert xarr.unit == unit_from


if __name__=="__main__":
    unit_from='GHz'
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_unit='GHz')
    xarr2=xarr.as_unit('km/s',quiet=False,debug=True)
    xarr3=xarr2.as_unit('GHz',quiet=False,debug=True)
        
