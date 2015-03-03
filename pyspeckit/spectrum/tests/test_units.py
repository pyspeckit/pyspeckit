import pyspeckit
from pyspeckit.spectrum import units
import numpy as np
import pytest
from astropy import units as u

convention = ['optical','radio','relativistic']
params = [(a,b,c,d) for a in units.unit_type_dict if a not in ('unknown',None)
                  for b in units.unit_type_dict if b not in ('unknown',None)
                  for c in convention
                  for d in ['GHz','cm']]
params = params[:5]


def test_equivalencies_1():
    x = SpectroscopicAxis(np.arange(5), unit=u.angstrom, velocity_convention='optical', equivalencies=u.doppler_optical(3*u.AA))
    assert x.equivalencies == u.doppler_optical(3*u.AA)

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
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_units=ref_unit,xtype=units.unit_type_dict[unit_from])
    xarr.convert_to_unit(unit_to,convention=convention)
    assert xarr.units == unit_to

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
    xarr = units.SpectroscopicAxis(xvals,unit=unit_from,refX=5,refX_units=ref_unit,xtype=units.unit_type_dict[unit_from])
    xarr.convert_to_unit(unit_to,convention=convention)
    xarr.convert_to_unit(unit_from,convention=convention)
    assert all(np.abs((xarr - xvals)/xvals) < threshold)
    assert xarr.units == unit_from

def test_from_quantity():
    arr = np.linspace(-50,50)*u.km/u.s
    xarr = units.SpectroscopicAxis.from_quantity(arr)

    assert np.all(xarr == arr.value)

if __name__=="__main__":
    unit_from='GHz'
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_units='GHz',xtype=units.unit_type_dict[unit_from])
    xarr2=xarr.as_unit('km/s',quiet=False,debug=True)
    xarr3=xarr2.as_unit('GHz',quiet=False,debug=True)
        
