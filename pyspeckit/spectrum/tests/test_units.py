import pyspeckit
from pyspeckit.spectrum import units
import numpy as np
import pytest

convention = ['optical','radio','relativistic']
params = [(a,b,c,d) for a in units.unit_type_dict if a not in ('unknown',None)
                  for b in units.unit_type_dict if b not in ('unknown',None)
                  for c in convention
                  for d in ['GHz','cm']]

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


if __name__=="__main__":
    unit_from='GHz'
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_units='GHz',xtype=units.unit_type_dict[unit_from])
    xarr2=xarr.as_unit('km/s',quiet=False,debug=True)
    xarr3=xarr2.as_unit('GHz',quiet=False,debug=True)
        
