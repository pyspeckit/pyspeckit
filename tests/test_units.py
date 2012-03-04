import pyspeckit
from pyspeckit.spectrum import units
import numpy as np
import pytest

params = [(a,b) for a in units.unit_type_dict if a not in ('unknown',None) for b in units.unit_type_dict if b not in ('unknown',None)]

@pytest.mark.parametrize(('unit_from','unit_to'),params)
def test_convert_units(unit_from,unit_to):
    xarr = units.SpectroscopicAxis(np.linspace(0,10,10),unit=unit_from,refX=5,refX_units='GHz',xtype=units.unit_type_dict[unit_from])
    xarr.convert_to_unit(unit_to)
    assert xarr.units == unit_to

