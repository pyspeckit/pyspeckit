import pyspeckit
from .. import units
from ..units import SpectroscopicAxis
import numpy as np
import pytest
from astropy import units as u

convention = ['optical','radio','relativistic']
convention_map = {'optical': u.doppler_optical,
                  'radio': u.doppler_radio,
                  'relativistic': u.doppler_relativistic}

# these are used to populate a matrix for unit conversion tests
# being complete is expensive, ~1+ minutes
unit_type_dict = {
    #'Hz' :'frequency', 'kHz':'frequency', 'MHz':'frequency',
    'GHz':'frequency',
    # not permitted by astropy 'hz' :'frequency', 'khz':'frequency', 'mhz':'frequency', 'ghz':'frequency',
    #'THz':'frequency',
    'meter/second':'velocity',
    #'m/s':'velocity',
    #'kilometer/s':'velocity',
    #'kilometer/second':'velocity',
    #'centimeter/second':'velocity',
    #'megameter/s':'velocity',
    #'megameter/second':'velocity','Mm/s':'velocity',
    'km/s':'velocity',
    #'kms':'velocity',
    #'centimeter/s':'velocity',
    #'kms-1':'velocity',
    #'km s-1':'velocity',
    #'m s-1':'velocity',
    'cm s-1':'velocity',
    #'ms-1':'velocity',
    #'cm/s':'velocity', #'cms':'velocity',
    #'meter':'wavelength','m':'wavelength',
    'centimeter':'wavelength',
    #'cm':'wavelength',
    #'millimeter':'wavelength','mm':'wavelength',
    #'nanometer':'wavelength',
    'nm':'wavelength',
    'micrometer':'wavelength',
    'micron':'wavelength',
    'um':'wavelength',
    #'kilometer':'wavelength','km':'wavelength',
    #'megameter':'wavelength','Mm':'wavelength',
    'angstrom':'wavelength',#'angstroms':'wavelength',
    'AA':'wavelength',
    'unknown':'pixels',
    None: 'pixels',
}

params = [(a,b,c,d) for a in unit_type_dict if a not in ('unknown',None)
          for b in unit_type_dict if b not in ('unknown',None)
          for c in convention
          for d in ['GHz','cm']]

class TestUnits(object):

    @classmethod
    def setup_class(cls):
        cls.x = units.SpectroscopicAxis(np.arange(5), unit=u.angstrom, velocity_convention='optical', equivalencies=u.doppler_optical(3*u.AA))
        cls.length = units.SpectroscopicAxis(np.arange(5), unit=u.nm)
        cls.frequency = units.SpectroscopicAxis(np.arange(5), unit=u.GHz)
        cls.speed = units.SpectroscopicAxis(np.arange(5), unit=u.km/u.s)
        return cls

    @pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
    def test_convert_to(self, unit_from, unit_to, convention, ref_unit):
        try:
            self.x.convert_to_unit(unit_to,
                                   velocity_convention=convention,
                                   make_dxarr=False,
                                  )
            assert self.x.unit == u.Unit(unit_to)
        except TypeError as ex:
            if str(ex.args[0]) == ('cannot convert value type to array type '
                                   'without precision loss'):
                print("{0}->{1} with ref {2} is too imprecise".format(
                    unit_from, unit_to, ref_unit))
            else:
                raise ex

    def test_equivalencies_1(self, ):
        assert self.x.equivalencies is not None # == u.doppler_optical(3*u.AA)

    @pytest.mark.parametrize(('velocity_convention'), convention)
    def test_equivalencies_2(self, velocity_convention):
        ref = 3*u.AA
        x = SpectroscopicAxis(np.arange(5), unit=u.angstrom, refX=ref,
                              velocity_convention=velocity_convention)
        expected = set(convention_map[velocity_convention](ref) + u.spectral())
        eq_units = set([item[:2] for item in x.equivalencies])
        eq_units_expected = set([item[:2] for item in expected])
        assert eq_units == eq_units_expected

    @pytest.mark.parametrize(('velocity_convention'), convention)
    def test_equivalencies_3(self, velocity_convention):
        # note that the ref has no unit; the unit is assumed to be the axis units
        x = SpectroscopicAxis(np.arange(5), unit=u.angstrom, refX=3,
                              refX_unit='angstrom',
                              velocity_convention=velocity_convention)
        expected = set(convention_map[velocity_convention](3*u.AA) + u.spectral())
        eq_units = set([item[:2] for item in x.equivalencies])
        eq_units_expected = set([item[:2] for item in expected])
        assert eq_units == eq_units_expected

    def test_initialize_units(self, ):
        xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=u.dimensionless_unscaled)
        assert xarr.unit == u.dimensionless_unscaled
        xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=u.m)
        assert xarr.unit == u.m

    @pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
    def test_convert_units(self, unit_from, unit_to, convention, ref_unit):
        xarr = units.SpectroscopicAxis(np.linspace(4.9, 5.1, 3), unit=unit_from,
                                       refX=5, refX_unit=ref_unit,
                                       velocity_convention=convention,
                                      )
        xarr.convert_to_unit(unit_to, make_dxarr=False,)
        assert xarr.unit == unit_to

    def test_convert_units2(self, ):
        xarr = units.SpectroscopicAxis(np.linspace(1,10,3),unit='Hz')

        velocity_arr = xarr.as_unit('km/s', velocity_convention='optical', refX=5*u.GHz)

        xarr2 = units.SpectroscopicAxis(np.linspace(1,10,3),unit='Hz',
                                        velocity_convention='optical',
                                        refX=5*u.GHz)
        velocity_arr2 = xarr2.to(u.km/u.s)

        assert np.all(velocity_arr.value==velocity_arr2.value)

    def test_in_range(self, ):
        xarr = units.SpectroscopicAxis(np.linspace(1,10,3),unit='Hz')

        assert not xarr.in_range(10*u.GHz)
        assert xarr.in_range(5*u.Hz)

        # TODO: make sure warning is raised
        assert xarr.in_range(5)


    @pytest.mark.parametrize(('unit_from','unit_to','convention','ref_unit'),params)
    def test_convert_back(self, unit_from, unit_to,convention,ref_unit):
        if unit_from in ('cms','cm/s','centimeter/second') or unit_to in ('cms','cm/s','centimeter/second'):
            xvals = np.linspace(1000,10000,3)
        else:
            xvals = np.linspace(1,10,3)
        assert not np.any(xvals==0.0)
        # all conversions include a * or / by speedoflight_ms
        threshold = np.spacing(units.speedoflight_ms) * 100
        if 'megameter' in unit_from or 'Mm' in unit_from:
            threshold *= 10
        if 'centimeter' in unit_from or 'cm' in unit_from:
            threshold *= 10

        unit_to = u.Unit(unit_to)
        unit_from = u.Unit(unit_from)
        # come up with a sane reference value
        if unit_from.physical_type in ('frequency','length'):
            refX = u.Quantity(5, unit_from).to(ref_unit, u.spectral())
        elif unit_to.physical_type in ('frequency','length'):
            refX = u.Quantity(5, unit_to).to(ref_unit, u.spectral())
        else:
            refX = u.Quantity(5, ref_unit)
        xarr = units.SpectroscopicAxis(np.copy(xvals), unit=unit_from, refX=refX,
                                       velocity_convention=convention)
        xarr.convert_to_unit(unit_to, velocity_convention=convention,
                             make_dxarr=False)
        xarr.convert_to_unit(unit_from, velocity_convention=convention,
                             make_dxarr=False)
        infinites = np.isinf(xarr.value)
        assert all(np.abs((xarr.value - xvals)/xvals)[~infinites] < threshold)
        assert xarr.unit == unit_from

    def test_x_to_pix(self):
        # regression test for issue 150

        xarr = pyspeckit.units.SpectroscopicAxis(np.linspace(-50,50)*u.km/u.s)
        assert xarr.x_to_pix(12*u.km/u.s) == 30
        assert xarr.x_to_pix(12000*u.m/u.s) == 30
        assert xarr.x_to_pix(12000,u.m/u.s) == 30
        assert xarr.x_to_pix(120,u.Mm/u.s) == 49

        with pytest.raises(u.UnitConversionError) as ex:
            xarr.x_to_pix(120,u.GHz)
        assert str(ex.value) == "'GHz' (frequency) and 'km / s' (speed/velocity) are not convertible"

        xarr = pyspeckit.units.SpectroscopicAxis(np.linspace(-50,50)*u.km/u.s, refX=100*u.GHz,
                                                 velocity_convention='radio')
        assert xarr.x_to_pix(100.001,u.GHz) == 23

    def test_comparison(self):
        # regression test for issue 153

        xx = np.linspace(-50,50)*u.km/u.s

        # test the quantity-only version first
        # (so if it crashes here, we know it's not our fault)
        print()
        print("Astropy quantity tests\n")
        mask1 = xx < 0*u.km/u.s
        print("Resulting type: {0}".format(type(mask1)))
        assert np.all(mask1[:25])
        assert np.all(~mask1[25:])

        print()
        print("Pyspeckit xarr tests\n")
        xarr = pyspeckit.units.SpectroscopicAxis(xx)
        mask = xarr < 0*u.km/u.s
        print("Resulting type: {0}".format(type(mask)))
        assert np.all(mask[:25])
        assert np.all(~mask[25:])

    def test_cdelt(self):
        # regression for issue 175

        xx = np.linspace(-50,50)*u.km/u.s
        dx = np.mean(np.diff(xx))
        xarr = pyspeckit.units.SpectroscopicAxis(xx, refX=5*u.GHz)

        np.testing.assert_almost_equal(xarr.cdelt().value, dx.value)
        assert len(xarr.dxarr) == len(xarr)

        xarr = pyspeckit.units.SpectroscopicAxis(xx, refX=5*u.GHz)

        # check both orders: dxarr is a property, so this could be different
        assert len(xarr.dxarr) == len(xarr)
        np.testing.assert_almost_equal(xarr.cdelt().value, dx.value)

        xarr2 = xarr[10:-10]
        np.testing.assert_almost_equal(xarr2.cdelt().value, dx.value)
        assert len(xarr2) == len(xarr2.dxarr)




if __name__=="__main__":
    unit_from='GHz'
    xarr = units.SpectroscopicAxis(np.linspace(1,10,10),unit=unit_from,refX=5,refX_unit='GHz')
    xarr2=xarr.as_unit('km/s',quiet=False,debug=True)
    xarr3=xarr2.as_unit('GHz',quiet=False,debug=True)
