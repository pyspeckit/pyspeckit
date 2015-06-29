"""
===========================
Units and SpectroscopicAxes
===========================

Unit parsing and conversion tool.
The SpectroscopicAxis class is meant to deal with unit conversion internally

Open Questions: Are there other FITS-valid projection types, unit types, etc.
that should be included?
What about for other fields (e.g., wavenumber?)

"""
from __future__ import print_function
import numpy as np
import warnings
from astropy import units as u
from astropy import log

# declare a case-insensitive dict class to return case-insensitive versions of each dictionary...
# this is just a shortcut so that units can be specified as, e.g., Hz, hz, HZ, hZ, etc.  3 of those 4 are "legitimate".  
# the templates for this code were grabbed from a few StackOverflow.com threads
# I changed this to a "SmartCase" dict, by analogy with VIM's smartcase, where anything with caps will
# be case sensitive, BUT if there is no unit matching the exact case, will search for the lower version too
# e.g., if you search for Mm, it will return Megameters.  If you search for mm or mM, it will get millimeters
class SmartCaseNoSpaceDict(dict):
    def __init__(self, inputdict=None):
        if inputdict:
            # Doesn't do keyword args
            if isinstance(inputdict, dict):
                self.CaseDict = inputdict
                for k,v in inputdict.items():
                    self[k] = v
                    #dict.__setitem__(self, k, v)
                    #if hasattr(k,'lower'):
                    #    dict.__setitem__(self, k.lower(), v)
            else:
                for k,v in inputdict:
                    self[k] = v
                    #dict.__setitem__(self, k, v)
                    #if hasattr(k,'lower'):
                    #    dict.__setitem__(self, k.lower(), v)

    def __getitem__(self, key):
        try: 
            return dict.__getitem__(self,key)
        except KeyError:
            return dict.__getitem__(self, key.lower().replace(" ",""))

    def __setitem__(self, key, value):
        dict.__setitem__(self,key,value)
        # only set lower-case version if there isn't one there already
        # prevents overwriting mm with Mm.lower()
        if hasattr(key,'lower'):
            if key.lower().replace(" ","") not in self:
                self.__setitem__(key.lower().replace(" ",""), value)

    def __contains__(self, key):
        cased = dict.__contains__(self,key)
        if hasattr(key,'lower') and cased is False:
            return dict.__contains__(self, key.lower().replace(" ",""))
        else:
            return cased

    def has_key(self, key):
        """ This is deprecated, but we're keeping it around """
        cased = dict.has_key__(self,key)
        if hasattr(key,'lower') and cased is False:
            return dict.has_key__(self, key.lower().replace(" ",""))
        else:
            return cased

    def get(self, key, def_val=None):
        cased = dict.get(self,key)
        if hasattr(key,'lower') and cased is None:
            return dict.get(self, key.lower().replace(" ",""),def_val)
        else:
            return cased

    def setdefault(self, key, def_val=None):
        """
        If key is in the dictionary, return its value. If not, insert key with
        a value of default and return default. default defaults to None.
        """
        cased = dict.setdefault(self,key)
        if hasattr(key,'lower') and cased is None:
            return dict.setdefault(self, key.lower().replace(" ",""),def_val)
        else:
            return cased

    def update(self, inputdict):
        for k,v in inputdict.items():
            # let self.__setitem__ hand the low/upperness
            self.__setitem__(self, k, v)

    def fromkeys(self, iterable, value=None):
        d = SmartCaseDict()
        for k in iterable:
            self.__setitem__(d, k, value)
        return d

    def pop(self, key, def_val=None):
        cased = dict.pop(self,key)
        if hasattr(key,'lower') and cased is None:
            return dict.pop(self, key.lower().replace(" ",""),def_val)
        else:
            return cased
    

length_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'megameters':1e6,'Mm':1e6,
        'angstrom':1e-10,'angstroms':1e-10,'A':1e-10,
        }
length_dict = SmartCaseNoSpaceDict(length_dict)

wavelength_dict = length_dict # synonym

frequency_dict = {
        'Hz':1.0,
        'kHz':1e3,
        'MHz':1e6,
        'GHz':1e9,
        'THz':1e12,
        }
frequency_dict = SmartCaseNoSpaceDict(frequency_dict)

velocity_dict = {'meters/second':1.0,'m/s':1.0,'m s-1':1.0,'ms-1':1.0,
                 'kilometers/second':1e3,'kilometers/s':1e3,'km/s':1e3,'kms':1e3,'km s-1':1e3,'kms-1':1e3,
                'centimeters/second':1e-2,'centimeters/s':1e-2,'cm/s':1e-2,'cms':1e-2,
                'megameters/second':1e6,'megameters/s':1e6,'Mm/s':1e6,'Mms':1e6,
        }
velocity_dict = SmartCaseNoSpaceDict(velocity_dict)
velocity_conventions = {'optical': u.doppler_optical, 'radio': u.doppler_radio, 'relativistic': u.doppler_relativistic}

pixel_dict = {'pixel':1,'pixels':1}
pixel_dict = SmartCaseNoSpaceDict(pixel_dict)

conversion_dict = {
        'VELOCITY':velocity_dict,  'Velocity':velocity_dict,  'velocity':velocity_dict,  'velo': velocity_dict, 'VELO': velocity_dict,
        'LENGTH':length_dict,      'Length':length_dict,      'length':length_dict, 
        'WAVELENGTH':length_dict,      'WAVELength':length_dict,      'WAVElength':length_dict, 
        'FREQUENCY':frequency_dict,'Frequency':frequency_dict,'frequency':frequency_dict, 'freq': frequency_dict, 'FREQ': frequency_dict,
        'pixels':pixel_dict,'PIXELS':pixel_dict,
        }
conversion_dict = SmartCaseNoSpaceDict(conversion_dict)


# to-do
# unit_prettyprint = dict([(a,a.replace('angstroms','\\AA')) for a in unit_type_dict.keys() if a is not None])

xtype_dict = {
        'VLSR':'velocity','VRAD':'velocity','VELO':'velocity',
        'VELO-LSR':'velocity',
        'VOPT':'velocity',
        'VHEL':'velocity',
        'VGEO':'velocity',
        'VREST':'velocity',
        'velocity':'velocity',
        'Z':'redshift',
        'FREQ':'frequency',
        'frequency':'frequency',
        'WAV':'wavelength',
        'WAVE':'wavelength',
        'wavelength':'wavelength',
        'lambda':'wavelength',
        # cm^-1 ? 'wavenumber':'wavenumber',
        }

frame_dict = SmartCaseNoSpaceDict({
        'VLSR':'LSRK','VRAD':'LSRK','VELO':'LSRK',
        'VOPT':'LSRK',
        'VELO-LSR':'LSRK',
        'LSRD':'LSRD',
        '-LSR':'LSRK',
        '-HEL':'heliocentric',
        '-BAR':'barycentric',
        'BAR':'barycentric',
        'BARY':'barycentric',
        '-OBS':'obs',
        'VHEL':'heliocentric',
        'VGEO':'geocentric',
        'topo':'geocentric',
        'topocentric':'geocentric',
        'velocity':'LSRK',
        'VREST':'rest',
        'Z':'rest',
        'FREQ':'rest',
        'frequency':'rest',
        'WAV':'rest',
        'WAVE':'rest',
        'wavelength':'rest',
        'lambda':'rest',
        'redshift':'obs'
        })

frame_type_dict = {'LSRK':'velocity','LSRD':'velocity','LSR':'velocity',
        'heliocentric':'velocity','topocentric':'velocity','geocentric':'velocity',
        'rest':'frequency','obs':'frequency','observed':'frequency'}

fits_frame = {'rest':'REST','LSRK':'-LSR','heliocentric':'-HEL','geocentric':'-GEO'}
fits_specsys = {'rest':'REST','LSRK':'LSRK','LSRD':'LSRD','heliocentric':'HEL','geocentric':'GEO'}
fits_type = {'velocity':'VELO','frequency':'FREQ','wavelength':'WAVE','length':'WAVE','redshift':'REDS',
        'Velocity':'VELO','Frequency':'FREQ','Wavelength':'WAVE','Length':'WAVE','Redshift':'REDS',
        u.Hz.physical_type:'FREQ', (u.km/u.s).physical_type:'VELO', u.m.physical_type:'WAVE'}
convention_suffix = {'radio':'RAD','optical':'OPT','relativistic':'REL','redshift':'RED'}

speedoflight_ms  = np.float64(2.99792458e8) # m/s
speedoflight_kms = np.float64(2.99792458e5) # km/s
c = speedoflight_cms = np.float64(2.99792458e10) # cm/s

astronomical_unit_cm = np.float64(1.49597870e13)

hplanck = h = 6.626068e-27    # Planck's Constant (erg s)
k = kb = kB = 1.3806503e-16   # Boltzmann's Constant (erg / K)
electronmass = me = 9.10938188e-28 # g
electroncharge = 4.80320427e-10 # esu

unitdict = {
        'cgs':{ 'h':6.626068e-27,
            'k':1.3806503e-16,
            'kb':1.3806503e-16,
            'c':2.99792458e10,
            'mh':1.67262158e-24 * 1.00794, # proton mass
            'me':9.10938188e-28,           # electron mass
            'e':4.80320427e-10, # electron charge, esu
            'speed':'cm/s',
            'length':'cm'},
        'mks':{ 'h':6.626068e-34,
            'k':1.3806503e-23,
            'kb':1.3806503e-23,
            'c':2.99792458e8,
            'mh':1.67262158e-27 * 1.00794,
            'me':9.10938188e-31,
            'speed':'m/s',
            'length':'m'}
        }


def generate_xarr(input_array, unit=None):
    """
    Create an 'xarr' from an array; can be a quantity

    unit is ignored unless the input is a simple ndarray
    """
    from astropy import units as u
    if hasattr(input_array, 'value') and hasattr(input_array, 'unit'):
        return SpectroscopicAxis(input_array.value,
                                 unit=input_array.unit.to_string())
    elif isinstance(input_array, SpectroscopicAxis):
        return input_array
    elif isinstance(input_array, np.ndarray):
        return SpectroscopicAxis(input_array, unit=unit)
    else:
        raise TypeError("Unrecognized input type. Input array of type: {0}"+\
            "is not a Quantity, SpectroscopicAxis or numpy.ndarray".format(type(input_array)))

class SpectroscopicAxis(u.Quantity):
    """
    A Spectroscopic Axis object to store the current units of the spectrum and
    allow conversion to other units and frames.  Typically, units are velocity,
    wavelength, frequency, or redshift.  Wavenumber is also hypothetically
    possible.

    WARNING: If you index a SpectroscopicAxis, the resulting array will be
    a SpectroscopicAxis without a dxarr attribute!  This can result in major problems;
    a workaround is being sought but subclassing numpy arrays is harder than I thought
    """
    def __new__(self, xarr, unit=None, refX=None, redshift=None,
                refX_unit=None, velocity_convention=None, use128bits=False, 
                bad_unit_response='raise', equivalencies=u.spectral(),
                center_frequency=None, center_frequency_unit=None):
        """
        Make a new spectroscopic axis instance
        Default units Hz
        
        Parameters
        ----------
        xarr : np.ndarray
            An array of X-axis values in whatever unit specified
        unit : str
            Any valid spectroscopic X-axis unit (km/s, Hz, angstroms, etc.).
            Spaces will be removed.
        refX : float
            Reference frequency/wavelength
        refX_unit : str | astropy.units.Unit
            Units of the reference frequency/wavelength
        center_frequency: float
            The reference frequency for determining a velocity. 
            Required for conversions between frequency/wavelength/energy and velocity.
        center_frequency_unit: string
            If converting between velocity and any other spectroscopic type,
            need to specify the central frequency around which that velocity is
            calculated.
        equivalencies : list
            astropy equivalencies list containing tuples of the form:
            (from_unit, to_unit, forward, backward)
            forward and backward are functions that convert values between those units
        bad_unit_response : 'raise' | 'pixel'
            What should pyspeckit do if the units are not recognized?  Default
            is to raise an exception.  Can make pixel units instead
        """        
        if use128bits:
            dtype='float128'
        else:
            dtype='float64'

        # Only need to convert xarr to array if it's not already one (e.g., if
        # it's a list)
        if not isinstance(xarr, np.ndarray):
            subarr = np.array(xarr, dtype=dtype)
        else:
            subarr = xarr
        
        subarr = subarr.view(self)

        # Only need to set xarr's unit if it's not a quantity
        # or if it is unitless
        if not isinstance(xarr, u.Quantity) or xarr.unit == u.dimensionless_unscaled:
            subarr._unit = self.validate_unit(unit, bad_unit_response)

        subarr.refX = refX

        if refX_unit is None:
            if hasattr(refX, 'unit'):
                refX_unit = refX.unit
            elif subarr._unit in frequency_dict:
                refX_unit = subarr.unit
            else:
                refX_unit = 'Hz'
        subarr.refX_unit = refX_unit

        subarr.redshift = redshift
        subarr.wcshead = {}
        subarr.velocity_convention = velocity_convention

        if not center_frequency:
            if refX:
                center_frequency = refX    
        if not center_frequency_unit:
            if refX_unit:
                center_frequency_unit = refX_unit
            else:
                center_frequency_unit = unit
        if center_frequency:
            subarr.center_frequency = u.Quantity(center_frequency, center_frequency_unit)
        else:
            subarr.center_frequency = None

        subarr.center_frequency, subarr._equivalencies = \
                self.find_equivalencies(subarr.velocity_convention,
                                        subarr.center_frequency, 
                                        equivalencies)
        return subarr

    def __getitem__(self, key):
        """
        We do *NOT* want to return a SpectroscopicAxis when indexed singly!
        """
        if self.isscalar:
            raise TypeError(
                "'{cls}' object with a scalar value does not support "
                "indexing".format(cls=self.__class__.__name__))

        out = super(u.Quantity, self).__getitem__(key)
        if np.isscalar(out):
            return u.Quantity(out, unit=self.unit)
        else:
            return self._new_view(out)

    def __repr__(self):
        if self.shape is ():
            rep = ("SpectroscopicAxis([%r], unit=%r, refX=%r, refX_unit=%r, frame=%r, redshift=%r, xtype=%r, velocity convention=%r)" %
                (self.value, self.unit, self.refX, self.refX_unit, self.frame, self.redshift, self.xtype, self.velocity_convention))
        else:
            rep = ("SpectroscopicAxis([%r,...,%r], unit=%r, refX=%r, refX_unit=%r, frame=%r, redshift=%r, xtype=%r, velocity convention=%r)" %
                (self[0].value, self[-1].value, self.unit, self.refX, self.refX_unit, self.frame, self.redshift, self.xtype, self.velocity_convention))
        return rep

    def __str__(self):
        selfstr =  "SpectroscopicAxis with units %s and range %g:%g." % (
                self.unit,self.umin().value,self.umax().value)
        if self.refX is not None:
            if not hasattr(self.refX, 'unit'):
                selfstr += "Reference is %g %s" % (self.refX, self.refX_unit)
            else:
                selfstr += "Reference is %s" % (self.refX)
        return selfstr

    @property
    def units(self):
        log.warn("'units' is deprecated; please use 'unit'", DeprecationWarning)
        return self._unit

    @units.setter
    def units(self, value):
        log.warn("'units' is deprecated; please use 'unit'", DeprecationWarning)
        self.set_unit(value)

    @property
    def refX_units(self):
        log.warn("'refX_units' is deprecated; please use 'refX_unit'", DeprecationWarning)
        return self.refX_unit

    @refX_units.setter
    def refX_units(self, value):
        log.warn("'refX_units' is deprecated; please use 'refX_unit'", DeprecationWarning)
        self.refX_unit = value

    def __array_finalize__(self,obj):
        """
        Array finalize allows you to take subsets of the array but keep units
        around, e.g.:
        xarr = self[1:20]
        """
        self._unit = getattr(obj, 'unit', u.dimensionless_unscaled)
        self.frame = getattr(obj, 'frame', None)
        self.xtype = getattr(obj, 'xtype', None)
        self.refX = getattr(obj, 'refX', None)
        self.refX_unit = getattr(obj, 'refX_unit', None)
        self.velocity_convention = getattr(obj, 'velocity_convention', None)
        self.redshift = getattr(obj, 'redshift', None)
        self.wcshead = getattr(obj, 'wcshead', None)
        self.center_frequency = getattr(obj, 'center_frequency', None)
        self.center_frequency_unit = getattr(obj, 'center_frequency_unit', None)
        self._equivalencies = getattr(obj, 'equivalencies', [])
        # moved from __init__ - needs to be done whenever viewed
        # (this is slow, though - may be better not to do this)
        if self.shape: # check to make sure non-scalar
            # can't use make_dxarr here! infinite recursion =(
            self.dxarr = np.diff(np.array(self))

    def __array_wrap__(self,out_arr,context=None):
        """
        Do this when calling ufuncs
        """
        ret = np.ndarray.__array_wrap__(self, out_arr, context)
        try:
            if len(ret) == 1: # I think this is a bad hack; it used to be unnecessary but a recent update to numpy changed that?
                ret = float(ret) # FORCE scalar
        except TypeError: # sometimes ret has no len
            pass
        return ret

    def _check_consistent_type(self):
        """
        Make sure self.xtype is unit_type_dict[unit]
        if this is NOT true, can cause significant errors
        """
        if self.xtype is None:
            OK = False
        else:
            OK = self.xtype.lower() == self.unit

        if not OK:
            raise InconsistentTypeError("Unit: %s Type[unit]: %s in %r" % (self.unit,unit_type_dict[self.unit],self))

    def _fix_type(self):
        try:
            self._check_consistent_type()
        except InconsistentTypeError:
            self.xtype = unit_type_dict[self.unit]

    @classmethod
    def validate_unit(self, unit, bad_unit_response='raise'):
        try:
            if unit is None or unit == 'unknown':
                unit = u.dimensionless_unscaled
            if unit == 'angstroms': unit = 'angstrom'
            unit = u.Unit(unit)
        except ValueError:
            if bad_unit_response == "pixel":
                unit = u.dimensionless_unscaled
            elif bad_unit_response == "raise":
                raise ValueError('Unit %s not recognized.' % unit)
            else:
                raise ValueError('Unit %s not recognized. Invalid bad_unit_response, valid options: [%s/%s]' 
                                % (unit, "raise", "pixel"))
        return unit

    def set_unit(self, unit, bad_unit_response='raise'):
        self._unit = self.validate_unit(unit, bad_unit_response)

    def umax(self, unit=None):
        """
        Return the maximum value of the SpectroscopicAxis.  If units specified,
        convert to those units first
        """

        if unit is not None:
            # could be reversed
            return np.max([self.coord_to_x(self.max(), units),self.coord_to_x(self.min(),units)])
        else: 
            return self.max()

    def umin(self, unit=None):
        """
        Return the minimum value of the SpectroscopicAxis.  If units specified,
        convert to those units first
        """

        if unit is not None:
            # could be reversed
            return np.min([self.coord_to_x(self.max(), units),self.coord_to_x(self.min(),units)])
        else: 
            return self.min()

    def x_to_pix(self, xval, xval_units=None):
        """
        Given an X coordinate in SpectroscopicAxis' units, return the corresponding pixel number
        """
        if xval_units in pixel_dict:
            return xval
        else:
            try:
                if not xval_units:
                    xval_units = xval.unit
            except AttributeError:
                xval_units = u.dimensionless_unscaled
            xval = xval * u.Unit(xval_units)            
            nearest_pix = np.argmin(np.abs(self.value-xval.value))
            # nearest_pix = np.argmin(np.abs(self.value-xval))
            return nearest_pix

    def in_range(self, xval):
        """
        Given an X coordinate in SpectroscopicAxis' units, return whether the pixel is in range
        """
        return bool((xval > self.min()) * (xval < self.max()))

    def x_to_coord(self, xval, xunit, verbose=False):
        """
        Given a wavelength/frequency/velocity, return the value in the SpectroscopicAxis's units
        e.g.:
        xarr.unit = 'km/s'
        xarr.refX = 5.0
        xarr.refX_unit = GHz
        xarr.x_to_coord(5.1,'GHz') == 6000 # km/s
        """
        if not hasattr(xval, 'unit'):
             xval = xval * u.Unit(xunit)
        return xval.to(self.unit, self.equivalencies)

    def coord_to_x(self, xval, xunit):
        """
        Given an X-value assumed to be in the coordinate axes, return that
        value converted to xunit
        e.g.:
        xarr.units = 'km/s'
        xarr.refX = 5.0 
        xarr.refX_unit = GHz
        xarr.coord_to_x(6000,'GHz') == 5.1 # GHz
        """
        return self.as_unit(xunit)

    def convert_to_unit(self, unit, **kwargs):
        """
        Return the X-array in the specified units without changing it
        Uses as_unit for the conversion, but changes internal values rather
        than returning them.
        """
        unit = self.validate_unit(unit)

        self.flags.writeable=True

        # Hack:
        # You can only set self[:] to a quantity
        # The quantity will automatically be converted to present units
        # Therefore, get the new array to set the *values* correctly,
        # then input the values and let the unit conversion take care of the
        # numerical changes, THEN set the unit.
        new_values = self.as_unit(unit, **kwargs)
        self[:] = new_values.value * self.unit
        self.set_unit(unit)

        self.flags.writeable=False
        self.make_dxarr()

    def as_unit(self, unit, equivalencies=[], velocity_convention=None, refX=None,
                refX_unit=None, center_frequency=None, center_frequency_unit=None,
                **kwargs):
        """
        Convert the spectrum to the specified units.  This is a wrapper function
        to convert between frequency/velocity/wavelength and simply change the 
        units of the X axis.  Frame conversion is... not necessarily implemented.

        Parameters
        ----------
        unit : string
            What unit do you want to 'view' the array as?
            None returns the x-axis unchanged (NOT a copy!)
        frame : string
            NOT IMPLEMENTED.  When it is, it will allow you to convert between
            LSR, topocentric, heliocentric, rest, redshifted, and whatever other
            frames we can come up with.  Right now the main holdup is finding a 
            nice python interface to an LSR velocity calculator... and motivation.
        center_frequency: float
            The reference frequency for determining a velocity. 
            Required for conversions between frequency/wavelength/energy and velocity.
        center_frequency_unit: string
            If converting between velocity and any other spectroscopic type,
            need to specify the central frequency around which that velocity is
            calculated.
            I think this can also accept wavelengths....
        """
        if not velocity_convention:
            velocity_convention = self.velocity_convention
        if not equivalencies:
            equivalencies = kwargs.get('equivalencies', self.equivalencies)
        if not refX:
            refX = self.refX
        if not refX_unit:
            refX_unit = self.refX_unit
        if not center_frequency:
            if self.center_frequency:
                center_frequency = self.center_frequency
            elif refX:
                center_frequency = refX
        if not center_frequency_unit:
            if refX_unit:
                center_frequency_unit = refX_unit
        try:
            if not hasattr(center_frequency, 'unit'):
                center_frequency = u.Quantity(center_frequency, center_frequency_unit)
        except:
            center_frequency = None

        self.center_frequency, self._equivalencies = \
            self.find_equivalencies(velocity_convention, 
                                    center_frequency,
                                    equivalencies)
        
        if isinstance(self.unit, str):
            self._unit = u.Unit(self.unit)

        return self.to(unit, equivalencies=self.equivalencies)

    def make_dxarr(self, coordinate_location='center'):
        """
        Create a "delta-x" array corresponding to the X array.

        Parameters
        ----------
        coordinate_location : [ 'left', 'center', 'right' ]
            Does the coordinate mark the left, center, or right edge of the
            pixel?  If 'center' or 'left', the *last* pixel will have the same
            dx as the second to last pixel.  If right, the *first* pixel will
            have the same dx as the second pixel.

        """
        dxarr = np.diff(self)
        if self.size <= 2:
            self.dxarr = np.ones(self.size)*dxarr
        elif coordinate_location in ['left','center']:
            self.dxarr = np.concatenate([dxarr,dxarr[-1:]])
        elif coordinate_location in ['right']:
            self.dxarr = np.concatenate([dxarr[:1],dxarr])
        else:
            raise ValueError("Invalid coordinate location.")

    def cdelt(self, tolerance=1e-8, approx=False):
        """
        Return the channel spacing if channels are linear

        Parameters
        ----------
        tolerance : float
            Tolerance in the difference between pixels that determines
            how near to linear the xarr must be
        approx : bool
            Return the mean DX even if it is inaccurate
        """
        if not hasattr(self,'dxarr'): # if cropping happens...
            self.make_dxarr()
        if len(self) <= 1:
            raise ValueError("Cannot have cdelt of length-%i array" % len(self))
        if approx or abs(self.dxarr.max()-self.dxarr.min())/abs(self.dxarr.min()) < tolerance:
            return self.dxarr.mean().flat[0]

    def _make_header(self, tolerance=1e-8):
        """
        Generate a set of WCS parameters for the X-array
        """
        self.make_dxarr()

        if self.wcshead is None:
            self.wcshead = {}

        self.wcshead['CUNIT1'] = self.unit
        if fits_type[self.unit.physical_type] == 'VELO' and self.velocity_convention is not None:
            ctype = 'V'+convention_suffix[self.velocity_convention]
        else:
            ctype = fits_type[self.unit.physical_type]
        #self.wcshead['CTYPE1'] = ctype+fits_frame[self.frame]
        try:
            from spectral_cube import spectral_axis
            self.wcshead['CTYPE1'] = spectral_axis.determine_ctype_from_vconv(ctype,
                                                                              self.unit,
                                                                              self.velocity_convention)
        except ImportError:
            self.wcshead['CTYPE1'] = ctype
        # not needed ? self.wcshead['SPECSYS'] = fits_specsys[self.frame]
        self.wcshead['RESTFRQ'] = self.refX

        # check to make sure the X-axis is linear
        cdelt = self.cdelt(tolerance=tolerance)
        if cdelt is None:
            self.wcshead['CDELT1'] = None
            self.wcshead['CRPIX1'] = None
            self.wcshead['CRVAL1'] = None
            return False
        else:
            self.wcshead['CDELT1'] = cdelt
            self.wcshead['CRVAL1'] = self[0]
            self.wcshead['CRPIX1'] = 1.0
            return True

    @classmethod
    def find_equivalencies(self, velocity_convention=None,
                           center_frequency=None, equivalencies=[]):
        """
        Utility function to add equivalencies from the velocity_convention 
        and the center_frequency

        Parameters
        ----------
        velocity_convention : str
            'optical', 'radio' or 'relativistic'
        center_frequency : float | astropy.units.Quantity 
            The reference frequency for determining a velocity. 
            Required for conversions between frequency/wavelength/energy and velocity.
        equivalencies : list
            astropy equivalencies list containing tuples of the form:
            (from_unit, to_unit, forward, backward)
            forward and backward are functions that convert values between those units
        """
        if velocity_convention and center_frequency:
            new_equivalencies = velocity_conventions[velocity_convention](center_frequency)
            return center_frequency, merge_equivalencies(new_equivalencies, equivalencies)
        else:
            return center_frequency, merge_equivalencies(equivalencies, u.spectral())

def merge_equivalencies(old_equivalencies, new_equivalencies):
    """
    Utility method to merge two equivalency lists
    Uses a dict with concatenated units as keys
    """
    seen = {}
    result = []
    total_equivalencies = old_equivalencies+new_equivalencies

    for equivalency in total_equivalencies:
        equivalency_id = equivalency[0].to_string()+equivalency[1].to_string()
        if equivalency_id in seen: continue
        seen[equivalency_id] = 1
        result.append(equivalency)
    return result


class SpectroscopicAxes(SpectroscopicAxis):
    """
    Counterpart to Spectra: takes a list of SpectroscopicAxis's and
    concatenates them while checking for consistency and maintaining
    header parameters
    """

    def __new__(self, axislist, frame='rest', xtype=None, refX=None,
            redshift=None):

        if type(axislist) is not list:
            raise TypeError("SpectroscopicAxes must be initiated with a list of SpectroscopicAxis objects")

        unit = axislist[0].unit
        xtype = axislist[0].xtype
        redshift = axislist[0].redshift
        velocity_convention = axislist[0].velocity_convention
        for ax in axislist:
            if ax.xtype != xtype:
                try: 
                    ax.change_xtype(xtype)
                except:
                    ValueError("Axis had wrong xtype and could not be converted.")
            if ax.unit != unit:
                try:
                    ax.convert_to_unit(unit)
                except:
                    ValueError("Axis had wrong units and could not be converted.")

        subarr = np.concatenate([ax for ax in axislist])
        subarr = subarr.view(self)
        subarr._unit = unit
        subarr.xtype = xtype
        subarr.redshift = redshift
        subarr.velocity_convention = velocity_convention

        # if all the spectra have the same reference frequency, there is one common refX
        # else, refX is undefined and velocity transformations should not be done
        refXs_diff = np.sum([axislist[0].refX != ax.refX for ax in axislist])
        if refXs_diff > 0:
            subarr.refX = None
            subarr.refX_unit = None
        else:
            subarr.refX = axislist[0].refX
            subarr.refX_unit = axislist[0].refX_unit

        return subarr

class EchelleAxes(SpectroscopicAxis):
    """
    Counterpart to Spectra: takes a list of SpectroscopicAxis's and
    stacks them while checking for consistency and maintaining
    header parameters
    """

    def __new__(self, axislist, frame='rest', xtype=None, refX=None,
            redshift=None):

        if type(axislist) is not list:
            raise TypeError("SpectroscopicAxes must be initiated with a list of SpectroscopicAxis objects")

        unit = axislist[0].unit
        xtype = axislist[0].xtype
        redshift = axislist[0].redshift
        velocity_convention = axislist[0].velocity_convention
        for ax in axislist:
            if ax.xtype != xtype:
                try: 
                    ax.change_xtype(xtype)
                except:
                    ValueError("Axis had wrong xtype and could not be converted.")
            if ax.unit != unit:
                try:
                    ax.convert_to_unit(unit)
                except:
                    ValueError("Axis had wrong units and could not be converted.")

        subarr = np.array(axislist)
        subarr = subarr.view(self)
        subarr._unit = unit
        subarr.xtype = xtype
        subarr.redshift = redshift
        subarr.velocity_convention = velocity_convention

        # if all the spectra have the same reference frequency, there is one common refX
        # else, refX is undefined and velocity transformations should not be done
        refXs_diff = np.sum([axislist[0].refX != ax.refX for ax in axislist])
        if refXs_diff > 0:
            subarr.refX = None
            subarr.refX_unit = None
        else:
            subarr.refX = axislist[0].refX
            subarr.refX_unit = axislist[0].refX_unit

        return subarr

def velocity_to_frequency(velocities, input_units, center_frequency=None,
        center_frequency_units=None, frequency_units='Hz',
        convention='radio'):
    """

    Parameters
    ----------
    velocities : np.ndarray
        An array of velocity values
    inpput_units : str
        A string representing the units of the velocities array
    center_frequency : float
        The reference frequency (i.e., the 0-m/s freq)
    center_frequency_units : str
        A string representing the units of the reference frequency
    frequency_units : str
        A string representing the desired output units
    convention : ['radio','optical','relativistic':
        Conventions defined here:
        http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
        * Radio 	V = c (f0 - f)/f0 	f(V) = f0 ( 1 - V/c )
        * Optical 	V = c (f0 - f)/f 	f(V) = f0 ( 1 + V/c )^-1
        * Redshift 	z = (f0 - f)/f 	f(V) = f0 ( 1 + z )-1
        * Relativistic 	V = c (f02 - f 2)/(f02 + f 2) 	f(V) = f0 { 1 - (V/c)2}1/2/(1+V/c) 

    """
    if input_units in frequency_dict:
        #print("Already in frequency units (%s)" % input_units)
        return velocities
    if center_frequency is None:
        raise ValueError("Cannot convert velocity to frequency without specifying a central frequency.")
    if frequency_units not in frequency_dict:
        raise ValueError("Bad frequency units: %s" % (frequency_units))

    velocity_ms = velocities / velocity_dict['m/s'] * velocity_dict[input_units]
    if convention == 'radio':
        freq = (velocity_ms / speedoflight_ms - 1.0) * center_frequency * -1
    elif convention == 'optical':
        freq = (velocity_ms / speedoflight_ms + 1.0)**-1 * center_frequency
    elif convention == 'relativistic':
        freq = (-(velocity_ms / speedoflight_ms)**2 + 1.0)**0.5 / (1.0 + velocity_ms/speedoflight_ms) * center_frequency
    else:
        raise ValueError('Convention "%s" is not allowed.' % (convention))
    frequencies = freq / frequency_dict[frequency_units] * frequency_dict[center_frequency_units]
    return frequencies

def frequency_to_velocity(frequencies, input_units, center_frequency=None,
                          center_frequency_units=None, velocity_units='m/s',
                          convention='radio'):
    """
    Conventions defined here:
    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
     
     * Radio 	V = c (f0 - f)/f0 	f(V) = f0 ( 1 - V/c )
     * Optical 	V = c (f0 - f)/f 	f(V) = f0 ( 1 + V/c )-1
     * Redshift 	z = (f0 - f)/f 	f(V) = f0 ( 1 + z )-1
     * Relativistic 	V = c (f02 - f 2)/(f02 + f 2) 	f(V) = f0 { 1 - (V/c)2}1/2/(1+V/c)
    """
    if input_units in velocity_dict:
        print("Already in velocity units (%s)" % input_units)
        return frequencies
    if center_frequency is None:
        raise ValueError("Cannot convert frequency to velocity without specifying a central frequency.")
    if center_frequency_units not in frequency_dict:
        raise ValueError("Bad frequency units: %s" % (center_frequency_units))
    if velocity_units not in velocity_dict:
        raise ValueError("Bad velocity units: %s" % (velocity_units))

    frequency_hz = frequencies / frequency_dict['Hz'] * frequency_dict[input_units]
    center_frequency_hz = center_frequency / frequency_dict['Hz'] * frequency_dict[center_frequency_units]

    # the order is very ugly because otherwise, if scalar, the spectroscopic axis attributes won't be inherited
    if convention == 'radio':
        velocity = ( frequency_hz - center_frequency_hz ) / center_frequency_hz * speedoflight_ms * -1
    elif convention == 'optical':
        velocity = ( frequency_hz - center_frequency_hz ) / frequency_hz * speedoflight_ms * -1
    elif convention == 'relativistic':
        velocity = ( frequency_hz**2 - center_frequency_hz**2 ) / ( center_frequency_hz**2 + frequency_hz )**2 * speedoflight_ms * -1
    else:
        raise ValueError('Convention "%s" is not allowed.' % (convention))
    velocities = velocity * velocity_dict['m/s'] / velocity_dict[velocity_units]

    return velocities

def frequency_to_wavelength(frequencies, input_units, wavelength_units='um'):
    """
    Simple conversion from frequency to wavelength:
    lambda = c / nu
    """
    if input_units in wavelength_dict:
        print("Already in wavelength units (%s)" % input_units)
        return
    if wavelength_units not in length_dict:
        raise ValueError("Wavelength units %s not valid" % wavelength_units)

    if input_units not in frequency_dict:
        raise AttributeError("Cannot convert from frequency unless units are already frequency.")

    wavelengths = speedoflight_ms / ( frequencies * frequency_dict[input_units] ) / length_dict[wavelength_units]

    return wavelengths

def wavelength_to_frequency(wavelengths, input_units, frequency_units='GHz'):
    """
    Simple conversion from frequency to wavelength:
    nu = c / lambda
    """
    if input_units in frequency_dict:
        #print "Already in frequency units (%s)" % input_units
        return wavelengths
    if frequency_units not in frequency_dict:
        raise ValueError("Frequency units %s not valid" % frequency_units)

    if input_units not in length_dict:
        raise AttributeError("Cannot convert from wavelength unless units are already wavelength.")

    frequencies = speedoflight_ms / ( wavelengths * length_dict[input_units] ) / frequency_dict[frequency_units]

    return frequencies

def velocity_to_wavelength(velocities, input_units, center_wavelength=None,
        center_wavelength_units=None, wavelength_units='meters',
        convention='optical'):
    """
    Conventions defined here:
    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html

     * Radio 	V = c (c/l0 - c/l)/(c/l0) 	f(V) = (c/l0) ( 1 - V/c )
     * Optical 	V = c ((c/l0) - (c/l))/(c/l) 	f(V) = (c/l0) ( 1 + V/c )^-1
     * Redshift 	z = ((c/l0) - (c/l))/(c/l) 	f(V) = (c/l0) ( 1 + z )-1
     * Relativistic 	V = c ((c/l0)^2 - (c/l)^2)/((c/l0)^2 + (c/l)^2) 	f(V) = (c/l0) { 1 - (V/c)2}1/2/(1+V/c)
    """
    if input_units in wavelength_dict:
        print("Already in wavelength units (%s)" % input_units)
        return velocities
    if center_wavelength is None:
        raise ValueError("Cannot convert velocity to wavelength without specifying a central wavelength.")
    if wavelength_units not in wavelength_dict:
        raise ValueError("Bad wavelength units: %s" % (wavelength_units))

    velocity_ms = velocities / velocity_dict['m/s'] * velocity_dict[input_units]
    center_frequency = speedoflight_ms / center_wavelength
    if convention == 'radio':
        wav = (velocity_ms / speedoflight_ms - 1.0) * center_frequency * -1
    elif convention == 'optical':
        wav = (velocity_ms / speedoflight_ms + 1.0)**-1 * center_frequency
    elif convention == 'relativistic':
        wav = (-(velocity_ms / speedoflight_ms)**2 + 1.0)**0.5 / (1.0 + velocity_ms/speedoflight_ms) * center_frequency
    else:
        raise ValueError('Convention "%s" is not allowed.' % (convention))
    wavelengths = wav / wavelength_dict[wavelength_units] * wavelength_dict[center_wavelength_units]
    return wavelengths

def wavelength_to_velocity(wavelengths, input_units, center_wavelength=None,
        center_wavelength_units=None, velocity_units='m/s',
        convention='optical'):
    """
    Conventions defined here:
    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html

     * Radio 	V = c (c/l0 - c/l)/(c/l0) 	f(V) = (c/l0) ( 1 - V/c )
     * Optical 	V = c ((c/l0) - f)/f 	f(V) = (c/l0) ( 1 + V/c )^-1
     * Redshift 	z = ((c/l0) - f)/f 	f(V) = (c/l0) ( 1 + z )-1
     * Relativistic 	V = c ((c/l0)^2 - f^2)/((c/l0)^2 + f^2) 	f(V) = (c/l0) { 1 - (V/c)2}1/2/(1+V/c)
    """
    if input_units in velocity_dict:
        print("Already in velocity units (%s)" % input_units)
        return wavelengths
    if center_wavelength is None:
        raise ValueError("Cannot convert wavelength to velocity without specifying a central wavelength.")
    if center_wavelength_units not in wavelength_dict:
        raise ValueError("Bad wavelength units: %s" % (center_wavelength_units))
    if velocity_units not in velocity_dict:
        raise ValueError("Bad velocity units: %s" % (velocity_units))

    wavelength_m = wavelengths / wavelength_dict['meters'] * wavelength_dict[input_units]
    center_wavelength_m = center_wavelength / wavelength_dict['meters'] * wavelength_dict[center_wavelength_units]
    frequency_hz = speedoflight_ms / wavelength_m
    center_frequency_hz = speedoflight_ms / center_wavelength_m

    # the order is very ugly because otherwise, if scalar, the spectroscopic axis attributes won't be inherited
    if convention == 'radio':
        velocity = ( frequency_hz - center_frequency_hz ) / center_frequency_hz * speedoflight_ms * -1
    elif convention == 'optical':
        velocity = ( frequency_hz - center_frequency_hz ) / frequency_hz * speedoflight_ms * -1
    elif convention == 'relativistic':
        velocity = ( frequency_hz**2 - center_frequency_hz**2 ) / ( center_frequency_hz**2 + frequency_hz )**2 * speedoflight_ms * -1
    else:
        raise ValueError('Convention "%s" is not allowed.' % (convention))
    velocities = velocity * velocity_dict['m/s'] / velocity_dict[velocity_units]

    return velocities

class InconsistentTypeError(Exception):
    pass

def parse_veldef(veldef):
    """
    Try to parse an 8-character FITS veldef
    """

    vconv_dict = {'OPTI':'optical','RADI':'radio','RELA':'relativistic'}

    vconv = veldef[:4]
    velocity_convention = vconv_dict[vconv]

    frame = veldef[4:]
    frame_type = frame_dict[frame]

    return velocity_convention,frame_type

