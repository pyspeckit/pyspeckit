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
#.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
#Affiliation: University of Colorado at Boulder
#Created on March 9th, 2011


import numpy as np
import warnings

# declare a case-insensitive dict class to return case-insensitive versions of each dictionary...
# this is just a shortcut so that units can be specified as, e.g., Hz, hz, HZ, hZ, etc.  3 of those 4 are "legitimate".  
# the templates for this code were grabbed from a few StackOverflow.com threads
class CaseInsensitiveDict(dict):
    def __init__(self, inputdict=None):
        if inputdict:
            # Doesn't do keyword args
            if isinstance(inputdict, dict):
                self.CaseDict = inputdict
                for k,v in inputdict.items():
                    if hasattr(k,'lower'):
                        dict.__setitem__(self, k.lower(), v)
            else:
                for k,v in inputdict:
                    if hasattr(k,'lower'):
                        dict.__setitem__(self, k.lower(), v)

    def __getitem__(self, key):
        return dict.__getitem__(self, key.lower())

    def __setitem__(self, key, value):
        dict.__setitem__(self, key.lower(), value)

    def __contains__(self, key):
        if hasattr(key,'lower'):
            return dict.__contains__(self, key.lower())
        else: # e.g., None will go here
            return dict.__contains__(self, key)

    def has_key(self, key):
        """ This is deprecated, but we're keeping it around """
        return dict.has_key(self, key.lower())

    def get(self, key, def_val=None):
        return dict.get(self, key.lower(), def_val)

    def setdefault(self, key, def_val=None):
        return dict.setdefault(self, key.lower(), def_val)

    def update(self, inputdict):
        for k,v in inputdict.items():
            dict.__setitem__(self, k.lower(), v)

    def fromkeys(self, iterable, value=None):
        d = CaseInsensitiveDict()
        for k in iterable:
            dict.__setitem__(d, k.lower(), value)
        return d

    def pop(self, key, def_val=None):
        return dict.pop(self, key.lower(), def_val)
    

length_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'angstroms':1e-10,'A':1e-10,
        }
length_dict = CaseInsensitiveDict(length_dict)

wavelength_dict = length_dict # synonym

frequency_dict = {
        'Hz':1.0,
        'kHz':1e3,
        'MHz':1e6,
        'GHz':1e9,
        'THz':1e12,
        }
frequency_dict = CaseInsensitiveDict(frequency_dict)

velocity_dict = {'meters/second':1.0,'m/s':1.0,
        'kilometers/s':1e3,'km/s':1e3,'kms':1e3,
        'centimeters/s':1e-2,'cm/s':1e-2,'cms':1e-2,
        }
velocity_dict = CaseInsensitiveDict(velocity_dict)

pixel_dict = {'pixel':1,'pixels':1}
pixel_dict = CaseInsensitiveDict(pixel_dict)

conversion_dict = {
        'VELOCITY':velocity_dict,  'Velocity':velocity_dict,  'velocity':velocity_dict,  'velo': velocity_dict, 'VELO': velocity_dict,
        'LENGTH':length_dict,      'Length':length_dict,      'length':length_dict, 
        'WAVELENGTH':length_dict,      'WAVELength':length_dict,      'WAVElength':length_dict, 
        'FREQUENCY':frequency_dict,'Frequency':frequency_dict,'frequency':frequency_dict, 'freq': frequency_dict, 'FREQ': frequency_dict,
        'pixels':pixel_dict,'PIXELS':pixel_dict,
        }
conversion_dict = CaseInsensitiveDict(conversion_dict)

unit_type_dict = {
    'Hz' :'frequency', 'kHz':'frequency', 'MHz':'frequency', 'GHz':'frequency',
    'HZ' :'frequency', 'KHZ':'frequency', 'MHZ':'frequency', 'GHZ':'frequency',
    'hz' :'frequency', 'khz':'frequency', 'mhz':'frequency', 'ghz':'frequency',
    'THz':'frequency', 
    'meters/second':'velocity', 'm/s':'velocity', 'kilometers/s':'velocity',
    'km/s':'velocity', 'kms':'velocity', 'centimeters/s':'velocity',
    'cm/s':'velocity', 'cms':'velocity', 
    'meters':'wavelength','m':'wavelength',
    'centimeters':'wavelength','cm':'wavelength',
    'millimeters':'wavelength','mm':'wavelength',
    'nanometers':'wavelength','nm':'wavelength',
    'micrometers':'wavelength','micron':'wavelength','microns':'wavelength','um':'wavelength',
    'kilometers':'wavelength','km':'wavelength',
    'angstroms':'wavelength','A':'wavelength',
    'unknown':'pixels',
    None: 'pixels',
    }

unit_type_dict = CaseInsensitiveDict(unit_type_dict)

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

frame_dict = CaseInsensitiveDict({
        'VLSR':'LSRK','VRAD':'LSRK','VELO':'LSRK',
        'VOPT':'LSRK',
        'VELO-LSR':'LSRK',
        'LSRD':'LSRD',
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
        'Velocity':'VELO','Frequency':'FREQ','Wavelength':'WAVE','Length':'WAVE','Redshift':'REDS'}
convention_suffix = {'radio':'RAD','optical':'OPT','relativistic':'REL','redshift':'RED'}

speedoflight_ms = 2.99792458e8 # m/s
speedoflight_kms = 2.99792458e5 # km/s
speedoflight_cms = 2.99792458e10 # cm/s

class SpectroscopicAxis(np.ndarray):
    """
    A Spectroscopic Axis object to store the current units of the spectrum and
    allow conversion to other units and frames.  Typically, units are velocity,
    wavelength, frequency, or redshift.  Wavenumber is also hypothetically
    possible.

    WARNING: If you index a SpectroscopicAxis, the resulting array will be
    a SpectroscopicAxis without a dxarr attribute!  This can result in major problems;
    a workaround is being sought but subclassing numpy arrays is harder than I thought
    """

    def __new__(self, xarr, unit="Hz", frame='rest', frame_offset=0.0,
            frame_offset_units='Hz', xtype=None, refX=None, redshift=None,
            refX_units=None):
        """
        Make a new spectroscopic axis instance
        Default units Hz
        
        *xarr* [ np.ndarray ]
            An array of X-axis values in whatever unit specified

        *unit* [ string ]
            Any valid spectroscopic X-axis unit (km/s, Hz, angstroms, etc.)

        *xtype* [ string | None ]
            irrelevant?

        *refX* [ float ]
            Reference frequency/wavelength

        *refX_units* [ string ]
            Units of the reference frequency/wavelength
        
        *frame* [ frame_dict key ]
            **NOT IMPLEMENTED**
            The frame of the axis, e.g. 
                * Topocentric (observatory-center) [up to 0.5 km/s offset from Earth-center]
                * Geocentric (earth-center) 
                * Barycentric (Solar System center) 
                * Heliocentric (Sun-center)
                * Kinematic Local Standard of Rest (LSRK)
                * Dynamic Local Standard of Rest (LSRD)
                * Galactocentric
                * Supergalactic
                * CMB?
                * "Rest" (relative to the rest frequency of the line)
                * Redshifted
            Unfortunately, there don't appear to be any implementations of this
            available in python or wrapped in python.  I'll work on incorporating
            GBTIDL's http://www.gb.nrao.edu/GBT/DA/gbtidl/release/user/toolbox/frame_velocity.html
            but it's not easy!
            http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
            BUT THIS MIGHT BE IT: https://github.com/scottransom/pyslalib
        """
        subarr = np.array(xarr)
        subarr = subarr.view(self)
        if unit in unit_type_dict:
            subarr.units = unit
        else:
            raise ValueError('Unit %s not recognized.' % unit)
        if subarr.units is None:
            subarr.units = 'none'
        subarr.frame = frame
        if xtype in xtype_dict:
            subarr.xtype = xtype_dict[xtype]
            subarr.frame = frame_dict[xtype]
        elif subarr.units in unit_type_dict:
            subarr.xtype = unit_type_dict[subarr.units]
        elif type(xtype) is str:
            warnings.warn("Unknown X-axis type in header: %s" % xtype)
            subarr.xtype = xtype
        else:
            if xtype is not None:
                warnings.warn("WARNING: xtype has been specified but was not recognized as a string.")
            subarr.xtype = 'unknown'
        subarr.refX = refX
        if refX_units is None:
            if subarr.units in frequency_dict:
                subarr.refX_units = subarr.units
            else:
                subarr.refX_units = 'Hz'
        else:
            subarr.refX_units = refX_units
        subarr.redshift = redshift
        subarr.wcshead = {}
        if subarr.xtype is not None:
            if 'RAD' in subarr.xtype:
                subarr.velocity_convention = 'radio'
            elif 'OPT' in subarr.xtype:
                subarr.velocity_convention = 'optical'
            elif 'REL' in subarr.xtype:
                subarr.velocity_convention = 'relativistic'
            elif subarr.xtype is 'unknown':
                subarr.velocity_convention = 'radio' # default
            else:
                subarr.velocity_convention = 'radio' # default
        else:
            subarr.velocity_convention = 'radio' # default

        return subarr

    def __array_finalize__(self,obj):
        """
        Array finalize allows you to take subsets of the array but keep units
        around, e.g.:
        xarr = self[1:20]
        """
        self.units = getattr(obj, 'units', None)
        self.frame = getattr(obj, 'frame', None)
        self.xtype = getattr(obj, 'xtype', None)
        self.refX = getattr(obj, 'refX', None)
        self.refX_units = getattr(obj, 'refX_units', None)
        self.velocity_convention = getattr(obj, 'velocity_convention', None)
        self.redshift = getattr(obj, 'redshift', None)
        self.wcshead = getattr(obj, 'wcshead', None)
        # moved from __init__ - needs to be done whenever viewed
        # (this is slow, though - may be better not to do this)
        if self.shape: # check to make sure non-scalar
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

    def __repr__(self):
        if self.shape is ():
            rep = ("SpectroscopicAxis([%r], units=%r, refX=%r, refX_units=%r, frame=%r, redshift=%r, xtype=%r)" %
                (self.__array__(), self.units, self.refX, self.refX_units, self.frame, self.redshift, self.xtype))
        else:
            rep = ("SpectroscopicAxis([%r,...,%r], units=%r, refX=%r, refX_units=%r, frame=%r, redshift=%r, xtype=%r)" %
                (self[0], self[-1], self.units, self.refX, self.refX_units, self.frame, self.redshift, self.xtype))
        return rep

    def __str__(self):
        selfstr =  "SpectroscopicAxis with units %s and range %g:%g." % (
                self.units,self.umin(),self.umax())
        if self.refX is not None:
            selfstr += "Reference is %g %s" % (self.refX, self.refX_units)
        return selfstr

    def _check_consistent_type(self):
        """
        Make sure self.xtype is unit_type_dict[units]
        if this is NOT true, can cause significant errors
        """
        OK = self.xtype == unit_type_dict[self.units]
        if not OK:
            raise InconsistentTypeError("Units: %s Type[units]: %s in %r" % (self.units,unit_type_dict[self.units],self) )

    """ OBSOLETE use convert_to_unit
    def change_xtype(self,new_xtype,**kwargs):
        if new_xtype not in conversion_dict:
            raise ValueError("There is no X-axis type %s" % new_xtype)
        if conversion_dict[self.xtype] is velocity_dict:
            if conversion_dict[new_xtype] is frequency_dict:
                self.velocity_to_frequency(**kwargs)
            elif conversion_dict[new_xtype] is wavelength_dict:
                self.velocity_to_frequency()
                self.frequency_to_wavelength(**kwargs)
        elif conversion_dict[self.xtype] is frequency_dict:
            if conversion_dict[new_xtype] is velocity_dict:
                self.frequency_to_velocity(**kwargs)
            elif conversion_dict[new_xtype] is wavelength_dict:
                self.frequency_to_wavelength(**kwargs)
        elif conversion_dict[self.xtype] is wavelength_dict:
            if conversion_dict[new_xtype] is velocity_dict:
                self.wavelength_to_frequency()
                self.frequency_to_velocity(**kwargs)
            elif conversion_dict[new_xtype] is wavelength_dict:
                self.wavelength_to_frequency(**kwargs)
        else:
            raise ValueError("Conversion to xtype %s was not recognized." % xtype)
    """

    def umax(self, units=None):
        """
        Return the maximum value of the SpectroscopicAxis.  If units specified,
        convert to those units first
        """

        if units is not None:
            # could be reversed
            return np.max([self.coord_to_x(self.max(), units),self.coord_to_x(self.min(),units)])
        else: 
            return self.max()

    def umin(self, units=None):
        """
        Return the minimum value of the SpectroscopicAxis.  If units specified,
        convert to those units first
        """

        if units is not None:
            # could be reversed
            return np.min([self.coord_to_x(self.max(), units),self.coord_to_x(self.min(),units)])
        else: 
            return self.min()

    def x_to_pix(self, xval):
        """
        Given an X coordinate in SpectroscopicAxis' units, return the corresponding pixel number
        """
        nearest_pix = np.argmin(np.abs(self-xval))
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
        xarr.units = 'km/s'
        xarr.refX = 5.0 
        xarr.refX_units = GHz
        xarr.x_to_coord(5.1,'GHz') == 6000 # km/s
        """
        if xunit in wavelength_dict: # shortcut - convert wavelength to freq
            xval = wavelength_to_frequency(xval, xunit, 'Hz')
            xunit = 'Hz'

        if unit_type_dict[self.units] == unit_type_dict[xunit]:
            if verbose: print "Input units of same type as output units."
            return xval * conversion_dict[unit_type_dict[xunit]][xunit] / conversion_dict[unit_type_dict[self.units]][self.units]
        elif xunit in velocity_dict:
            if verbose: print "Requested units are Velocity"
            if self.units in frequency_dict:
                return velocity_to_frequency(xval, xunit,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        frequency_units=self.units,
                        convention=self.velocity_convention)
            elif self.units in wavelength_dict:
                FREQ = velocity_to_frequency(xval, xunit,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        frequency_units='Hz',
                        convention=self.velocity_convention)
                return frequency_to_wavelength(FREQ, 'Hz',
                        wavelength_units=self.units)
        elif xunit in frequency_dict:
            if verbose: print "Requested units are Frequency"
            if self.units in velocity_dict:
                return frequency_to_velocity(xval, xunit,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        velocity_units=self.units,
                        convention=self.velocity_convention)
            elif self.units in wavelength_dict:
                return frequency_to_wavelength(xval, xunit,
                        wavelength_units=self.units)
        else:
            raise ValueError("Units not recognized.")

    def coord_to_x(self, xval, xunit):
        """
        Given an X-value assumed to be in the coordinate axes, return that
        value converted to xunit
        e.g.:
        xarr.units = 'km/s'
        xarr.refX = 5.0 
        xarr.refX_units = GHz
        xarr.coord_to_x(6000,'GHz') == 5.1 # GHz
        """
        if unit_type_dict[self.units] == unit_type_dict[xunit]:
            return xval / conversion_dict[unit_type_dict[xunit]][xunit] * conversion_dict[unit_type_dict[self.units]][self.units]


        if xunit in velocity_dict:
            if self.units in frequency_dict:
                return frequency_to_velocity(xval, self.units,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        velocity_units=xunit,
                        convention=self.velocity_convention)
            elif self.units in wavelength_dict:
                FREQ = wavelength_to_frequency(xval, self.units, frequency_units='Hz')
                return frequency_to_velocity(FREQ, 'Hz',
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        velocity_units=xunit,
                        convention=self.velocity_convention)
        elif xunit in frequency_dict:
            if self.units in velocity_dict:
                return velocity_to_frequency(xval, self.units,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        frequency_units=xunit,
                        convention=self.velocity_convention)
            elif self.units in wavelength_dict:
                return wavelength_to_frequency(xval, self.units, frequency_units=xunit)
        elif xunit in wavelength_dict:
            if self.units in velocity_dict:
                FREQ = velocity_to_frequency(xval, self.units,
                        center_frequency=self.refX,
                        center_frequency_units=self.refX_units,
                        frequency_units='Hz',
                        convention=self.velocity_convention)
                return frequency_to_wavelength(FREQ, 'Hz', wavelength_units=xunit)
            elif self.units in frequency_dict:
                return frequency_to_wavelength(xval, self.units, wavelength_units=xunit)
        else:
            raise ValueError("Units not recognized.")

    def convert_to_unit(self, unit, **kwargs):
        """
        Return the X-array in the specified units without changing it
        Uses as_unit for the conversion, but changes internal values rather
        than returning them.
        """
        self[:] = self.as_unit(unit, **kwargs)
        
        if unit in velocity_dict:
            self.xtype = "Velocity"
        elif unit in frequency_dict:
            self.xtype = "Frequency"
        elif unit in wavelength_dict:
            self.xtype = "Wavelength"

        self.units = unit
        self.dxarr = np.diff(self)

    def as_unit(self, unit, frame=None, quiet=True, center_frequency=None,
            center_frequency_units=None, **kwargs):
        """
        Convert the spectrum to the specified units.  This is a wrapper function
        to convert between frequency/velocity/wavelength and simply change the 
        units of the X axis.  Frame conversion is... not necessarily implemented.

        *unit* [ string ] 
            What unit do you want to 'view' the array as?

        *frame* [ None ]
            NOT IMPLEMENTED.  When it is, it will allow you to convert between
            LSR, topocentric, heliocentric, rest, redshifted, and whatever other
            frames we can come up with.  Right now the main holdup is finding a 
            nice python interface to an LSR velocity calculator... and motivation.
 
        *center_frequency* [ None | float ]
        *center_frequency_units* [ None | string ]
            If converting between velocity and any other spectroscopic type,
            need to specify the central frequency around which that velocity is
            calculated.
            I think this can also accept wavelengths....
        """

        if unit is None:
            return self
        elif unit in pixel_dict:
            return np.arange(self.shape[0])
        elif unit not in conversion_dict[self.xtype]:
            change_xtype = True
            change_units = True
        elif unit != self.units: 
            change_xtype = False
            change_units = True
        else: 
            change_xtype = False
            change_units = False
        if frame is not None and frame != self.frame: change_frame = True
        else: change_frame = False

        if center_frequency is None:
            center_frequency = self.refX
        if center_frequency_units is None:
            center_frequency_units = self.refX_units

        if change_xtype:
            if unit in velocity_dict:
                if conversion_dict[self.xtype] is frequency_dict:
                    newxarr = frequency_to_velocity(self,self.units,
                            center_frequency=center_frequency,
                            center_frequency_units=center_frequency_units,
                            velocity_units=unit, convention=self.velocity_convention)
                    newxtype = "Velocity"
                    newunit = unit
                elif conversion_dict[self.xtype] is wavelength_dict:
                    freqx = wavelength_to_frequency(self, self.units) 
                    cf = wavelength_to_frequency(center_frequency, center_frequency_units)
                    cfu = 'GHz'
                    newxarr = frequency_to_velocity(freqx, 'GHz',
                            center_frequency=cf,
                            center_frequency_units=cfu,
                            velocity_units=unit, convention=self.velocity_convention)
                    newxtype = "Velocity"
                    newunit = unit
            elif unit in frequency_dict:
                if conversion_dict[self.xtype] is velocity_dict:
                    newxarr = velocity_to_frequency(self,self.units,
                            center_frequency=center_frequency,
                            center_frequency_units=center_frequency_units,
                            frequency_units=unit, convention=self.velocity_convention)
                    newxtype = "Frequency"
                    newunit = unit
                elif conversion_dict[self.xtype] is wavelength_dict:
                    newxarr = wavelength_to_frequency(self, self.units,
                            frequency_units=unit)
                    newxtype = "Frequency"
                    newunit = unit
            elif unit in wavelength_dict:
                if conversion_dict[self.xtype] is velocity_dict:
                    freqx = velocity_to_frequency(self, self.units,
                            center_frequency=center_frequency,
                            center_frequency_units=center_frequency_units,
                            frequency_units=unit, convention=self.velocity_convention)
                    newxarr = frequency_to_wavelength(freqx, 'Hz',
                            wavelength_units=unit)
                    newxtype = "Wavelength"
                    newunit = unit
                elif conversion_dict[self.xtype] is frequency_dict:
                    newxarr = frequency_to_wavelength(self, self.units,
                            wavelength_units=unit)
                    newxtype = "Wavelength"
                    newunit = unit
            else:
                warnings.warn("Could not convert from %s to %s" % (self.units,unit))
        else:
            newxtype = self.xtype
            newxarr = self
            newunit = self.units

        # re-check whether units need to be changed; it is possible that change_xtype left you 
        # with the correct units
        if unit != newunit and change_units:
            conversion_factor = conversion_dict[newxtype][newunit] / conversion_dict[newxtype][unit] 
            if not quiet: print "Converting units from %s to %s" % (newunit,unit)
            newxarr = newxarr*conversion_factor

        if change_frame and not quiet:
            if not quiet: print "Conversion from frame %s to %s is not yet supported" % (self.frame,frame)

        if not change_units and not change_xtype and not change_frame and not quiet:
            if not quiet: print "Already in desired units, X-type, and frame"

        if not quiet:
            print "Converted to %s (%s)" % (newxtype,newunit)

        # this should be implemented but requires a callback to spectrum...
        #if replot:
        #    self.spectrum.plotter(reset_xlimits=True)

        return newxarr

    def change_frame(self, frame):
        """
        Change velocity frame
        """
        if self.frame == frame:
            return
        
        if frame in frame_type_dict:
            self[:] = self.in_frame(frame)
            self.dxarr = np.diff(self)
            self.frame = frame

    def make_dxarr(self, coordinate_location='center'):
        """
        Create a "delta-x" array corresponding to the X array.

        *coordinate_location* [ 'left', 'center', 'right' ]
            Does the coordinate mark the left, center, or right edge of the
            pixel?  If 'center' or 'left', the *last* pixel will have the same
            dx as the second to last pixel.  If right, the *first* pixel will
            have the same dx as the second pixel.

        """

    def in_frame(self, frame):
        """
        Return a shifted xaxis
        """
        if self.frame in ('obs','observed') and frame in ('rest',):
            if self.redshift is not None:
                return self/(1.0+self.redshift)
            else:
                warnings.warn("WARNING: Redshift is not specified, so no shift will be done.")
        elif self.frame in ('rest',) and frame in ('obs','observed'):
            if self.redshift is not None:
                return self*(1.0+self.redshift)
            else:
                warnings.warn("WARNING: Redshift is not specified, so no shift will be done.")
        else:
            print "Frame shift from %s to %s is not defined or implemented." % (self.frame, frame)
            return self

    def x_in_frame(self, xx, frame):
        """
        Return the value 'x' shifted to the target frame
        """
        return self.in_frame(frame)[self.x_to_pix(xx)]

    def cdelt(self, tolerance=1e-8):
        """
        Return the channel spacing if channels are linear
        """
        if not hasattr(self,'dxarr'): # if cropping happens...
            self.dxarr = np.diff(self)
        if abs(self.dxarr.max()-self.dxarr.min())/abs(self.dxarr.min()) < tolerance:
            return self.dxarr.mean().flat[0]

    def _make_header(self, tolerance=1e-8):
        """
        Generate a set of WCS parameters for the X-array
        """
        self.dxarr = np.diff(self)

        if self.wcshead is None:
            self.wcshead = {}

        self.wcshead['CUNIT1'] = self.units
        if fits_type[self.xtype] == 'VELO' and self.velocity_convention is not None:
            ctype = 'V'+convention_suffix[self.velocity_convention]
        else:
            ctype = fits_type[self.xtype]
        self.wcshead['CTYPE1'] = ctype+fits_frame[self.frame]
        self.wcshead['SPECSYS'] = fits_specsys[self.frame]
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

        units = axislist[0].units
        xtype = axislist[0].xtype
        frame = axislist[0].frame
        redshift = axislist[0].redshift
        velocity_convention = axislist[0].velocity_convention
        for ax in axislist:
            if ax.xtype != xtype:
                try: 
                    ax.change_xtype(xtype)
                except:
                    ValueError("Axis had wrong xtype and could not be converted.")
            if ax.units != units or ax.frame != frame:
                try:
                    ax.convert_to_unit(units,frame=frame)
                except:
                    ValueError("Axis had wrong units and could not be converted.")

        subarr = np.concatenate([ax for ax in axislist])
        subarr = subarr.view(self)
        subarr.units = units
        subarr.xtype = xtype
        subarr.frame = frame
        subarr.redshift = redshift
        subarr.velocity_convention = velocity_convention

        # if all the spectra have the same reference frequency, there is one common refX
        # else, refX is undefined and velocity transformations should not be done
        refXs_diff = np.sum([axislist[0].refX != ax.refX for ax in axislist])
        if refXs_diff > 0:
            subarr.refX = None
            subarr.refX_units = None
        else:
            subarr.refX = axislist[0].refX
            subarr.refX_units = axislist[0].refX_units

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

        units = axislist[0].units
        xtype = axislist[0].xtype
        frame = axislist[0].frame
        redshift = axislist[0].redshift
        velocity_convention = axislist[0].velocity_convention
        for ax in axislist:
            if ax.xtype != xtype:
                try: 
                    ax.change_xtype(xtype)
                except:
                    ValueError("Axis had wrong xtype and could not be converted.")
            if ax.units != units or ax.frame != frame:
                try:
                    ax.convert_to_unit(units,frame=frame)
                except:
                    ValueError("Axis had wrong units and could not be converted.")

        subarr = np.array(axislist)
        subarr = subarr.view(self)
        subarr.units = units
        subarr.xtype = xtype
        subarr.frame = frame
        subarr.redshift = redshift
        subarr.velocity_convention = velocity_convention

        # if all the spectra have the same reference frequency, there is one common refX
        # else, refX is undefined and velocity transformations should not be done
        refXs_diff = np.sum([axislist[0].refX != ax.refX for ax in axislist])
        if refXs_diff > 0:
            subarr.refX = None
            subarr.refX_units = None
        else:
            subarr.refX = axislist[0].refX
            subarr.refX_units = axislist[0].refX_units

        return subarr

def velocity_to_frequency(velocities, input_units, center_frequency=None,
        center_frequency_units=None, frequency_units='Hz',
        convention='radio'):
    """
    Conventions defined here:
    http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
    * Radio 	V = c (f0 - f)/f0 	f(V) = f0 ( 1 - V/c )
    * Optical 	V = c (f0 - f)/f 	f(V) = f0 ( 1 + V/c )-1
    * Redshift 	z = (f0 - f)/f 	f(V) = f0 ( 1 + z )-1
    * Relativistic 	V = c (f02 - f 2)/(f02 + f 2) 	f(V) = f0 { 1 - (V/c)2}1/2/(1+V/c) 
    """
    if input_units in frequency_dict:
        print "Already in frequency units"
        return
    if center_frequency is None:
        raise ValueError("Cannot convert velocity to frequency without specifying a central frequency.")
    if frequency_units not in frequency_dict:
        raise ValueError("Bad frequency units: %s" % (frequency_units))

    velocity_ms = velocities / velocity_dict['m/s'] * velocity_dict[input_units]
    if convention == 'radio':
        freq = center_frequency * (1.0 - velocity_ms / speedoflight_ms)
    elif convention == 'optical':
        freq = center_frequency * (1.0 + velocity_ms / speedoflight_ms) - 1.0
    elif convention == 'relativistic':
        freq = center_frequency * (1.0 - (velocity_ms / speedoflight_ms)**2)**0.5 / (1.0 + velocity_ms/speedoflight_ms)
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
        print "Already in velocity units"
        return
    if center_frequency is None:
        raise ValueError("Cannot convert frequency to velocity without specifying a central frequency.")
    if center_frequency_units not in frequency_dict:
        raise ValueError("Bad frequency units: %s" % (center_frequency_units))
    if velocity_units not in velocity_dict:
        raise ValueError("Bad velocity units: %s" % (velocity_units))

    frequency_hz = frequencies / frequency_dict['Hz'] * frequency_dict[input_units]
    center_frequency_hz = center_frequency / frequency_dict['Hz'] * frequency_dict[center_frequency_units]

    if convention == 'radio':
        velocity = speedoflight_ms * ( center_frequency_hz - frequency_hz ) / center_frequency_hz
    elif convention == 'optical':
        velocity = speedoflight_ms * ( frequency_hz - center_frequency_hz ) / frequency_hz
    elif convention == 'relativistic':
        velocity = speedoflight_ms * ( center_frequency_hz**2 - frequency_hz**2 ) / ( center_frequency_hz**2 + frequency_hz )**2
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
        print "Already in wavelength units"
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
        print "Already in frequency units"
        return
    if frequency_units not in frequency_dict:
        raise ValueError("Frequency units %s not valid" % wavelength_units)

    if input_units not in length_dict:
        raise AttributeError("Cannot convert from wavelength unless units are already wavelength.")

    frequencies = speedoflight_ms / ( wavelengths * length_dict[input_units] ) / frequency_dict[frequency_units]

    return frequencies

class InconsistentTypeError(Exception):
    pass
