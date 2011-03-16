"""
units.py

Author: Adam Ginsburg
Affiliation: University of Colorado at Boulder
Created on March 9th, 2011

Unit parsing and conversion tool.  
The SpectroscopicAxis class is meant to deal with unit conversion internally
"""


length_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'angstroms':1e-10,'A':1e-10,
        }

wavelength_dict = length_dict # synonym

frequency_dict = {
        'Hz':1.0,
        'kHz':1e3,
        'MHz':1e6,
        'GHz':1e9,
        'THz':1e12,
        }

velocity_dict = {'meters/second':1.0,'m/s':1.0,
        'kilometers/s':1e3,'km/s':1e3,'kms':1e3,
        'centimeters/s':1e-2,'cm/s':1e-2,'cms':1e-2,
        }

conversion_dict = {
        'VELOCITY':velocity_dict,  'Velocity':velocity_dict,  'velocity':velocity_dict,  'velo': velocity_dict, 'VELO': velocity_dict,
        'LENGTH':length_dict,      'Length':length_dict,      'length':length_dict, 
        'FREQUENCY':frequency_dict,'Frequency':frequency_dict,'frequency':frequency_dict, 'freq': frequency_dict, 'FREQ': frequency_dict,
        }

unit_type_dict = {
    'Hz' :'frequency', 'kHz':'frequency', 'MHz':'frequency', 'GHz':'frequency',
    'THz':'frequency', 
    'meters/second':'velocity', 'm/s':'velocity', 'kilometers/s':'velocity',
    'km/s':'velocity', 'kms':'velocity', 'centimeters/s':'velocity',
    'cm/s':'velocity', 'cms':'velocity', 
    'meters':'length','m':'length',
    'centimeters':'length','cm':'length',
    'millimeters':'length','mm':'length',
    'nanometers':'length','nm':'length',
    'micrometers':'length','micron':'length','microns':'length','um':'length',
    'kilometers':'length','km':'length',
    'angstroms':'length','A':'length',
    None: None,
    }

xtype_dict = {
        'VLSR':'velocity','VRAD':'velocity','VELO':'velocity',
        'VOPT':'velocity',
        'VHEL':'velocity',
        'VGEO':'velocity',
        'VREST':'velocity',
        'velocity':'velocity',
        'Z':'redshift',
        'FREQ':'frequency',
        'frequency':'frequency',
        'WAV':'length',
        'WAVE':'length',
        'wavelength':'length',
        }

frame_dict = {
        'VLSR':'LSR','VRAD':'LSR','VELO':'LSR',
        'VOPT':'LSR',
        'VHEL':'heliocentric',
        'VGEO':'geocentric',
        'velocity':'LSR',
        'VREST':'rest',
        'Z':'rest',
        'FREQ':'rest',
        'frequency':'rest',
        'WAV':'rest',
        'WAVE':'rest',
        'wavelength':'rest',
        }

speedoflight_ms = 2.99792458e8 # m/s

import numpy as np

class SpectroscopicAxis(np.ndarray):
    """
    A Spectroscopic Axis object to store the current units of the spectrum and
    allow conversion to other units and frames.  Typically, units are velocity,
    wavelength, frequency, or redshift.  Wavenumber is also hypothetically
    possible.
    """

    def __new__(self, xarr, unit, frame='rest', xtype=None, reffreq=None,
            redshift=None):
        subarr = np.array(xarr)
        subarr = subarr.view(self)
        subarr.units = unit
        subarr.frame = frame
        if xtype in xtype_dict:
            subarr.xtype = xtype_dict[xtype]
            subarr.frame = frame_dict[xtype]
        elif subarr.units in unit_type_dict:
            subarr.xtype = unit_type_dict[subarr.units]
        else:
            subarr.xtype = 'unknown'
        subarr.reffreq = reffreq
        subarr.redshift = redshift

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
        self.reffreq = getattr(obj, 'reffreq', None)
        self.redshift = getattr(obj, 'redshift', None)

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

    def convert_to_unit(self,unit,frame='rest',**kwargs):
        if unit == self.units and frame == self.frame:
            print "Already in desired units and frame"
        elif frame == self.frame:
            conversion_factor = conversion_dict[self.xtype][self.units] / conversion_dict[self.xtype][unit] 
            print "Converting units from %s to %s" % (self.units,unit)
            self *= conversion_factor
            self.units = unit
        else:
            print "(not actually) Converting frames from %s to %s" % (self.frame,frame)

    def velocity_to_frequency(self,center_frequency=None,frequency_units='Hz',convention='radio'):
        """
        Conventions defined here:
        http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
        * Radio 	V = c (f0 - f)/f0 	f(V) = f0 ( 1 - V/c )
        * Optical 	V = c (f0 - f)/f 	f(V) = f0 ( 1 + V/c )-1
        * Redshift 	z = (f0 - f)/f 	f(V) = f0 ( 1 + z )-1
        * Relativistic 	V = c (f02 - f 2)/(f02 + f 2) 	f(V) = f0 { 1 - (V/c)2}1/2/(1+V/c) 
        """
        if center_frequency is None and self.reffreq is None:
            raise ValueError("Cannot convert velocity to frequency without specifying a central frequency.")
        elif self.reffreq is not None:
            center_frequency = self.reffreq
        if frequency_units not in frequency_dict:
            raise ValueError("Bad frequency units: %s" % (frequency_units))

        velocity_ms = self / velocity_dict['m/s'] * velocity_dict[self.units]
        if convention == 'radio':
            freq = center_frequency * (1.0 - velocity_ms / speedoflight_ms)
        elif convention == 'optical':
            freq = center_frequency * (1.0 + velocity_ms / speedoflight_ms) - 1.0
        elif convention == 'relativistic':
            freq = center_frequency * (1.0 - (velocity_ms / speedoflight_ms)**2)**0.5 / (1.0 + velocity_ms/speedoflight_ms)
        else:
            raise ValueError('Convention "%s" is not allowed.' % (convention))
        self[:] = freq / frequency_dict[frequency_units]
        self.units = frequency_units
        self.xtype = 'Frequency'

    def frequency_to_velocity(self,center_frequency=None,center_frequency_units='Hz',velocity_units='m/s',convention='radio'):
        """
        Conventions defined here:
        http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
        * Radio 	V = c (f0 - f)/f0 	f(V) = f0 ( 1 - V/c )
        * Optical 	V = c (f0 - f)/f 	f(V) = f0 ( 1 + V/c )-1
        * Redshift 	z = (f0 - f)/f 	f(V) = f0 ( 1 + z )-1
        * Relativistic 	V = c (f02 - f 2)/(f02 + f 2) 	f(V) = f0 { 1 - (V/c)2}1/2/(1+V/c) 
        """
        if center_frequency is None and self.reffreq is None:
            raise ValueError("Cannot convert frequency to velocity without specifying a central frequency.")
        elif self.reffreq is not None:
            center_frequency = self.reffreq
        if center_frequency_units not in frequency_dict:
            raise ValueError("Bad frequency units: %s" % (frequency_units))
        if velocity_units not in velocity_dict:
            raise ValueError("Bad velocity units: %s" % (velocity_units))

        frequency_hz = self * frequency_dict['Hz'] / frequency_dict[self.units]
        center_frequency_hz = center_frequency * frequency_dict['Hz'] / frequency_dict[center_frequency_units]

        if convention == 'radio':
            velocity = speedoflight_ms * ( frequency_hz - center_frequency_hz ) / center_frequency_hz
        elif convention == 'optical':
            velocity = speedoflight_ms * ( frequency_hz - center_frequency_hz ) / frequency_hz
        elif convention == 'relativistic':
            velocity = speedoflight_ms * ( center_frequency_hz**2 - frequency_hz**2 ) / ( center_frequency_hz**2 + frequency_hz )**2
        else:
            raise ValueError('Convention "%s" is not allowed.' % (convention))
        self[:] = velocity * velocity_dict['m/s'] / velocity_dict[velocity_units]
        self.units = velocity_units
        self.xtype = 'Velocity'

    def frequency_to_wavelength(self,wavelength_units='um'):
        """
        Simple conversion from frequency to wavelength:
        lambda = c / nu
        """
        if wavelength_units not in length_dict:
            raise ValueError("Wavelength units %s not valid" % wavelength_units)

        if self.xtype not in frequency_dict:
            raise AttributeError("Cannot convert from frequency unless units are already frequency.  Current xtype is %s." % self.xtype)

        self[:] = speedoflight_ms / ( self * frequency_dict[self.units] ) / length_dict[wavelength_units]
        self.xtype = 'Wavelength'
        self.units = wavelength_units

    def wavelength_to_frequency(self,frequency_units='GHz'):
        """
        Simple conversion from frequency to wavelength:
        nu = c / lambda
        """
        if frequency_units not in frequency_dict:
            raise ValueError("Frequency units %s not valid" % wavelength_units)

        if self.xtype not in length_dict:
            raise AttributeError("Cannot convert from wavelength unless units are already wavelength.  Current xtype is %s." % self.xtype)

        self[:] = speedoflight_ms / ( self * length_dict[self.units] ) / frequency_dict[frequency_units]
        self.xtype = 'Frequency'
        self.units = frequency_units

class SpectroscopicAxes(SpectroscopicAxis):
    """
    Counterpart to Spectra: takes a list of SpectroscopicAxis's and
    concatenates them while checking for consistency and maintaining
    header parameters
    """

    def __new__(self, axislist, frame='rest', xtype=None, reffreq=None,
            redshift=None):

        if type(axislist) is not list:
            raise TypeError("SpectroscopicAxes must be initiated with a list of SpectroscopicAxis objects")

        units = axislist[0].units
        xtype = axislist[0].xtype
        frame = axislist[0].frame
        redshift = axislist[0].redshift
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

        subarr.reffreq = [ax.reffreq for ax in axislist]

        return subarr

