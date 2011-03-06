"""
obstools.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-01-13.

Description: Container for functions useful in data mining / observational work.
"""

def RAtoDeg(ra):
    """
    Convert right acsension from hhmmss.ss format to degrees.
    
        ra: 
            Type: str
            Right ascension in hhmmss.ss format
            Ex: ra = '123456.78'
    """
    
    ret = float(ra[0:2]) * 360. / 24.
    ret += float(ra[2:4]) * 360. / 24. / 60.
    ret += float(ra[4:]) * 360. / 24. / 3600.
    
    return ret
    
def RAtoStr(ra):
    """
    Convert right ascension from degrees to hhmmss.ss format. 
    
       ra: 
            Type: float
            Right ascension in degrees
    """
    
    hrs = ra * 24. / 360.
    mins = (hrs % 1.) * 60.
    secs = (mins % 1.) * 60.
    
    hh = "{0:02g}".format(int(hrs))
    mm = "{0:02g}".format(int(mins))
    ss = ('%.02f' % (round(secs, 2))).zfill(5)
    
    return "{0}{1}{2}".format(hh, mm, ss)
        
def DECtoDeg(dec):
    """
    Convert declination from ddmmss.ss format to degrees.
    
        dec: 
            Type: str
            Declination in ddmmss.ss format
            Ex: dec = '+123456.78'
    """
    
    if dec[0] == '+': ret = float(dec[1:3])
    elif dec[0] == '-': ret = -float(dec[1:3])
    else: print 'Declination must include +/- as first character!'
    
    if ret > 0:
        ret += float(dec[3:5]) / 60.
        ret += float(dec[5:]) / 3600.
    if ret < 0:
        ret -= float(dec[3:5]) / 60.
        ret -= float(dec[5:]) / 3600.
        
    return ret
    
def DECtoStr(dec):
    """
    Convert declination from degrees to ddmmss.ss format.
    
        dec: 
            Type: float
            Declination in degrees
    """
    
    deg = int(abs(dec))
    mins = (abs(dec) % 1.) * 60.
    secs = (mins % 1) * 60.
    
    dd = "{0:02g}".format(int(deg))
    mm = "{0:02g}".format(int(mins))
    ss = ('%.02f' % (round(secs, 2))).zfill(5)
    
    if dec < 0: return "-{0}{1}{2}".format(dd, mm, ss)
    else: return "+{0}{1}{2}".format(dd, mm, ss)
    
        
        
        