from __future__ import print_function
import numpy as np
import units

def header_to_vlsr(header,**kwargs):
    """
    Attempt to determine the vlsr offset using default header keywords
    """

    if 'RA' in header:
        ra = sexagesimal_to_decimal(header['RA'],ra=True)
        dec = sexagesimal_to_decimal(header['DEC'],ra=False)
    elif 'CRVAL1' in header and 'RA' in header['CTYPE1']:
        ra = header['CRVAL1']
        dec = header['CRVAL2']
    latitude = header['LATITUDE']
    if 'UT1' in header:
        datestr =  header['UT1']
    elif 'DATE-OBS' in header:
        datestr =  header['DATE-OBS']
    else:
        print("Not sure what date keywords exist...")
    jd = date_to_JD(datestr)
    if 'EQUINOX' in header:
        epoch = header['EQUINOX']
    elif 'EPOCH' in header:
        epoch = header['EPOCH']
    else:
        epoch = 2000.

    lst = sexagesimal_to_decimal(header['LST'])

    print(ra,dec,latitude,lst,jd,epoch)
    return topo_to_lsr(ra,dec,latitude,lst,jd,epoch)

def date_to_JD(datestr, **kwargs):
    import pyslalib
    import dateutil
    date = dateutil.parser.parse(datestr, **kwargs)
    day = date.day + (date.hour/24.) + (date.minute/60.)/24. + (date.second/3600.)/24.
    JD = pyslalib.slalib.sla_cldj(date.year,date.month,day)
    return JD

def sexagesimal_to_decimal(string, delimiter=":", ra=False):
    hh,mm,ss = [float(x) for x in string.split(delimiter)]
    decim = hh+mm/60.+ss/3600.
    if ra:
        decim *= 15.
    return decim

def topo_to_geo(ra,dec,latitude,lst,height=None):
    """
    :ra: Decimal right ascension (degrees)
    :dec: Decimal declination (degrees)
    :latitude: Observatory latitude in degrees
    :lst: Local sidereal time (hours)

    .. warning: this option was implemented but returns inconsistent results
        and has therefore been disabled
    :height: (optional) Height above the reference spheroid in meters
        (if height is specified, a more precise version of the code is used)
    """
    import pyslalib
    ra_rad = ra*np.pi/180.
    dec_rad = dec*np.pi/180.
    latitude_rad = latitude*np.pi/180.
    lst_rad = lst*15.*np.pi/180.
    if height is None:
        return pyslalib.slalib.sla_rverot(ra_rad,dec_rad,latitude_rad,lst_rad)
    else:
        return pyslalib.slalib.sla_rverot(ra_rad,dec_rad,latitude_rad,lst_rad)
        vearth_6d = pyslalib.slalib.sla_pvobs(latitude_rad, height, lst_rad)
        vradec_3d = pyslalib.slalib.sla_dcs2c(ra_rad, dec_rad)
        projected_vearth = -1*pyslalib.slalib.sla_dvdv(vearth_6d[3:],vradec_3d)
        return projected_vearth * units.astronomical_unit_cm*units.length_dict['cm']/units.length_dict['km']


def geo_to_bary(ra,dec,jd,epoch):
    """
    For a given ra/dec, return the conversion from geocentric velocity to heliocentric

    :ra: Decimal right ascension (degrees)
    :dec: Decimal declination (degrees)
    :jd: Modified (2000 = 51544) julian date
    :epoch: Epoch of observations (e.g., 2000, 2008.202)

    follows instructions given here:
    http://star-www.rl.ac.uk/docs/sun67.htx/node230.html

    *  Star vector, J2000
        CALL sla_DCS2C(RM,DM,V)

    *  Earth/Sun velocity and position, J2000
        CALL sla_EVP(TDB,2000D0,DVB,DPB,DVH,DPH)

    *  Radial velocity correction due to Earth orbit (km/s)
        VCORB = -sla_DVDV(V,DVH)*149.597870D6

    """
    import pyslalib
    velocity_3d = pyslalib.slalib.sla_dcs2c(ra*np.pi/180.,dec*np.pi/180.)
    dvb, dpb, dvh, dph = pyslalib.slalib.sla_evp(jd,epoch) 

    vcorb = -pyslalib.slalib.sla_dvdv(velocity_3d,dvb)*149.597870e6

    return vcorb

def geo_to_helio(ra,dec,jd,epoch):
    """
    For a given ra/dec, return the conversion from geocentric velocity to heliocentric

    :ra: Decimal right ascension (degrees)
    :dec: Decimal declination (degrees)
    :jd: Modified (2000 = 51544) julian date
    :epoch: Epoch of observations (e.g., 2000, 2008.202)

    follows instructions given here:
    http://star-www.rl.ac.uk/docs/sun67.htx/node230.html

    *  Star vector, J2000
        CALL sla_DCS2C(RM,DM,V)

    *  Earth/Sun velocity and position, J2000
        CALL sla_EVP(TDB,2000D0,DVB,DPB,DVH,DPH)

    *  Radial velocity correction due to Earth orbit (km/s)
        VCORB = -sla_DVDV(V,DVH)*149.597870D6

    """
    import pyslalib
    velocity_3d = pyslalib.slalib.sla_dcs2c(ra*np.pi/180.,dec*np.pi/180.)
    dvb, dpb, dvh, dph = pyslalib.slalib.sla_evp(jd,epoch) 

    vcorb = -pyslalib.slalib.sla_dvdv(velocity_3d,dvh)*149.597870e6

    return vcorb

def helio_to_lsr(ra, dec):
    """
    """

    import pyslalib
    return pyslalib.slalib.sla_rvlsrk(ra/180.*np.pi, dec/180.*np.pi)

def topo_to_lsr(ra, dec, latitude, lst, jd, epoch, height=None):
    print("helio->lsr: ",helio_to_lsr(ra,dec))
    print("geo->helio: ",geo_to_helio(ra,dec,jd,epoch))
    print("topo->geo: ",topo_to_geo(ra,dec,latitude,lst,height=height))
    return helio_to_lsr(ra,dec) + geo_to_helio(ra,dec,jd,epoch) + topo_to_geo(ra,dec,latitude,lst,height=height)

def frame_grid(ra=0.0,dec=0.0,latitude=0.0,jd=51544,lst=0,vtopo=0.0,epoch=2000,):
    frame_names = ['TOPO','GEO','BARY','HELI','LSRK']
    frame1_to_frame2 = dict([(n1,None) for n1 in frame_names])
    for frame1 in frame_names:
        frame1_to_frame2[frame1] = dict([(n2,None) for n2 in frame_names])
        for frame2 in frame_names:
            if frame1 == frame2:
                frame1_to_frame2[frame1][frame2] = 0.0
            elif frame1 == 'TOPO' and frame2 == 'GEO':
                frame1_to_frame2[frame1][frame2] = topo_to_geo(ra,dec,latitude,lst)


if __name__ == "__main__":
    import pytest 
    @pytest.mark.parametrize
    def test_vlsr(ra, dec, utyr, utmon, utday, epoch, latitude, elevation):
        pass

