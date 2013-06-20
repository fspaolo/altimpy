import os
import re
import numpy as np
import scipy as sp
import tables as tb
import datetime as dt
import pandas as pd

### time conversion functions

class SecsToDateTime(object):
    """
    Converts `seconds since epoch` to `datetime` (i.e., year, month, day).

    secs : 1D array, decimal seconds.
    since_year : int, ref_epoch = <since_year>-Jan-1 00:00:00 is assumed.
    since_epoch : tuple, especifies ref_epoch as (YYYY, MM, DD, hh, mm, ss).

    Notes
    -----
    1. Matlab uses as time reference the year 0000, and Python 
       `datetime` uses the year 0001.
    2. utc85 (or ESA-time) is seconds since 1985-1-1 00:00:00,
       ICESat-time is seconds since 2000-1-1 12:00:00,
       secs00 is seconds since 2000-1-1 00:00:00.

    """
    def __init__(self, secs=0, since_year=1985, since_epoch=None):
        if np.ndim(secs) > 0:
            self.secs = np.asarray(secs)
        else:
            self.secs = secs  

        if since_epoch is None:
            # <since_year>-Jan-1 00:00:00
            ref_epoch = dt.date(since_year, 1, 1)
        else:
            # not working !!!!!!!!!!!!!!!!!!!!!
            ref_epoch = dt.datetime(since_epoch)    

        # ref_epoch in days since 0001-Jan-1 00:00:00
        ref_epoch_in_days = mpl.date2num(ref_epoch)  

        # secs/86400 -> frac days -> date
        frac_days = self.secs / (24*60*60.)
        self._datenum = ref_epoch_in_days + frac_days
        self._dates = mpl.num2date(self._datenum)

    def datenum(self, matlab=False):
        if matlab:
            # frac days since 0000-Jan-1 00:00:00
            return self._datenum + 366.
        else:
            # frac days since 0001-Jan-1 00:00:00
            return self._datenum

    def dates(self):
        return self._dates

    def years(self):
        return np.array([d.year for d in self._dates])

    def months(self):
        return np.array([d.month for d in self._dates])
        
    def days(self):
        return np.array([d.day for d in self._dates])

    def ymdhms(self):
        return np.array([(d.year, d.month, d.day, d.hour,
            d.minute, d.second) for d in self._dates])


def lon_180_360(lon, region):
    """
    Convert lon from 180 to 360 (or vice-verse) according to `region`.
    """
    l, r, b, t = region
    if l < 0:
        lon[lon>180] -= 360  # 360 -> 180
    elif r > 180:
        lon[lon<0] += 360    # 180 -> 360
    return lon


def ll2xy(lon, lat, slat=71, slon=0, hemi='s', units='km'):
    """
    Convert from 'lon/lat' to polar stereographic 'x/y'.
 
    This function converts from geodetic latitude and longitude to
    polar stereographic 'x/y' coordinates for the polar regions. The 
    equations are from Snyder, J.P., 1982, Map Projections Used by 
    the U.S. Geological Survey, Geological Survey Bulletin 1532, U.S. 
    Government Printing Office. See JPL Technical Memorandum 
    3349-85-101 for further details.
    
    Parameters
    ----------
    lon, lat : array_like (rank-1 or 2) or float 
        Geodetic longitude and latitude (degrees, -/+180 or 0/360 and -/+90).
    slat : float
        Standard latitude (e.g., 71 S), see Notes.
    slon : float
        Standard longitude (e.g., -70), see Notes.
    hemi : string
        Hemisphere: 'n' or 's' (not case-sensitive).
    units : string
        Polar Stereographic x,y units: 'm' or 'km' (not case-sensitive).
    
    Returns
    -------
    x, y : ndarray (rank-1 or 2) or float
        Polar stereographic x and y coordinates (in 'm' or 'km').

    Notes
    -----
    SLAT is is the "true" latitude in the plane projection 
    (the map), so there is no deformation over this latitude; 
    e.g., using the same SLON but changing SLAT from 70 to 71 
    degrees, will move things in polar stereo. The goal is to 
    locally preserve area and angles. Most users use 71S but 
    the sea ice people use 70S.
    
    SLON provides a "vertical" coordinate for plotting and for 
    rectangle orientation. E.g., for Arctic sea ice, NSIDC use 
    SLON=45 in order to make a grid that is optimized for where 
    sea ice occurs. CATS2008a has SLON=-70 (AP roughly up), so 
    that the grid can be long enough to include South Georgia.

    MOA Image Map (the GeoTIFF): SLAT=-71, SLON=0
    MOA mask grid (from Laurie): SLAT=-71, SLON=-70
    Scripps mask grid (from GeoTIFF): SLAT=-71, SLON=0

    History
    -------
    Written in Fortran by C.S. Morris - Apr 29, 1985
    Revised by C.S. Morris - Dec 11, 1985
    Revised by V.J. Troisi - Jan 1990
        SGN - provides hemisphere dependency (+/- 1)
    Revised by Xiaoming Li - Oct 1996
        Corrected equation for RHO
    Converted to Matlab by L. Padman - Oct 25, 2006
    Updated for slon by L. Padman - Nov 21, 2006
    Converted to Python by F.S. Paolo - Mar 23, 2010
    
    Example
    -------
    >>> lon = [-150.3, 66.2, 5.3]
    >>> lat = [70.2, 75.5, 80.3]
    >>> x, y = ll2xy(lon, lat, slat=71, slon=-70, hemi='s', units='m')

    Original (Matlab) documentation
    -------------------------------
    ARGUMENTS:                                                         
                                                                       
    Variable     I/O    Description                          
                                                                        
    lat           I     Geodetic Latitude (degrees, +90 to -90)
    lon           I     Geodetic Longitude (degrees, 0 to 360)
    SLAT          I     Standard latitude (typ. 71, or 70)
    SLON          I  
    HEMI          I     Hemisphere (char*1: 'N' or 'S' (not
                                    case-sensitive)
    x             O     Polar Stereographic X Coordinate (km)
    y             O     Polar Stereographic Y Coordinate (km)

    """
    if units != 'm':
        units = 'km'
    print 'parameters:'
    print 'standard lat:', slat
    print 'standard lon:', slon
    print 'hemisphere:', hemi
    print 'units of x,y:', units
    print "converting lon,lat -> x,y ..."

    # if sequence convert to ndarray double
    if type(lon).__name__ in ['list', 'tuple']:
        lon = np.array(lon, 'f8') 
        lat = np.array(lat, 'f8')        

    # if ndarray convert to double if it isn't
    if type(lon).__name__ == 'ndarray' and lon.dtype != 'float64':
        lon = lon.astype(np.float64)
        lat = lat.astype(np.float64)
 
    # convert longitude
    if type(lon).__name__ == 'ndarray':  # is numpy array
        lon[lon<0] += 360.               # -/+180 -> 0/360
    elif lon < 0:                        # is scalar
        lon += 360.                    
 
    if (str.lower(hemi) == 's'):
        SGN = -1
    else:
        SGN = 1
    if (np.abs(slat) == 90):
        RHO = 2. * RE / ((1 + E)**(1 + E) * (1 - E)**(1 - E))**(E/2.)
    else:
        SL  = np.abs(slat) / CDR
        TC  = np.tan(PI/4. - SL/2.) / ((1 - E * np.sin(SL)) \
            / (1 + E * np.sin(SL)))**(E/2.)
        MC  = np.cos(SL) / np.sqrt(1 - E2 * (np.sin(SL)**2))
        RHO = RE * MC / TC
 
    lat = np.abs(lat) / CDR
    T = np.tan(PI/4. - lat/2.) / ((1 - E * np.sin(lat)) \
      / (1 + E * np.sin(lat)))**(E/2.)
 
    lon2 = -(lon - slon) / CDR
    x = -RHO * T * np.sin(lon2)  # global vars
    y =  RHO * T * np.cos(lon2)

    if units == 'm':            # computations are done in km
        x *= 1000.
        y *= 1000.

    print 'done.'
    return [x, y]

 
def xy2ll(x, y, slat=71, slon=0, hemi='s', units='km'):
    """
    Convert from polar stereographic 'x,y' to 'lon,lat'.
 
    This subroutine converts from Polar Stereographic 'x,y' coordinates 
    to geodetic longitude and latitude for the polar regions. The 
    equations are from Snyder, J.P., 1982, Map Projections Used by the 
    U.S. Geological Survey, Geological Survey Bulletin 1532, U.S. 
    Government Printing Office. See JPL Technical Memorandum 
    3349-85-101 for further details.  
 
    Parameters
    ----------
    x, y : array_like (rank-1 or 2) or float
        Polar stereographic x and y coordinates (in 'm' or 'km').
    slat : float
        Standard latitude (e.g., 71 S), see Notes.
    slon : float
        Standard longitude (e.g., -70), see Notes.
    hemi : string
        Hemisphere: 'n' or 's' (not case-sensitive).
    units : string
        Polar Stereographic x,y units: 'm' or 'km' (not case-sensitive).
 
    Returns
    -------
    lon, lat : ndarray (rank-1 or 2) or float
        Geodetic longitude and latitude (degrees, 0/360 and -/+90).
 
    Notes
    -----
    SLAT is is the "true" latitude in the plane projection 
    (the map), so there is no deformation over this latitude; 
    e.g., using the same SLON but changing SLAT from 70 to 71 
    degrees, will move things in polar stereo. The goal is to 
    locally preserve area and angles. Most users use 71S but 
    the sea ice people use 70S.
    
    SLON provides a "vertical" coordinate for plotting and for 
    rectangle orientation. E.g., for Arctic sea ice, NSIDC use 
    SLON=45 in order to make a grid that is optimized for where 
    sea ice occurs. CATS2008a has SLON=-70 (AP roughly up), so 
    that the grid can be long enough to include South Georgia.

    MOA Image Map (the GeoTIFF): SLAT=-71, SLON=0
    MOA mask grid (from Laurie): SLAT=-71, SLON=-70
    Scripps mask grid (from GeoTIFF): SLAT=-71, SLON=0

    History
    -------
    Written in Fortran by C.S. Morris - Apr 29, 1985
    Revised by C.S. Morris - Dec 11, 1985
    Revised by V.J. Troisi - Jan 1990
        SGN - provides hemisphere dependency (+/- 1)
    Converted to Matlab by L. Padman - Oct 25, 2006
    Updated for slon by L. Padman - Nov 21, 2006
    Converted to Python by F.S. Paolo - Mar 23, 2010
 
    Example
    -------
    >>> x = [-2141.06767831  1096.06628549  1021.77465469]
    >>> y = [  365.97940112 -1142.96735458   268.05756254]
    >>> lon, lat = xy2ll(x, y, slat=71, slon=-70, hemi='s', units='km')

    Original (Matlab) documentation
    -------------------------------
    ARGUMENTS:                                                           
                                                                         
    Variable     I/O    Description                          
                                                                      
    X             I     Polar Stereographic X Coordinate (km) 
    Y             I     Polar Stereographic Y Coordinate (km)
    SLAT          I     Standard latitude (typ. 71, or 70)
    SLON          I     Standard longitude
    HEMI          I     Hemisphere (char*1, 'S' or 'N', 
                                    not case-sensitive)
    lat           O     Geodetic Latitude (degrees, +90 to -90)
    lon           O     Geodetic Longitude (degrees, 0 to 360) 

    """
    if units != 'm':
        units = 'km'
    print 'parameters:'
    print 'standard lat:', slat
    print 'standard lon:', slon
    print 'hemisphere:', hemi
    print 'units of x,y:', units
    print "converting 'x,y' -> 'lon,lat' ..."

    # if sequence convert to ndarray
    if type(x).__name__ in ['list', 'tuple']:
        x = np.array(x, 'f8')
        y = np.array(y, 'f8')
    
    # if ndarray convert to double if it isn't
    if type(x).__name__ == 'ndarray' and x.dtype != 'float64':
        x = x.astype(np.float64)
        y = y.astype(np.float64)
 
    if units == 'm':      # computations are done in km !!!
        x *= 1e-3
        y *= 1e-3

    if(str.lower(hemi) == 's'):
        SGN = -1.
    else:
        SGN = 1.
    slat = np.abs(slat)
    SL = slat / CDR
    RHO = np.sqrt(x**2 + y**2)    # if scalar, is numpy.float64
    if np.alltrue(RHO < 0.1):     # Don't calculate if "all points" on the equator
        lat = 90.0 * SGN
        lon = 0.0
        return lon, lat
    else:
        CM = np.cos(SL) / np.sqrt(1 - E2 * (np.sin(SL)**2))
        T = np.tan((PI/4.) - (SL/2.)) / ((1 - E * np.sin(SL)) \
            / (1 + E * np.sin(SL)))**(E/2.)
        if (np.abs(slat - 90.) < 1.e-5):
            T = ((RHO * np.sqrt((1 + E)**(1 + E) * (1 - E)**(1 - E))) / 2.) / RE
        else:
            T = RHO * T / (RE * CM)
 
        a1 =  5 * E2**2 / 24.
        a2 =  1 * E2**3 / 12.
        a3 =  7 * E2**2 / 48.
        a4 = 29 * E2**3 / 240.
        a5 =  7 * E2**3 / 120.
        
        CHI = (PI/2.) - 2. * np.arctan(T)
        lat = CHI + ((E2/2.) + a1 + a2) * np.sin(2. * CHI) \
              + (a3 + a4) * np.sin(4. * CHI) + a5 * np.sin(6. * CHI)
        lat = SGN * lat * CDR
        lon = -(np.arctan2(-x, y) * CDR) + slon

        #lon = SGN * (np.arctan2(SGN * x, -SGN * y) * CDR) + slon # original !!!
        #lon[lon<-180] += 360; lon[lon>180] -= 360                # keep lon to -/+180
    print 'done.'
    return [lon, lat]


### time conversion function

def sec2dt(secs, since_year=1985):
    """
    Convert seconds since_year to datetime objects.
    secs : float [array], decimal seconds.
    """
    dt_ref = dt.datetime(since_year, 1, 1, 0, 0)
    return [dt_ref + dt.timedelta(seconds=s) for s in secs]


def year2dt(year):
    """
    Convert decimal year to `datetime` object.
    year : float [array], decimal years.
    """
    if not np.iterable(year):
        year = np.asarray([year])
    ym = np.asarray([year2ym(y) for y in year])
    dt = np.asarray([pd.datetime(y, m, 15) for y, m in ym])
    return dt


def num2dt(times):
    """
    Convert a numeric representation of time to datetime.
    times : int/float array_like representing YYYYMMDD.
    """
    return np.asarray([pd.datetime.strptime(str(int(t)), '%Y%m%d') for t in times])


def ym2dt(year, month):
    """
    Convert year and month to `datetime` object.
    year, month : int arrays.
    """
    return [dt.datetime(y, m, 15) for y, m in zip(year, month)]


# NOT SURE THIS FUNC IS OK. REVIEW THE ALGORITHM!
def num2year(iyear):
    """Integer representation of year to decimal year."""
    if not np.iterable(iyear):
        iyear = np.asarray([iyear])
    iyear = np.asarray([int(y) for y in iyear])
    fyear = lambda y, m, d: y + (m - 1)/12. + d/365.25
    ymd = [int2ymd(iy) for iy in iyear]
    return [fyear(y,m,d) for y,m,d in ymd]


# NOT SURE THIS FUNC IS OK. REVIEW THE ALGORITHM!
def ym2year(year, month):
    """Year, month -> decimal year."""
    year = np.asarray(year)
    month = np.asarray(month)
    fyear = year + (month - 1)/12. + 15.22/365.25  # decimal years (month-centered)
    return fyear 


def year2ym(fyear):
    """Decimal year -> year, month."""
    fy, y  = np.modf(fyear)
    m, y = int(np.ceil(fy*12)), int(y)
    if (m == 0): m = 1
    return [y, m]


def num2ymd(iyear):
    f, y = np.modf(iyear/10000.)
    d, m = np.modf(f*100)
    return (int(y), int(m), int(d*100))


def year2num(year, day=15):
    """Decimal year to integer representation -> YYYMMDD."""
    if not np.iterable(year):
        year = np.asarray([year])
    ym = [year2ym(y) for y in year]
    return [int(y*10000 + m*100 + day) for y,m in ym]


