#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
================================
Module with conversion functions
================================

A note on 'time' vs 'datetime' modules
--------------------------------------

The time module was intended to match the functionality of the C standard 
library time.h kit, and named accordingly. The datetime module came much later.

The time module is actually more used for file dates and times. That would 
explain the epoch and 1970 boundaries as there are not much files before the 
pre-PC era to keep timestamps for. The functions in this module do not handle 
dates and times before the epoch or far in the future

The datetime module is more suited to general data processing than the time 
module.

Tips on using python 'datetime' module
--------------------------------------

http://www.enricozini.org/2009/debian/using-python-datetime/

"""
print(__doc__)

import os
import re
import time
import datetime as dt
import numpy as np
import scipy as sp
import tables as tb
import pandas as pd

from const import *

### time conversion functions

# DEPRECATED?
class SecsToDateTime(object):
    """Seconds since epoch -> datetime (year, month, day).

    secs : 1D array, decimal seconds.
    since_year : int, ref_epoch = <since_year>-Jan-1 00:00:00 is assumed.
    since_epoch : tuple, especifies ref_epoch as (YYYY, MM, DD, hh, mm, ss).

    Notes
    -----
    1. Matlab uses as time reference the year 0000, and Python 
       datetime uses the year 0001.
    2. utc85 (or ESA-time) is seconds since 1985-1-1 00:00:00 UTC,
       ICESat-time is seconds since 2000-1-1 12:00:00 UTC,
       secs00 is seconds since 2000-1-1 00:00:00 UTC.

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


#--------------------------------------------
#       coordinates transformation
#--------------------------------------------


def lon_180_360(lon, region=None, inverse=False):
    """longitude -/+180 -> 0/360 (or vice-versa). 
    
    Converts according `region` if given, otherwise converts
    from 180 to 360 if `inverse` is `False` or from 360 to 180 
    if `inverse` is True.
    """
    if region:
        l, r, b, t = region
        if l < 0:
            lon[lon>180] -= 360  # 360 -> 180
        elif r > 180:
            lon[lon<0] += 360    # 180 -> 360
    elif inverse:
        lon[lon>180] -= 360
    else:
        lon[lon<0] += 360
    return lon


def ll2xy(lon, lat, slat=71, slon=0, hemi='s', units='km'):
    """ Spherical lon/lat -> Polar Steregraphic x/y.
 
    This function converts from geodetic latitude and longitude to
    polar stereographic 'x/y' coordinates for the polar regions. The 
    equations are from Snyder, J.P., 1982, Map Projections Used by 
    the U.S. Geological Survey, Geological Survey Bulletin 1532, U.S. 
    Government Printing Office. See JPL Technical Memorandum 
    3349-85-101 for further details.
    
    Parameters
    ----------
    lon, lat : array-like (1d or 2d) or float 
        Geodetic longitude and latitude (degrees, -/+180 or 0/360 and -/+90).
    slat : float
        Standard latitude (e.g., 71 S), see Notes.
    slon : float
        Standard longitude (e.g., -70), see Notes.
    hemi : string
        Hemisphere: 'n' or 's' (not case-sensitive).
    units : string
        Polar Stereographic x/y units: 'm' or 'km' (not case-sensitive).
    
    Returns
    -------
    x, y : ndarray (1d or 2d) or float
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
    sea ice occurs. Ex: CATS2008a has SLON=-70 (AP roughly up), 
    so that the grid can be long enough to include South Georgia.
    Other examples are:

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
    print 'units of x/y:', units
    print "converting lon/lat -> x/y ..."

    # if sequence, convert to ndarray double
    if type(lon).__name__ in ['list', 'tuple']:
        lon = np.array(lon, 'f8') 
        lat = np.array(lat, 'f8')        

    # if ndarray, convert to double if it isn't
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
    """Polar Stereographic x/y -> Spherical lon/lat.
 
    This subroutine converts from Polar Stereographic 'x,y' coordinates 
    to geodetic longitude and latitude for the polar regions. The 
    equations are from Snyder, J.P., 1982, Map Projections Used by the 
    U.S. Geological Survey, Geological Survey Bulletin 1532, U.S. 
    Government Printing Office. See JPL Technical Memorandum 
    3349-85-101 for further details.  
 
    Parameters
    ----------
    x, y : array-like (1d or 2d) or float
        Polar stereographic x and y coordinates (in 'm' or 'km').
    slat : float
        Standard latitude (e.g., 71 S), see Notes.
    slon : float
        Standard longitude (e.g., -70), see Notes.
    hemi : string
        Hemisphere: 'n' or 's' (not case-sensitive).
    units : string
        Polar Stereographic x/y units: 'm' or 'km' (not case-sensitive).
 
    Returns
    -------
    lon, lat : ndarray (1d or 2d) or float
        Geodetic longitude and latitude (degrees, 0/360 and -/+90).
 
    Notes
    -----
    SLAT is the "true" latitude in the plane projection 
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

    # if sequence, convert to ndarray
    if type(x).__name__ in ['list', 'tuple']:
        x = np.array(x, 'f8')
        y = np.array(y, 'f8')
    
    # if ndarray, convert to double if it isn't
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


def sph2xyz(lon, lat, radius=1):
    """Spherical lon/lat[/r] -> Cartesian x/y/z (3d)."""
    lat *= D2R 
    lon *= D2R
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return [x, y, z]


#--------------------------------------------
#       time conversion functions
#--------------------------------------------

# OK
def sec2date(secs, since=(1985, 1, 1, 0, 0, 0)):
    """Seconds since epoch -> datetime object.

    Parameters
    ----------
    secs : scalar or array-like
        [Decimal] Seconds.
    since : tuple, (year, month, day, hour, min, sec)
        The reference time for the elapsed seconds. If only (year, month, day) 
        is provided, the following is assumed (year, month, day, 0, 0, 0).

    See also
    --------
    sec2year

    """
    assert len(since) in [3, 6], "'since' must be (y,m,d) or (y,m,d,h,m,s)"
    print 'elapsed seconds since', since, '-> date'
    if not np.iterable(secs):
        secs = np.asarray([secs], 'f8')
    else:
        secs = np.asarray(secs, 'f8')  # cast type for timedelta()
    if len(since) == 3:
        since = list(since) + [0, 0, 0]
    year, month, day, hour, minute, second = since
    dt_epoch = dt.datetime(year, month, day, hour, minute, second)
    dates = [dt_epoch + dt.timedelta(seconds=s) for s in secs]  # conversion
    return np.asarray(dates)


def sec2year(secs, since=(1985, 1, 1, 0, 0, 0)):
    """Seconds since epoch -> decimal year.

    Parameters
    ----------
    secs : scalar or array-like
        [Decimal] Seconds.
    since : tuple, (year, month, day, hour, min, sec)
        The reference time for the elapsed seconds. If only (year, month, day) 
        is provided, the following is assumed (year, month, day, 0, 0, 0).

    See also
    --------
    sec2date

    """
    return date2year(sec2date(secs, since=since))


def day2date(days, since=(1985, 1, 1, 0, 0, 0)):
    """Days since epoch -> datetime object.

    Parameters
    ----------
    days : scalar or array-like
        [Decimal] Days.
    since : tuple, (year, month, day, hour, min, sec)
        The reference time for the elapsed days. If only (year, month, day) 
        is provided, the following is assumed (year, month, day, 0, 0, 0).

    See also
    --------
    day2year

    """
    assert len(since) in [3, 6], "'since' must be (y,m,d) or (y,m,d,h,m,s)"
    print 'elapsed days since', since, '-> date'
    if not np.iterable(days):
        days = np.asarray([days], 'f8')
    else:
        days = np.asarray(days, 'f8')  # cast type for timedelta()
    if len(since) == 3:
        since = list(since) + [0, 0, 0]
    year, month, day, hour, minute, second = since
    dt_epoch = dt.datetime(year, month, day, hour, minute, second)
    dates = [dt_epoch + dt.timedelta(days=d) for d in days]  # conversion
    return np.asarray(dates)


def day2year(days, since=(1985, 1, 1, 0, 0, 0)):
    """Days since epoch -> decimal year.

    Parameters
    ----------
    days : scalar or array-like
        [Decimal] Days.
    since : tuple, (year, month, day, hour, min, sec)
        The reference time for the elapsed days. If only (year, month, day) 
        is provided, the following is assumed (year, month, day, 0, 0, 0).

    See also
    --------
    day2date

    """
    return date2year(day2date(days, since=since))


# OK
def num2date(dnum):
    """Date number as YYYYMMDD -> datetime object.

    Parameters
    ----------
    dnum : scalar or array-like
        Int or float representing date as YYYYMMDD.

    See also
    --------
    date2num

    """
    if not np.iterable(dnum):
        dnum = np.asarray([dnum])
    dates = [dt.datetime.strptime(str(int(t)), '%Y%m%d') for t in dnum]
    return np.asarray(dates)


# OK
def date2num(date):
    """Datetime object -> date number as YYYYMMDD.
    
    Parameters
    ----------
    date : object or array-like
        Datetime object(s).

    See also
    --------
    num2date

    """
    if not np.iterable(date):
        date = np.asarray([date])
    date = [int(''.join(d.date().isoformat().split('-'))) for d in date]
    return np.asarray(date)


# DEPRECATED
def ym2date(year, month):
    """Year and month -> datetime object.

    year, month : int array-like.
    """
    return np.asarray([dt.datetime(y, m, 15) for y, m in zip(year, month)])


# DEPRECATED
def year2ymd(yearfrac):
    """Decimal year -> year, month, day.
    
    It uses the Julian Year and defines months as 12
    equal-size blocks (will not necessarily coincide with 
    the Gregorian Calendar).

    The output is `year` and `months and days` *past* since 
    `year`. So in this case `days` is not a calendar day.
    """
    frac, year = np.modf(yearfrac)
    year = int(year)
    frac, month = np.modf(frac*12)
    month = int(month + 1)
    frac, day = np.modf(frac*MONTH_IN_DAYS)
    day = int(day)
    if day < 30: day += 1  # avoids using day 31
    return [year, month, day]


# OK
def year2num(year):
    """Decimal year -> date number as YYYMMDD.
    
    Parameters
    ----------
    year : scalar or array-like
        Decimal years.

    See also
    --------
    num2year
    """
    return date2num(year2date(year))


# OK
def num2year(dnum):
    """Date number as YYYYMMDD -> decimal year.
    
    Parameters
    ----------
    dnum : scalar or array-like
        Int or float representing date as YYYYMMDD.

    See also
    --------
    year2num
    """
    return date2year(num2date(dnum))


# OK
def year2date(year):
    """Decimal year -> datetime object.

    This method is probably accurate to within the second (or the 
    hour if daylight savings or other strange regional things are 
    in effect). It also works correctly during leap years.

    Parameters
    ----------
    year : scalar or array-like
        Decimal years.

    Notes
    -----
    The function can handle dates between years 0001--9999.

    See also
    --------
    date2year

    """
    if not np.iterable(year):
        year = np.asarray([year])
    def y2d(yearfrac):
        # returns seconds elapsed since 0001-01-01 00:00:00 local
        secs = lambda date: (date - dt.datetime(1,1,1)).total_seconds()
        frac, year = np.modf(yearfrac)
        year = int(year)
        start_of_this_year = dt.datetime(year=year, month=1, day=1)
        start_of_next_year = dt.datetime(year=year+1, month=1, day=1)
        year_duration = secs(start_of_next_year) - secs(start_of_this_year)
        year_elapsed = frac * year_duration 
        return start_of_this_year + dt.timedelta(seconds=year_elapsed)
    return np.asarray([y2d(y) for y in year])


# OK
def date2year(date):
    """Datetime object -> decimal year.

    This method is probably accurate to within the second (or the 
    hour if daylight savings or other strange regional things are 
    in effect). It also works correctly during leap years.

    Parameters
    ----------
    date : single or array-like datetime object(s)
        The input date can be any time zone.

    Notes
    -----
    Unlike the original function(1) using the `time` module, now the 
    time is platform-independent, and the function can handle dates 
    between years 0001--9999 (rather than 1900--2038).

    (1) Modified from http://stackoverflow.com/questions/6451655/
    python-how-to-convert-datetime-dates-to-decimal-years

    See also
    --------
    year2date

    """
    if not np.iterable(date):
        date = np.asarray([date])
    def d2y(date):
        # returns seconds elapsed since 0001-01-01 00:00:00 local
        secs = lambda date: (date - dt.datetime(1,1,1)).total_seconds()
        start_of_this_year = dt.datetime(year=date.year, month=1, day=1)
        start_of_next_year = dt.datetime(year=date.year+1, month=1, day=1)
        year_elapsed = secs(date) - secs(start_of_this_year)
        year_duration = secs(start_of_next_year) - secs(start_of_this_year)
        fraction = year_elapsed/year_duration
        return date.year + fraction
    return np.asarray([d2y(d) for d in date])
