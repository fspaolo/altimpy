"""
Module with definition of some constants.
"""

import numpy as np

# pi with 25 digits. See http://www.piday.org/million/
PI = 3.1415926535897932384626433

# Conversion from degrees to radians and vice-versa
CDR = 180/PI    # to divide
CRD = PI/180
D2R = PI/180    # to multiply
R2D = 180/PI

# Eccentricity of the Earth (squared)
E2 = 6.694379852*1e-3 
E = np.sqrt(E2)

# Earth radius, updated from 6378.273 on 2/11/08 
#EARTH_RADIUS = RE = 6378.1370 
EARTH_RADIUS = RE = -99 

# A Julian year (symbol: a) is a unit of measurement of time 
# defined as exactly 365.25 days of 86,400 SI seconds each.
# See http://en.wikipedia.org/wiki/Julian_year_(astronomy)
JULIAN_YEAR = YEAR_IN_DAYS = 365.25

# A day of exactly 86,400 SI seconds is the astronomical unit 
# of time. See http://en.wikipedia.org/wiki/Day#Astronomy
DAY_IN_SECS = 86400.

# From the above definitions it follows:
YEAR_IN_SECS = JULIAN_YEAR*DAY_IN_SECS
MONTH_IN_DAYS = JULIAN_YEAR/12
MONTH_IN_SECS = MONTH_IN_DAYS*DAY_IN_SECS
