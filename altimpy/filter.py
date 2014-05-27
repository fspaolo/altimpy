"""
Module with some high level filter functions.

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# May 20, 2014

import os
import re
import numpy as np
import pandas as pd
import scipy.ndimage as ni
import statsmodels.api as sm


def stdfilt(arr, n=3, per_field=False):
    """Filter out pts greater than n std.

    Parameters
    ----------
    arr : array-like (2d or 3d)
        Array containing the field(s) to filter, with or without NaNs.
    n : integer, optional
        Number of standar deviations to use, default n=3.
    per_field : boolean, default False
        If 'False', uses the distribution of the whole dataset. If 'True',
        filters each 2d field independently (one distribution per field).

    """
    nstd = lambda arr, n: n * np.std(arr[~np.isnan(arr)])
    npts_before = len(arr[~np.isnan(arr)])
    if per_field:
        for i, field in enumerate(arr):
            arr[i, np.abs(field) > nstd(field, n)] = np.nan
    else:
        arr[np.abs(arr) > nstd(arr, n)] = np.nan
    npts_after = float(len(arr[~np.isnan(arr)]))
    print 'filter out pts: %d (%.1f%%)' % \
          (npts_before - npts_after, 100 * (1 - npts_after/npts_before))
    return arr


def medfilt(arr, size=(3,3), min_pixels=3, **kw):
    """Median filter 2d array. 
    
    It handleds NaNs.
    It uses a minimum number of valid pixels.
    """
    def median(x, min_pixels=min_pixels):
        central_pixel = x[(len(x)-1)/2]
        valid_pixels = x[~np.isnan(x)]
        if np.isnan(central_pixel) or len(valid_pixels) < min_pixels:
            m = central_pixel
        else:
            m = np.median(valid_pixels)
        return m
    return ni.generic_filter(arr, median, size=size, **kw)


def hpfilt(y, lamb=7):
    """Hodrick-Prescott filter 1d array."""
    return sm.tsa.filters.hpfilter(y, lamb=lamb)[1]


def timefilt(t, y, from_time=1991, to_time=2013):
    """Filter an array based on time interval."""
    k, = np.where((t > from_time) & (t < to_time))
    return t[k], y[k,...]


def percfilt(x, min_perc=0.7):
    """Filter vector with a min percentage of non-null entries."""
    if np.isnan(x).all():
        pass 
    elif ( len(x[~np.isnan(x)]) / float(len(x)) ) < min_perc:
        x[:] = np.nan
    return x


def stepfilt(x, delta=3, window=7):
    """Filter step-changes in a verctor.
    
    Detects level-shifts in a time series and corrects them by levelling both
    sides of the record. 
    
    It discrimitaes steps from peaks using a moving-median approach.
    """
    assert window % 2 != 0, 'window size must be odd'
    n = window / 2
    v = np.r_[np.repeat(x[0], n), x, np.repeat(x[-1], n)] # expanded arr
    m = pd.rolling_median(v, window, center=True)         # filtered arr
    for i in range(len(m)-1):
        diff = m[i+1] - m[i]
        if np.abs(diff) > delta:
            v[i+1:] -= diff
    return v[n:-n], m[n:-n]
 

