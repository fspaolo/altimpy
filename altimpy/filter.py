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


def std_filt(arr, n=3, per_field=False):
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


def med_filt(arr, size=(3,3), min_pixels=3, **kw):
    """Median filter with constrain for 2d array. 
    
    Supports NaNs.
    Uses a minimum number of non-null pixels.
    """
    def _median(x, min_pixels=min_pixels):  # not an ordinary median!
        central_pixel = x[(len(x)-1)/2]
        valid_pixels = x[~np.isnan(x)]
        if np.isnan(central_pixel) or len(valid_pixels) < min_pixels:
            m = central_pixel
        else:
            m = np.median(valid_pixels)
        return m
    return ni.generic_filter(arr, _median, size=size, **kw)


def hp_filt(y, lamb=7, nan=False):
    """Hodrick-Prescott filter for 1d array.
    
    Supports NaNs, nan=True.
    Assumes an evenly spaced array.
    """
    if nan:
        y2 = y.copy()
        i_isnan, = np.where(np.isnan(y))
        i_notnan, = np.where(~np.isnan(y))
        x = np.arange(len(y))
        y2[i_isnan] = np.interp(i_isnan, x[i_notnan], y[i_notnan])
        y2 = sm.tsa.filters.hpfilter(y2, lamb=lamb)[1]
        y2[i_isnan] = np.nan
    else:
        y2 = sm.tsa.filters.hpfilter(y, lamb=lamb)[1]
    return y2


def time_filt(t, y, from_time=1991, to_time=2013):
    """Filter an array based on time interval."""
    k, = np.where((t > from_time) & (t < to_time))
    return t[k], y[k,...]


def perc_filt(x, min_perc=0.7):
    """Filter vector with a min percentage of non-null entries."""
    if np.isnan(x).all():
        pass 
    elif ( len(x[~np.isnan(x)]) / float(len(x)) ) < min_perc:
        x[:] = np.nan
    return x


def step_filt(x, delta=3, window=7):
    """Filter step-changes in a verctor.
    
    Detects level-shifts in a time series and corrects them by levelling both
    sides of the record. 
    
    Discriminates steps from peaks using a moving-median approach.
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


def _peak_filt(x, y, n_std=3):
    """Filter spikes in a vector.
    
    See peak_filt()
    """
    assert not np.isnan(y).all(), 'empty array'
    y2 = y.copy()
    i_notnan, = np.where(~np.isnan(y))
    # detrend
    p = ap.lasso_cv(x[i_notnan], y[i_notnan], max_deg=3)
    y2[i_notnan] = y[i_notnan] - p
    # filter
    i_peaks, = np.where(np.abs(y2) > n_std * np.nanstd(y2))
    y[i_peaks] = np.nan
    return len(i_peaks)


def peak_filt(x, y, n_std=3, iterative=True):
    """Filter spikes in a vector (iteratively).

    Remove values greater than n*std from the trend.
    """
    n_peaks = 0
    if not np.isnan(y).all():
        n_peaks = _peak_filt(x, y, n_std=n_std)
        if iterative and n_peaks != 0:
            while n_peaks != 0 and not np.isnan(y).all():
                n_peaks = _peak_filt(x, y, n_std=n_std)
    return y


