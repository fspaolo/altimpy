"""
Module with some high level filtering functions.

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# May 20, 2014

import os
import re
import numpy as np
import pandas as pd
import scipy.ndimage as ni
import statsmodels.api as sm
import matplotlib.pyplot as plt

from altimpy import lasso_cv


def std_filter(arr, n=3, per_field=False):
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


def median_filter(arr, size=(3,3), min_pixels=3, **kw):
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


def hp_filter(y, lamb=7, nan=False, return_series='trend'):
    """Hodrick-Prescott filter for 1d array.
    
    lamb (int): smoothing paramenter (e.g., 1600 for trend in quarterly data)
    nan=True|False: supports NaNs.
    return_series='trend'|'cycle'|'both': returns the trend and/or cycle.

    Assumes an evenly spaced array.

    Detrend with HP-filter
    ----------------------
    
    Rule of thumb for smoothing is (empirical)[1]: 
    
    Lambda = 100 * (number of periods in a year)^2
    
    Annual data = 100*1^2 = 100
    Quarterly data = 100*4^2 = 1,600
    Monthly data = 100*12^2 = 14,400
    Weekly data = 100*52^2 = 270,400
    
    There is additional research that suggests using a power of 4 instead
    of 2. See Ravn and Uhlig (2002)[2]
    
    Harvey and Trimbur (2008) explain the risk in using a too-small smoothing
    parameter (lambda), though I have yet to find research explaining the risk
    of using too-large of a smooth parameter, other than the trend becomes
    increasingly linear and less sensitive to recent data[3]
    
    [1] http://forums.eviews.com/viewtopic.php?f=4&t=870
    [2] https://ideas.repec.org/a/tpr/restat/v84y2002i2p371-375.html
    [3] http://www.terrapub.co.jp/journals/jjss/pdf/3801/38010041.pdf
        
    """
    if nan:
        y2 = y.copy()
        i_isnan, = np.where(np.isnan(y))
        i_notnan, = np.where(~np.isnan(y))
        x = np.arange(len(y))
        y2[i_isnan] = np.interp(i_isnan, x[i_notnan], y[i_notnan])
        c2, y2 = sm.tsa.filters.hpfilter(y2, lamb=lamb)  # cycle and trend
        y2[i_isnan] = np.nan
        c2[i_isnan] = np.nan
    else:
        c2, y2 = sm.tsa.filters.hpfilter(y, lamb=lamb)   # cycle and trend
    if return_series == 'trend':
        res = y2
    elif return_series == 'cycle':
        res = c2
    else:
        res = [y2, c2]  # trend and cycle
    return res


def time_filter(t, y, from_time=1991, to_time=2013):
    """Filter an array based on time interval."""
    k, = np.where((t > from_time) & (t < to_time))
    return t[k], y[k,...]


def percent_filter(y, min_perc=0.7):
    """Filter vector with a min percentage of non-null entries."""
    y2 = y.copy()
    percent = len(y2[~np.isnan(y2)]) / float(len(y2))
    if percent < min_perc:
        y2[:] = np.nan
    return y2, percent


def step_filter(x, delta=3, window=7):
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
            #plt.plot(v)
            v[i+1:] -= diff
            #plt.plot(v)
            #plt.show()
    return v[n:-n], m[n:-n]


def _peak_filter(x, y, n_std=3, max_deg=3):
    """Filter spikes in a vector.
    
    See peak_filter()
    """
    assert not np.isnan(y).all(), 'empty array'
    y2 = y.copy()
    i_notnan, = np.where(~np.isnan(y))
    # detrend
    poly = lasso_cv(x[i_notnan], y[i_notnan], max_deg=max_deg)
    y2[i_notnan] = y[i_notnan] - poly
    # filter
    i_peaks, = np.where(np.abs(y2) > n_std * np.nanstd(y2))
    y[i_peaks] = np.nan
    return len(i_peaks)


def peak_filter(x, y, n_std=3, iterative=True, max_deg=3):
    """Filter spikes in a vector (iteratively).

    Remove values greater than n*std from the trend.
    """
    n_peaks = _peak_filter(x, y, n_std=n_std, max_deg=3)
    if iterative and n_peaks != 0:
        while n_peaks != 0 and not np.isnan(y).all():
            n_peaks = _peak_filter(x, y, n_std=n_std)
    return y


def std_series_filter(x, y, max_std=2, max_deg=3):
    """Filter entire vector if the detrended std > max_std."""
    std = None
    if not np.isnan(y).all():
        i_notnan, = np.where(~np.isnan(y))
        poly = lasso_cv(x[i_notnan], y[i_notnan], max_deg=max_deg)
        std = np.nanstd(y[i_notnan] - poly) # std of detrended series
        if std > max_std:
            y[:] = np.nan
    return y, std
