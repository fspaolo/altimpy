"""
Module with some high level utility functions.

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# December 15, 2011

import os
import re
import numpy as np
import scipy as sp
import tables as tb
import datetime as dt


class CircularList(list):
    """A list that wraps around instead of throwing an index error.

    Simple, perhaps incomplete implementation of a Circular List in 
    Python. Subclasses list and overrides __getitem__. The only 
    special behavior is that attempts to access indices which are out 
    of bounds will wrap around - accessing mylist(len(mylist)) should 
    return the first item in the list instead of an IndexError Slice 
    operations must still be 'in bounds' First tries list's 
    __getitem__. If that is not successful, it converts the index key 
    to an integer, then calculates the appropriate 'in bounds' index 
    and returns whatever is stored there. If converting the key to 
    integer fails, TypeError is raised.
    
    Works like a regular list:
    >>> cl = CircularList([1,2,3])
    >>> cl
    [1, 2, 3]
    >>> cl[0]
    1
    >>> cl[-1]
    3
    >>> cl[2]
    3

    Except wraps around:
    >>> cl[3]
    1
    >>> cl[-4]
    3
    
    Slices work
    >>> cl[0:2]
    [1, 2]
    
    but only in range.
    """
    def __getitem__(self, key):
        # try normal list behavior
        try:
            return super(CircularList, self).__getitem__(key)
        except IndexError:
            pass
        # key can be either integer or slice object,
        # only implementing int now.
        try:
            index = int(key)
            index = index % self.__len__()
            return super(CircularList, self).__getitem__(index)
        except ValueError:
            raise TypeError


def linear_fit(x, y, return_coef=False):
    """
    Fit a straight-line by Ordinary Least Squares.

    If `return_coef=True` returns the slope (m) and intercept (c).
    """
    ind, = np.where((~np.isnan(x)) & (~np.isnan(y)))
    if len(ind) < 2: 
        return [np.nan, np.nan]
    x, y = x[ind], y[ind]
    A = np.ones((len(x), 2))
    A[:,0] = x
    m, c = np.linalg.lstsq(A, y)[0]
    if return_coef:
        return (m, c)
    else:
        x_val = np.linspace(x.min(), x.max(), 200)
        y_fit = m*x_val + c
        return (x_val, y_fit)


def linear_fit_robust(x, y, return_coef=False):
    """
    Fit a straight-line by robust regression (M-estimate).

    If `return_coef=True` returns the slope (m) and intercept (c).
    """
    import scikits.statsmodels.api as sm
    ind, = np.where((~np.isnan(x)) & (~np.isnan(y)))
    if len(ind) < 2: 
        return [np.nan, np.nan]
    x, y = x[ind], y[ind]
    X = sm.add_constant(x, prepend=False)
    y_model = sm.RLM(y, X, M=sm.robust.norms.HuberT())
    y_fit = y_model.fit()
    if return_coef:
        if len(y_fit.params) < 2: return (y_fit.params[0], 0.)
        else: return y_fit.params[:]
    else:
        return (x, y_fit.fittedvalues)


def poly_fit(x, y, order=3, return_coef=False, npts=200):
    """
    Fit a polynomial of order `order` to data points `x,y`.
    """
    ind, = np.where((~np.isnan(x)) & (~np.isnan(y)))
    if len(ind) < 3: 
        return [np.nan, np.nan]
    x, y = x[ind], y[ind]
    coef = np.polyfit(x, y, order)
    if return_coef:
        return coef
    else:
        x_val = np.linspace(x.min(), x.max(), npts)
        y_fit = np.polyval(coef, x_val)
        return (x_val, y_fit)


def spline_interp(x, y, smooth=0.01):
    """
    Interpolate data using cubic splines of given smoothness.
    smooth : smoothness factor
    """
    from scipy.interpolate import splrep, splev
    ind, = np.where((~np.isnan(x)) & (~np.isnan(y)))
    x, y = x[ind], y[ind]
    # find the knot points
    tck = splrep(x, y, s=smooth)
    # evaluate spline on interpolated points
    x_val = np.linspace(x.min(), x.max(), 200)
    y_fit = splev(x_val, tck)
    return (x_val, y_fit)


def get_size(arr):
    """
    Get the size in MB of a Numpy or PyTables object.

    parameters
    ----------
    arr : 1D/2D Numpy or PyTables Array.
    """
    try:
        m, n = arr.shape
    except:
        m, n = arr.shape, 1
    num_elem = m*n
    item_size = arr.dtype.itemsize
    return (item_size*num_elem/1e6)


def can_be_loaded(data, max_size=512):
    """
    Check if a PyTables Array can be loaded in memory.
    """
    if get_size(data) > max_size:
        msg = 'data is larger than %d MB, not loading in-memory!' \
            % max_size
        raise MemoryError(msg)
    else:
        return True


def need_to_save(data, max_size=128):
    """
    Check when data in memory need to be flushed on disk.
    """
    data = np.asarray(data)
    if get_size(data) > max_size:
        return True
    else:
        return False


def get_season(year, month, return_month=2):
    """
    Apply `_get_season()` to a scalar or sequence. See `_get_season()`.
    """
    if not np.iterable(year) or not np.iterable(month):
        year = np.asarray([year])
        month = np.asarray([month])
    ym = np.asarray([_get_season(y, m, return_month) for y, m in zip(year, month)])
    return ym[:,0], ym[:,1]


def _get_season(year, month, return_month=2):
    """
    Returns the first, second or third month of the 3-month 
    season-block, and update the `year` when needed.

    year, month : int
    return_month : 1, 2 or 3
    """
    if return_month != 1 and return_month != 2 and return_month != 3:
        raise IOError('`return_month` must be: 1, 2 or 3')
    MAM = [3, 4, 5]      # Mar/Apr/May -> Fall SH 
    JJA = [6, 7, 8]      # Jun/Jul/Aug -> winter SH
    SON = [9, 10, 11]    # Sep/Oct/Nov -> Spring SH
    DJF = [12, 1, 2]     # Dec/Jan/Feb -> summer SH
    return_month -= 1
    if month in MAM:
        return year, MAM[return_month]
    elif month in JJA:
        return year, JJA[return_month]
    elif month in SON:
        return year, SON[return_month]
    elif month in DJF:
        if month == 12 and return_month > 0:
            year += 1
        return year, DJF[return_month]
    else:
        print 'not a valid month from 1 to 12!'
        return None, None


def get_box(region, npts=None):
    """
    Generate a box given the region coords: (L,R,T,B).
    """
    west, east, south, north = region
    if npts:
        n = int(npts/4.)
        x = np.empty(n*4, 'f8')
        y = np.empty(n*4, 'f8')
        lons = np.linspace(west, east, n)
        lats = np.linspace(south, north, n)
        x[:n] = lons[:]
        y[:n] = north
        x[n:n*2] = east
        y[n:n*2] = lats[::-1]
        x[n*2:n*3] = lons[::-1]
        y[n*2:n*3] = south
        x[n*3:] = west 
        y[n*3:] = lats[:]
    else:
        x = np.array([west, east, east, west, west])
        y = np.array([north, north, south, south, north])
    return [x, y]


# find the n-tuples correspondent to each element of a list
ntuples = lambda lst, n: zip(*[lst[i:]+lst[:i] for i in range(n)])


# CHECK
def first_non_null(arr):
    """Return index of first non-null element."""
    return (i for i,elem in enumerate(arr) if ~np.isnan(elem)).next()


def first_non_null2(arr):
    ind, = np.where(~np.isnan(arr))
    if len(ind) > 0:
        return ind[0]
    else:
        return None


def find_nearest(arr, val, return_value=False):
    """Find index or value of array entry "nearest" to val.
    
    Parameters
    ----------
    arr : array_like
        The array to search in (nd).
    val : scalar or array_like
        Value(s) to find.
    return_value : bool, optional
        To return the actual value instead of index.
    Returns
    -------
    output : scalar or tuple of ndarray 
        The index (tuple) or value of nearest entry found. 
        If `val` is a list of values then a tuple of ndarray
        with the indices of each value is return, or a list
        with the nearest values if `return_value` is True.
    """
    shape = arr.shape
    if np.ndim(val) == 0:                 # scalar
        idx = (np.abs(arr-val)).argmin()  # index of flat array
    else:
        idx = []
        for v in val:
            idx.append((np.abs(arr-v)).argmin())
    idx = np.unravel_index(idx, shape)
    if return_value: 
        return arr[idx]
    else:
        return idx

