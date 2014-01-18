"""
Module with some high level utility functions.

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# December 15, 2011

import os
import re
import numpy as np
import scipy as sp
import pandas as pd
import tables as tb
import datetime as dt
import scipy.spatial as sl
import scipy.ndimage as ni
import mpl_toolkits.basemap as bm
from scipy.interpolate import UnivariateSpline as Spline


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


def spline(x, y, x_eval=None, weights=None, smooth=None):
    """Smoothing spline fit.

    The fit and smoothness depend on (i) number of samples, (ii) variance and
    (iii) error of each data point.

    Parameters
    ----------
    x, y : array-like
        The data to fit: y(x).
    x_eval : array-like, optional
        The points to evaluate the fitted spline. If not given, then the spline
        is evaluated on the same original x's.
    weights : array-like
        Array with 1/std of each data point. Note: to reflect the 
        heterokedasticity use 1/moving-window-std as weight values.
    smooth : float, optional
        Smoothing parameter. If weights are given, then s = len(weights).

    """
    ind, = np.where((~np.isnan(x)) & (~np.isnan(y)))
    x2, y2 = x[ind], y[ind]
    if weights is not None:
        weights = weights[ind]
    if x_eval is None:
        x_eval = x
    return Spline(x2, y2, w=weights, s=smooth)(x_eval)


def spline2d(time, arr3d, window=4, min_pts=10):
    """Weighted smoothing spline fit of 2d time series (3d array)."""
    _, ny, nx = arr3d.shape
    splines = np.empty_like(arr3d) * np.nan
    for i in range(ny):
        for j in range(nx):
            ts = arr3d[:,i,j]
            ind, = np.where(~np.isnan(ts))
            if len(ind) >= min_pts:
                w = 1 / pd.rolling_std(pd.Series(ts, index=time), window,
                                       min_periods=2).values
                w[0] = w[1]
                splines[:,i,j] = spline(time, ts, weights=w)
    return splines


def gradient2d(arr3d, dx=1, min_pts=10):
    """Gradient of 2d time series (3d array)."""
    _, ny, nx = arr3d.shape
    grad = np.empty_like(arr3d) * np.nan
    for i in range(ny):
        for j in range(nx):
            ts = arr3d[:,i,j]
            ind, = np.where(~np.isnan(ts))
            if len(ind) >= min_pts:
                grad[:,i,j] = np.gradient(ts, dx)
    return grad


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
    """Check if a PyTables Array can be loaded in memory."""
    if get_size(data) > max_size:
        msg = 'data is larger than %d MB, not loading in-memory!' \
            % max_size
        raise MemoryError(msg)
    else:
        return True


def need_to_save(data, max_size=128):
    """Check when data in memory need to be flushed to disk."""
    data = np.asarray(data)
    if get_size(data) > max_size:
        return True
    else:
        return False


def get_season(year, month, return_month=2):
    """Apply `_get_season()` to a scalar or sequence. 
    
    See function `_get_season()`.
    """
    if not np.iterable(year) or not np.iterable(month):
        year = np.asarray([year])
        month = np.asarray([month])
    ym = np.asarray([_get_season(y, m, return_month) for y, m in zip(year, month)])
    return ym[:,0], ym[:,1]


def _get_season(year, month, return_month=2):
    """Return year/month for a season block.
    
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
    """Generate a box given the region coords: (L,R,T,B)."""
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


def first_valid_index(arr):
    """Return index for first non-null value."""
    ind, = np.where(~np.isnan(arr))
    if len(ind) > 0:
        return ind[0]
    else:
        return None


def find_nearest(arr, val):
    """Find index for "nearest" value.
    
    Parameters
    ----------
    arr : array_like
        The array to search in (nd). No need to be sorted.
    val : scalar or array_like
        Value(s) to find.

    Returns
    -------
    out : tuple
        The index (tuple) of nearest entry found. 
        If `val` is a list of values then a tuple of ndarray
        with the indices of each value is return.

    See also
    --------
    find_nearest2

    """
    if np.ndim(val) == 0:  # scalar
        val = np.array([val]) 
    idx = []
    for v in val:
        idx.append((np.abs(arr-v)).argmin())
    idx = np.unravel_index(idx, arr.shape)
    return idx


def find_nearest2(x, y, points):
    """Find nearest x/y coords to given points.
    
    Finds the indexes of nearest coords in the 2d x and y arrays 
    to the given list of points. It searches the nearest-neighbours 
    using a k-d tree.

    Parameters
    ----------
    x, y : 2d arrays (m,n)
        Arrays containing the spatial coordinates.
    points : 2d array_like (k,2)
        List of points to find, e.g., list of tuples.

    Returns
    -------
    out : tuple
        The indices (tuple of ndarray) of the nearest entries found. 

    See also
    --------
    find_nearest

    """
    xy = np.column_stack([x.ravel(), y.ravel()]) # 2x(m,n) -> (mxn,2)
    kdtree = sl.cKDTree(xy)               # construct k-d tree
    dist, indices = kdtree.query(points)  # search points in k-d tree
    indices = np.unravel_index(indices, x.shape)
    return indices


def referenced(x, to='first'):
    """Reference a series to its first (non-null) or mean value.
    
    Parameters
    ----------
    x : array-like
        Series to be referenced.
    to : str, default 'first'
        To reference the series to its 'first' or 'mean'.
    """
    assert to in ['first', 'mean'], "`to` must be 'first' or 'mean'"
    if (not isinstance(x, np.ndarray)) and (not isinstance(x, pd.Series)):
        x = np.asarray(x)
    if np.alltrue(np.isnan(x)):
        pass
    elif to == 'first':
        x -= x[first_valid_index(x)]
    else:
        x -= x[~np.isnan(x)].mean()
    return x


def referenced2d(arr3d, to='first'):
    """Reference a 2d series to its first (non-null) or mean values."""
    nt, ny, nx = arr3d.shape
    for i in range(ny):
        for j in range(nx):
            ts = arr3d[:,i,j]
            if not np.alltrue(np.isnan(ts)):
                arr3d[:,i,j] = referenced(ts, to=to)
    return arr3d


def filter_std(arr, n=3, per_field=False):
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


def get_mask(arr):
    """Get mask from 3d array by summing-up the 0 axis.

    output: 0 = invalid region, 1 = valid region.
    """
    mask = np.nansum(arr, axis=0)
    mask[:] = ~np.isnan(mask)
    return mask


def polyfit2d(time, arr3d, deg=2, min_pts=5):
    """Least squares polynomial fit of 2d time series (3d array)."""
    _, ny, nx = arr3d.shape
    poly = np.empty_like(arr3d) * np.nan
    for i in range(ny):
        for j in range(nx):
            ts = arr3d[:,i,j]
            ind, = np.where(~np.isnan(ts))
            if len(ind) >= min_pts:
                coef = np.polyfit(time[ind], ts[ind], deg=deg)
                poly[:,i,j] = np.polyval(coef, time)
    return poly


def polyderiv(x, y, deg=3):
    """Derivative of a fitted polynomial (x,y) of degree N."""
    coeff = np.polyfit(x, y, deg=deg)
    coeff *= np.arange(deg+1)[::-1]   # [N,N-1,..,0]
    return np.polyval(coeff[:-1], x)  # N coefficients


def polyderiv2d(time, arr3d, deg=3):
    """Derivative of fitted polynomials (x,y) of degree N on a 3d array."""
    nt, ny, nx = arr3d.shape
    deriv = np.empty_like(arr3d) * np.nan
    for i in range(ny):
        for j in range(nx):
            ts = arr3d[:,i,j]
            if not np.alltrue(np.isnan(ts)):
                deriv[:,i,j] = polyderiv(time, ts, deg=deg)
    return deriv


def gfilter(field, sigma):
    """Gaussian filter (smooth) a 2d field."""
    ind = np.where(np.isnan(field))
    field[ind] = 0
    field = np.c_[field[:,-1], field, field[:,0]]  # add crossed borders!
    field = ni.gaussian_filter(field, sigma, order=0)
    field = field[:,1:-1]                          # exclude borders
    field[ind] = np.nan
    return field


def gfilter2d(arr3d, sigma):
    """Gaussian filter of 2d time series (3d array)."""
    ind = np.where(np.isnan(arr3d))
    arr3d[ind] = 0
    for k, field in enumerate(arr3d):
        field = np.c_[field[:,-1], field, field[:,0]]  # add borders
        field = ni.gaussian_filter(field, sigma, order=0)
        arr3d[k] = field[:,1:-1]                       # exclude borders
    arr3d[ind] = np.nan
    return arr3d


def regrid2d(arr3d, x, y, inc_by=2):
    """Regrid 2d time series (3d array) increasing resolution."""
    nt, ny, nx = arr3d.shape
    out = np.empty((nt, inc_by * ny, inc_by * nx), 'f8')
    xi = np.linspace(x.min(), x.max(), inc_by * len(x))
    yi = np.linspace(y.min(), y.max(), inc_by * len(y))
    xx, yy = np.meshgrid(xi, yi)
    arr3d = np.ma.masked_invalid(arr3d)
    for k, field in enumerate(arr3d):
        field1 = bm.interp(field, x, y, xx, yy, order=0) # nearest neighb.
        field2 = bm.interp(field, x, y, xx, yy, order=1) # linear
        ######## soemthing "wierd" when the field is zero ########
        ind = np.where(field2 == 0) #<<<<< check!
        try:
            field2[ind] = field1[ind]
        except:
            pass
        ##########################################################
        out[k] = field2
    return [out, xx, yy]


def human_order(text):
    """Usage: list.sort(key=human_order), sorts strings in human order."""
    atoi = lambda text: int(text) if text.isdigit() else text
    return [atoi(c) for c in re.split('(\d+)', text)]


def read_climate_index(fname, from_year=1992, to_year=2012, missing_val=-9999,
                       comments='#'):
    """Load ascii climate-index table into a time series y(x) -> [x,y]."""
    table = np.loadtxt(fname, comments=comments)
    table[table==missing_val] = np.nan
    x = np.arange(table[0,0], table[-1,0]+1, 1/12.) 
    y = table[:,1:].flatten()
    ind, = np.where((x >= from_year) & (x <= to_year))
    return [x[ind], y[ind]]


def u2l(m, k=0, mult=1):
    """Coppy the upper triangle to the lower triangle.

    Parameters
    ----------
    m : array_like, shape (M, N)
        Input array.
    k : int, optional
        Diagonal above which to copy. `k = 0` (the default) is the main
        diagonal, `k < 0` is below it and `k > 0` is above.
    mult : float, optional
        Multiplying factor for the triangle to be copied.

    See Also
    --------
    l2u : same thing, but copy lower to upper triangle. 

    """
    nrows, ncols = m.shape
    for i in range(nrows):
        for j in range(i+k, ncols):
            m[j][i] = mult * m[i][j]
    return m


def l2u(m, k=0, mult=1):
    """Coppy the lower triangle to the upper triangle.

    Parameters
    ----------
    m : array_like, shape (M, N)
        Input array.
    k : int, optional
        Diagonal above which to copy. `k = 0` (the default) is the main
        diagonal, `k < 0` is below it and `k > 0` is above.
    mult : float, optional
        Multiplying factor for the triangle to be copied.

    See Also
    --------
    u2l : same thing, but copy upper to lower triangle. 

    """
    nrows, ncols = m.shape
    for j in range(ncols):
        for i in range(j+k, ncols):
            m[j][i] = mult * m[i][j]
    return m
