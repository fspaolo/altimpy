"""
Module with functions to form and handle time series. 

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# August 6, 2013 

import os
import re
import numpy as np
import scipy as sp
import tables as tb
import datetime as dt
from scipy import spatial

from const import *


def get_area_cells(grid, x, y):
    """Calculates the area of each grid-cell.

    Parameters
    ----------
    grid : 2d-array
        A rectangular grid.
    x, y : 1d-arrays
        The coordinates of the cells or nodes (edges). If x and y are
        of the same lengh as grid dimensions, then cell-centered 
        coordinates are assumed, otherwise edges are assumed.

    Returns
    -------
    out : 2d-array
        Same shape as 'grid' with the values of the area on each cell.

    """
    ny, nx = grid.shape
    # cells -> nodes
    if len(x) == nx and len(y) == ny:  
        dx, dy = x[1] - x[0], y[1] - y[0]
        x -= dx/2.
        y -= dy/2.
        x = np.append(x, x[-1]+dx)
        y = np.append(y, y[-1]+dy)
    C = EARTH_RADIUS**2*D2R  # deg -> rad
    A = np.empty_like(grid)
    for i in xrange(ny):
        for j in xrange(nx):
            lat1, lat2 = y[i]*D2R, y[i+1]*D2R
            lon1, lon2 = x[j]*D2R, x[j+1]*D2R
            A[i,j] = C * np.abs(np.sin(lat1) - np.sin(lat2))*np.abs(lon1 - lon2)
    return A


# REDO THIS FUNCTION USING PANDAS DATAFRAMES?
def area_weighted_mean(X, A):
    """Compute the area-weighted-average time series from a 3d array.

    Parameters
    ----------
    X : 3d-array 
        Array containing one time series per grid-cell, where the
        first dimension (i) is the time, and the second and thrid 
        dimensions (j and k) are the spatial coordinates (x,y).
    A : 2d-array
        A grid containing the area of each grid-cell on X, i.e.,
        the spatial coordinates.

    Returns
    -------
    out : 1d-array
        The average time series weighted by each grid-cell area.

    See also
    --------
    get_area_cells

    """
    nt, _, _ = X.shape
    ts = np.zeros(nt, 'f8')  # container for output time series
    X[X==0] = np.nan
    X = np.ma.masked_invalid(X)
    for k in range(nt):      # weight-average each 2d time step
        G = X[k,...]
        W = A.copy()
        W[np.isnan(G)] = 0
        s = W.sum()
        if s != 0:
            W /= s 
        else:
            W[:] = 0 
        ts[k] = np.sum(W*G)  # area-weighted average
    return ts


