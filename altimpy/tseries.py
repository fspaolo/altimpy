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


def get_area_cells(X, x_edges, y_edges):
    """ Calculates the area of each grid-cell.

    Parameters
    ----------
    X : 2d-array
        The rectangular grid.
    x_edges, y_edges : 1d-array
        The coordinates of the nodes.

    Returns
    -------
    out : 2d array
        Same shape as X with the values of the area on each cell.

    """
    C = EARTH_RADIUS**2*D2R  # deg -> rad
    A = np.empty_like(X)
    ny, nx = X.shape
    for i in xrange(ny):
        for j in xrange(nx):
            lat1, lat2 = y_edges[i]*D2R, y_edges[i+1]*D2R
            lon1, lon2 = x_edges[j]*D2R, x_edges[j+1]*D2R
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


