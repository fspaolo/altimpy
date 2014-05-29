#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Some interpolation functionalities.

Notes
-----
The Gaussian Process Regression (kriging) uses the GP pacakge from 'sklearn'
to interpolate spatial data points (2d fields).

"""

import sys
import numpy as np
from sklearn.gaussian_process import GaussianProcess

np.random.seed(1)


class Kriging2d(object):
    """Interpolate using Gaussian Process Regression (kriging).

    Uses the GP pacakge from 'sklearn' to interpolate spatial data points (2d
    fields).
    """

    def __init__(self, regr='constant', corr='squared_exponential',
                 theta0=[10] * 2, thetaL=[1e-1] * 2, thetaU=[20] * 2,
                 nugget=1e5 * sys.float_info.epsilon):
        self.regr = regr
        self.corr = corr
        self.theta0 = theta0
        self.thetaL = thetaL
        self.thetaU = thetaU
        self.nugget = nugget
        self.rand_start = 50

    def kriging(self, X, y, X_pred):
        """Interpolate using Gaussian Process Regression (kriging).

        Uses the GP pacakge from 'sklearn' to interpolate spatial data points
        (2d).

        Interpolation equal noiseless case, i.e., "almost" no uncertainty in
        the observations.

        Bounds are defined assuming anisotropy.

        """
        # instanciate a Gaussian Process model
        gp = GaussianProcess(regr=self.regr, corr=self.corr,
                             theta0=self.theta0, thetaL=self.thetaL,
                             thetaU=self.thetaU, random_start=self.rand_start,
                             nugget=self.nugget, verbose=True)

        # fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, y)

        # evaluate the prediction points (ask for MSE as well)
        y_pred, MSE = gp.predict(X_pred, eval_MSE=True)

        return [y_pred, np.sqrt(MSE)]

    def interpolate(self, field, mask, lon, lat):
        """Interpolate missing cells in 'field' defined on 'mask' (values=1).
        """
        if np.ndim(lon) == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        error = np.empty_like(field) * np.nan

        ij_samp = np.where(~np.isnan(field))  # data to generate the model
        ij_pred = np.where(mask == 1)         # points to interpolate
        X_samp = np.column_stack((lon2d[ij_samp], lat2d[ij_samp]))
        X_pred = np.column_stack((lon2d[ij_pred], lat2d[ij_pred]))
        y_samp = field[ij_samp]

        y_pred, sigma = self.kriging(X_samp, y_samp, X_pred)
        field[ij_pred] = y_pred
        error[ij_pred] = sigma

        return [field, error]


def med_impute(arr, ij, size=3, min_pixels=3):
    """Median imputation for 2d array.

    Fill in i,j elements using median of footprint.

    It handleds NaNs.
    It uses a minimum number of non-null pixels.
    ij = ([i0,i1,..], [j0,j1,..])
    """
    def median(x, min_pixels=3):
        x = x.ravel()
        central_pixel = x[(len(x)-1)/2]
        valid_pixels = x[~np.isnan(x)]
        if len(valid_pixels) < min_pixels:
            m = central_pixel
        else:
            m = np.median(valid_pixels)
        return m
    l = size / 2
    for i, j in zip(ij[0], ij[1]):
        try:  # 'try' to avoid border issues. Add borders and update indices. FIXME
            arr[i,j] = median(arr[i-l:i+l+1, j-l:j+l+1], min_pixels=min_pixels)
        except:
            print 'footprint outside the boundaries, need to add borders!'
            pass
    return arr
