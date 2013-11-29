#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Interpolation using Gaussian Process Regression (krigin).

Uses the GP pacakge from 'sklearn' to interpolate spatial data points (2d
fields).

"""

import sys
import numpy as np
from sklearn.gaussian_process import GaussianProcess

np.random.seed(1)


class Krigin2d(object):
    """Interpolate using Gaussian Process Regression (krigin).

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

    def krigin(self, X, y, X_pred):
        """Interpolate using Gaussian Process Regression (krigin).

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

    def interpolate(self, field, mask, lat2d, lon2d):
        """Interpolate missing points in 'field' according to 'mask'.
        """
        error = np.empty_like(field) * np.nan

        ij_samp = np.where((mask == 1) & (~np.isnan(field)))
        ij_pred = np.where((mask == 1) & (np.isnan(field)))
        X_samp = np.column_stack((lon2d[ij_samp], lat2d[ij_samp]))
        X_pred = np.column_stack((lon2d[ij_pred], lat2d[ij_pred]))
        y_samp = field[ij_samp]

        y_pred, sigma = self.krigin(X_samp, y_samp, X_pred)
        field[ij_pred] = y_pred
        error[ij_pred] = sigma

        return [field, error]
