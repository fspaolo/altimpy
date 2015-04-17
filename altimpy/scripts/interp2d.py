#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Interpolate using Gaussian Process Regression (kriging).

Uses the GP pacakge from 'sklearn' to interpolate spatial data points (2d).

"""

import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

import altimpy as ap

MASK_FILE = '/Users/fpaolo/data/masks/scripps/scripps_antarctica_mask1km_v1.h5'

np.random.seed(1)

try:
    fname = sys.argv[1]
except:
    fname = ''


# get fields to interpolate
f = tb.openFile(fname, 'a')
dhdt = f.root.dh_mean_xcal_short_const_dhdt[:]
#dgdt = f.root.dg_mean_xcal_dgdt[:]
lon = f.root.lon[:]
lat = f.root.lat[:]

# increase resolution
#lon, lat, dhdt = ap.regrid(lon, lat, dhdt, inc_by=2)

# generate 2d mask
m = ap.Mask()
m.get_mask(MASK_FILE)
m.gen_mask(lon, lat, dhdt)
m.arr_z[m.arr_z!=4] = 0
m.arr_z[m.arr_z==4] = 1
mask = m.arr_z.copy()

ij = np.where(~np.isnan(dhdt))  
mask[ij] = 2  # for plotting: higlight bins to interpolate
#mask[ij] = 0   # for interpolation: mask-out all valid entries
print 'bins to interpolate (%):', \
      (float(len(mask[mask==1])) / len(m.arr_z[m.arr_z==1])) * 1e2

plt.imshow(mask, origin='lower', interpolation='nearest')
plt.show()
sys.exit()

# filter out outliers (reduce variability) prior interpolating. Cut-off 3 std 
# for each time step independently. NOT NEEDED!
'''
print 'range before filtering:', np.nanmin(dhdt), np.nanmax(dhdt)
dhdt = ap.filter_std(dhdt, n=3, per_field=True)
print 'range after filtering:', np.nanmin(dhdt), np.nanmax(dhdt)
'''

# interpolate each field
krig = ap.Kriging2d()
field, error = krig.interpolate(dhdt, mask, lon, lat)

plt.figure()
plt.imshow(field, origin='lower', interpolation='nearest')
plt.figure()
plt.imshow(error, origin='lower', interpolation='nearest')
plt.show()
sys.exit()

# save interpolated fields
f.createArray('/', 'dh_mean_xcal_short_const_dhdt_interp', field)
f.createArray('/', 'dh_mean_xcal_short_const_dhdt_interp_err', error)
f.close()
