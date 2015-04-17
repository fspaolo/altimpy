#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generate a 2d grid with (linear) dh/dt per grid-cell from a 3d array.

"""

import os
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import altimpy as ap

# parameters
#---------------------------------------------------------------------

SAVE_TO_FILE = False

MOAFILE = '/Users/fpaolo/data/MOA/moa750_r1_hp1.tif'
MASKFILE = '/Users/fpaolo/data/masks/scripps/scripps_antarctica_mask1km_v1.h5'
MOA_RES = 10

MAX_DHDT = 0.5    # for trend
MAX_DH = 10       # for difference
MIN_NPTS = 6
FROM_YEAR, TO_YEAR = 1991, 2013
BBOX_REG = (-131, -60, 51, -58)       # Antarctica

width = 0.9  # gaussian width

cmap = ap.colormap('rgb')
legend = 'Elevation change rate 1992-2012'
ticklabels = '<%.2f m/yr' % -MAX_DHDT, '0', '>%.2f m/yr' % MAX_DHDT
ticks = -MAX_DHDT, 0, MAX_DHDT 
colorlim = ticks[0], ticks[-1]
colorexp = 1.0
inches = 6.4, 3.6
extent = (-121.1, -113.8), (32.1, 35.5)
ylim = 150.0
xlim = ylim * inches[0] / inches[1]
axis = -xlim, xlim, -ylim, ylim

#---------------------------------------------------------------------

fname = sys.argv[1]

def filter_npts(X, npts=6):
    _, ny, nx = X.shape
    for i in range(ny):
        for j in range(nx):
            ii, = np.where(~np.isnan(X[:,i,j]))
            if len(ii) < npts:
                X[:,i,j] = np.nan
    return X


def filter_abs(X, abs_val=2):
    X[np.abs(X)>abs_val] = np.nan
    return X

#---------------------------------------------------------------------

def plot_grid(m, lon, lat, dhdt, alpha=1, vmin=-MAX_DHDT, vmax=MAX_DHDT):
    print 'ploting grid...'
    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)  # 1d -> 2d
    xx, yy = m(lon, lat)
    dhdt = np.ma.masked_invalid(dhdt)
    m.pcolormesh(xx, yy, dhdt, cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha)
    print 'done.'
    return m


def plot_boundaries(m, region):
    l, r, b, t = region
    f = {}
    f[1] = tb.openFile('/Users/fpaolo/data/coastline/antarctica_gl_ll.h5')
    f[2] = tb.openFile('/Users/fpaolo/data/coastline/moa_coastfile_ll.h5')
    f[3] = tb.openFile('/Users/fpaolo/data/coastline/moa_islands_ll.h5')
    for ff in f.values():
        d = ff.root.data[::10]
        lat, lon = d[:,0], d[:,1]
        ii, = np.where((lon >= l) & (lat <= r) & (lat >= b) & (lat <= t))
        x, y = m(lon[ii], lat[ii])
        m.scatter(x, y, s=0.1, c='w', facecolor='0.5', lw = 0)
    return m

#---------------------------------------------------------------------

f = tb.openFile(fname, 'a')
#dh = f.root.dg_mean_xcal[:]
#dh = f.root.dh_mean_xcal[:]
dh = f.root.dh_mean_xcal_short_const[:]
time = ap.num2year(f.root.time_xcal[:])
lon = f.root.lon[:]
lat = f.root.lat[:]
xed = f.root.x_edges[:]
yed = f.root.y_edges[:]

# filter data
dh = filter_abs(dh, abs_val=MAX_DH)
dh = filter_npts(dh, npts=MIN_NPTS)

# get 2D grid with linear dh/dt
dhdt = ap.get_dydx(time, dh)

# gaussian smoothing
#dhdt = ap.gfilter(dhdt, width)

# filter values
#dhdt[np.abs(dhdt)>MAX_DHDT] = np.nan

# save data
if SAVE_TO_FILE:
    f.createArray('/', 'dh_mean_xcal_short_const_dhdt', dhdt)

#---------------------------------------------------------------------

# plot
fig = plt.figure(figsize=(10,9), frameon=False)
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('off')

# the MOA image 
m = ap.make_proj_stere(BBOX_REG)
x, y, moa, bbox_moa = ap.get_gtif(MOAFILE, units='m')
x, y, moa = x[::10], y[::10], moa[::10,::10]  # downsample

# mask-out ocean and islands
mask = ap.Mask()
mask.get_mask(MASKFILE)
mask.gen_mask(x, y, moa, coord='x/y')
ij = np.where((mask.arr_z == 0) | (mask.arr_z == 3))
moa[ij] = 0

m, _, _, _ = ap.plot_moa_subreg(m, x, y, moa, bbox_moa, 1)

m = plot_grid(m, lon, lat, dhdt, alpha=1, vmin=-.2, vmax=.2)
#m = plot_boundaries(m, (-180, 180, -90, -50))

# colorbar and text
rect = 0.17, 0.05, 0.24, 0.02
ap.colorbar(fig, cmap, colorlim, legend, rect, ticks, ticklabels, 
            size=14, weight='bold', color='k', edgecolor='w')
ap.rcparams()
#ap.intitle('Filchner-Ronne Ice Shelf', loc=1)
#ap.intitle('dh/dt 1992-2012', loc=3)
leg = fig.add_axes([0, 0, 1, 1])
leg.set_axis_off()

'''
plt.savefig('map20yr.pdf', dpi=150, bbox_inches='tight', pad_inches=-0.5)
os.system('cp map_fris2.pdf /Users/fpaolo/posters/agu12/figures/')
'''

plt.show()

ap.close_files()
