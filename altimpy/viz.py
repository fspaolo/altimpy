"""
Module with some utility functions for visualization/plotting.

Notes
-----
Some of the functionalities were taken/modified from:
http://earth.usc.edu/~gely/coseis/www/index.html
https://github.com/gely/coseis

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# December 15, 2011

__author__ = 'Fernando Paolo'

import os
import subprocess
import cStringIO
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from const import *

### Visualization utilities

def make_cmap(colors, position=None, name='newcmap', n=256):
    """Creates a custom colormap for Matplotlib.

    Parameters
    ----------
    colors : list of tuples
        Contains RGB values (r, g, b). The RGB values may either be in 8-bit
        [0 to 255] or arithmetic [0 to 1]. Arrange your tuples so that the
        first color is the lowest value for the colorbar and the last is the
        highest.
    position : list
        Contains values from 0 to 1 to dictate the location of each color.
        Useful for irregular color segments.
    name : string
        The name of the colormap.
    n : int
        The number of intervals in the colormap (e.g. discrete colorbar).

    Return
    ------
    cmap : a colormap with equally spaced colors.

    Example
    -------
    # uneven colormap with 20 intervals
    colors = [(0, 0, 1), (.5, .5, 1), (1, 1, 1)]
    position = [0, 0.3, 1]
    cmap = make_cmap(colors, position, name='BlueWhite', n=20)

    Credits
    -------
    Chris Slocum (initial version)
    Fernando Paolo (several additions/modifications)

    TODO
    -----
    Add support for alpha in the tuples (forth channel) => (r, g, b, a)

    """
    if position == None:
        position = np.linspace(0, 1, len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    # if color in range [0,255]
    if np.max(colors) > 1:
        bit_rgb = np.linspace(0, 1, 256)
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]],
                         )
    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict, n)
    plt.register_cmap(cmap=newcmap)
    return newcmap


def shift_cmap(cmap, start=0, midpoint=0.5, stop=1, name='shiftedcmap'):
    '''Offset the median value of a colormap.
    
    And scale the remaining color range. Useful for data with a negative
    minimum and positive maximum where you want the middle of the colormap's
    dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave 
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax) 
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin)) 
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave 
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin)) 

    Credits
    -------
    Paul H (initial version)
    Horea Christian (additions/modifications)
    Fernando Paolo (additions/modifications)

    TODO
    ----
    Set 'start' and 'stop' dynamically when negative/positive bounds.

    '''
    # if array given, find optimal value to center new cmap
    if np.ndim(midpoint) != 0:
        midpoint = np.asarray(midpoint)[~np.isnan(midpoint)]
        midpoint = abs(midpoint.min()) / float(abs(midpoint.max()) + \
                                               abs(midpoint.min())) 
    # regular index to compute the colors
    reg_index = np.hstack([
        np.linspace(start, 0.5, 128, endpoint=False), 
        np.linspace(0.5, stop, 129, endpoint=True)
    ])
    # shifted index to match the midpoint of the data
    new_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
    for ri, si in zip(reg_index, new_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict, 256)
    plt.register_cmap(cmap=newcmap)
    return newcmap


def create_cmap(cmap, colorexp=1.0, nmod=0, modlim=0.5, upsample=True, 
              invert=False):
    """Color map creator.

    Parameters
    ----------
    cmap: Either a named colormap from cmap_lib or a 5 x N array,
        with rows specifying: (value, red, green, blue, alpha) components.
    colorexp: Exponent applied to the values to shift the colormap.
    nmod: Number of brightness modulations applied to the colormap.
    modlim: Magnitude of brightness modulations.
    upsample: Increase the number of samples if non-linear map (colorexp != 1)
    invert: Invert the order of colors.
    """
    if type(cmap) is str:
        cmap = cmap_lib[cmap]
    cmap = np.array(cmap, 'f')
    if invert:
        cmap[1:,:] = cmap[1:,::-1]  # only invert the color rows
    cmap[1:] /= max(1.0, cmap[1:].max())
    v, r, g, b, a = cmap
    v /= v[-1]
    if upsample and colorexp != 1.0:
        n = 16
        x  = np.linspace(0.0, 1.0, len(v))
        xi = np.linspace(0.0, 1.0, (len(v) - 1) * n + 1)
        r = np.interp(xi, x, r)
        g = np.interp(xi, x, g)
        b = np.interp(xi, x, b)
        a = np.interp(xi, x, a)
        v = np.interp(xi, x, v)
        v = np.sign(v) * abs(v) ** colorexp
    v = (v - v[0]) / (v[-1] - v[0])
    if nmod > 0:
        if len(v) < 6 * nmod:
            vi = np.linspace(v[0], v[-1], 8 * nmod + 1)
            r = np.interp(vi, v, r)
            g = np.interp(vi, v, g)
            b = np.interp(vi, v, b)
            a = np.interp(vi, v, a)
            v = vi
        w1 = np.cos(np.pi * 2.0 * nmod * v) * modlim
        w1 = 1.0 - np.maximum(w1, 0.0)
        w2 = 1.0 + np.minimum(w1, 0.0)
        r = (1.0 - w2 * (1.0 - w1 * r))
        g = (1.0 - w2 * (1.0 - w1 * g))
        b = (1.0 - w2 * (1.0 - w1 * b))
        a = (1.0 - w2 * (1.0 - w1 * a))
    return np.array([v, r, g, b, a])


def gmt_cmap(*args, **kwargs):
    """
    GMT style colormap (cpt). See `create_cmap` for details.
    """
    v, r, g, b, a = create_cmap(*args, **kwargs)
    cmap = ''
    fmt = '%-10r %3.0f %3.0f %3.0f     %-10r %3.0f %3.0f %3.0f\n'
    for i in range(len(v) - 1):
        cmap += fmt % (
            v[i],   255 * r[i],   255 * g[i],   255 * b[i],
            v[i+1], 255 * r[i+1], 255 * g[i+1], 255 * b[i+1],
        )
    return cmap


def pv_cmap(fname, cname, cmap):
    """Create ParaView's XML colormap file.

    Parameters
    ----------
    fname : string
        XML file name.
    cname : string
        Colorbar name.
    cmap : 5 x N array
        Colormap with 'rows' specifying (value, red, green, blue, alpha), and
        'columns' the discrete values between 0 and 1 in the colorscale.

    Example
    -------
    import altimpy as ap
    cmap = ap.create_cmap(ap.cmap_lib['wrgb'])
    ap.pv_cmap('wrgb.xml', 'wrgb', cmap)

    """
    v, r, g, b, a = cmap
    f = open(fname, 'w')
    f.write('<ColorMap name="%s" space="RGB">\n' % cname)
    for i, j, k, l, m in zip(v, a, r, g, b):
        f.write('<Point x="%f" o="%f" r="%f" g="%f" b="%f"/>\n' 
                % (i, j, k, l, m))
    f.write('</ColorMap>\n')
    f.close()


def mayavi_cmap(*args, **kwargs):
    """
    Mayavi colormap. See 'viz.cmap' for details.
    """
    cmap = create_cmap(*args, **kwargs)
    v, r, g, b, a = cmap
    if len(v) < 1001:
        vi = np.linspace(v[0], v[-1], 2001)
        r = np.interp(vi, v, r)
        g = np.interp(vi, v, g)
        b = np.interp(vi, v, b)
        a = np.interp(vi, v, a)
        cmap = np.array([r, g, b, a])
    return 255 * cmap.T

"""
5 x N cmap array is:

row 1 - the discretre position of each defined color
row 2 - the amount of red at each position (0 to 2)
row 3 - the amount of green at each position (0 to 2)
row 4 - the amount of blue at each position (0 to 2)
row 5 - the amount of transparency at each position (0 to 2)

Note 1: Each column represents discrete positions to be interpolated to create
color intervals (range of colors).

Note 2: Progressive colormaps (0,..,N) cannot be inverted, only divergent
colormaps (-N,..,N) can be inverted.
"""

cmap_lib = {
    'wwwwbgr': [
        (0, 4, 5, 7, 8, 9, 11, 12),
        (2, 2, 0, 0, 0, 2, 2, 2),
        (2, 2, 1, 2, 2, 2, 1, 0),
        (2, 2, 2, 2, 0, 0, 0, 0),
        (2, 2, 2, 2, 2, 2, 2, 2),
    ],
    'wwwbgr': [
        (0, 2, 3, 5, 6, 7, 9, 10),
        (2, 2, 0, 0, 0, 2, 2, 2),
        (2, 2, 1, 2, 2, 2, 1, 0),
        (2, 2, 2, 2, 0, 0, 0, 0),
        (2, 2, 2, 2, 2, 2, 2, 2),
    ],
    'wwbgr': [
        (0, 1, 2, 4, 5, 6, 8, 9),
        (2, 2, 0, 0, 0, 2, 2, 2),
        (2, 2, 1, 2, 2, 2, 1, 0),
        (2, 2, 2, 2, 0, 0, 0, 0),
        (2, 2, 2, 2, 2, 2, 2, 2),
    ],
    'wbgr': [
        (0,  1,  3,  4,  5,  7,  8),
        (2,  0,  0,  0,  2,  2,  2),
        (2,  1,  2,  2,  2,  1,  0),
        (2,  2,  2,  0,  0,  0,  0),
        (2,  2,  2,  2,  2,  2,  2),
    ],
    'wrgb': [
        (0,  1,  3,  4,  5,  7,  8),
        (2,  2,  2,  0,  0,  0,  0),
        (2,  1,  2,  2,  2,  1,  0),
        (2,  0,  0,  0,  2,  2,  2),
        (2,  2,  2,  2,  2,  2,  2),
    ],
    'bgr': [
        (0,  1,  3,  4,  5,  7,  8),
        (0,  0,  0,  0,  2,  2,  2),
        (0,  1,  2,  2,  2,  1,  0),
        (2,  2,  2,  0,  0,  0,  0),
        (2,  2,  2,  2,  2,  2,  2),
    ],
    'rgb': [
        (0,  1,  3,  4,  5,  7,  8),
        (2,  2,  2,  0,  0,  0,  0),
        (0,  1,  2,  2,  2,  1,  0),
        (0,  0,  0,  0,  2,  2,  2),
        (2,  2,  2,  2,  2,  2,  2),
    ],
    'r2b': [
        (-3, -1.5, -0.1,   0, 0.5,  1.5,  3),
        ( 2,    2,    2,   2,   0,    0,  0),
        ( 0,    1,    2,   2,   2,    1,  0),
        ( 0,    0,    0,   2,   2,    2,  2),
        ( 2,    2,    2,   2,   2,    2,  2),
    ],
    'bwr': [
        (-4, -3, -1,  0,  1,  3,  4),
        ( 0,  0,  0,  2,  2,  2,  1),
        ( 0,  0,  2,  2,  2,  0,  0),
        ( 1,  2,  2,  2,  0,  0,  0),
        ( 2,  2,  2,  2,  2,  2,  2),
    ],
    'cwy': [
        (-2, -1,  0,  1,  2),
        ( 0,  0,  1,  1,  1),
        ( 1,  0,  1,  0,  1),
        ( 1,  1,  1,  0,  0),
        ( 1,  1,  1,  1,  1),
    ],
    'w2r': [
        (0,  1,  3,  4),
        (2,  2,  2,  1),
        (2,  2,  0,  0),
        (2,  0,  0,  0),
        (2,  2,  2,  2),
    ],
    'w2b': [
        (0,  1,  3,  4),
        (2,  0,  0,  0),
        (2,  2,  1,  0),
        (2,  2,  2,  2),
        (2,  2,  2,  2),
    ],
}

### Map projection utilities

# DEPRECATED - do it directly
def align_data_with_fig(x, y, data):
    """Align map grid with figure frame.
    
    Map proj origin (x0,y0) is at the center of the grid 
    (polar stereo) and y-dim increases up, figure origin is at 
    llcorner (cartesian) and y-dim increases down. 
    """
    x = np.asarray(x)
    y = np.asarray(y)
    x = x - x.min()             # shift
    y = y - y.min()             # shift
    #y = y[::-1] - y.min()      # reverse y-dim and shift
    #data = data[::-1,:]         # reverse y-dim
    return [x, y, data]


def plot_img_subreg(m, x, y, img, bbox, **kwargs):
    """Plot image subregion defined by projection 'm'.
    
    m : Basemap projection defining subregion
    x, y : 1d arrays of coordinates
    img : 2d array with image
    bbox : low-left and upp-right lon/lat

    """
    bbox_sub = (m.llcrnrlon, m.llcrnrlat, m.urcrnrlon, m.urcrnrlat) 
    # shift coords and flip img to align projection and figure
    x = x - x.min()
    y = y - y.min()
    # img coords m -> km
    x /= 1e3
    y /= 1e3
    # img fig domain
    m1 = make_proj_stere(bbox)
    # subregion corners in fig coords 
    x0, y0 = m1(bbox_sub[0], bbox_sub[1])
    x1, y1 = m1(bbox_sub[2], bbox_sub[3])
    # select img subregion
    j, = np.where((x0 < x) & (x < x1))
    i, = np.where((y0 < y) & (y < y1))
    img2 = img[i[0]:i[-1], j[0]:j[-1]]
    x2 = x[j[0]:j[-1]]
    y2 = y[i[0]:i[-1]]
    # plot img
    m.imshow(img2, **kwargs)
    return [m, x2, y2, img2]


def make_proj_stere(bbox, ax=None, lat_ts=-71, lon_0=0, lat_0=-90):
    """Make a `basemap` polar stereographic projection.

    bbox : is (lon0, lat0, lon1, lat1), the ll and ur corners.
    lat_ts : is standard lat (true scale in the projection).
    lon_0, lat_0 : is the proj center (NOT grid center!).

    Default values are the MOA/LIMA parameters: 
    http://nsidc.org/data/moa/users_guide.html

    """
    # Ellipsoid: http://nsidc.org/data/polar_stereo/ps_grids.html
    a = RE
    b = a * np.sqrt(1.0 - E2) 
    m = Basemap(
        ax=ax,
        projection='stere', 
        lat_ts=lat_ts,         
        lon_0=lon_0, lat_0=lat_0,
        # corners in lon,lat
        llcrnrlon=bbox[0], llcrnrlat=bbox[1],
        urcrnrlon=bbox[2], urcrnrlat=bbox[3], 
        # ellipsoid for NSIDC grids:
        rsphere=(a, b),
        # ellipsoid WGS84:
        #rsphere=(6378137.00, 6356752.3142),  
        # sphere (mean radius):
        #rsphere=6371200.,                    
        )
    return m


def plot_grid_proj(m, lon, lat, grid, shift=True, masked=True, 
                   contourf=False, **kwargs):
    """Plot a rectangular grid on a given map projection.

    Parameters
    ----------
    m : Basemap object
        The projection.
    lon, lat : ndarray (1d or 2d)
        The grid coordinates.
    grid : 2d ndarray
        The grid to plot.
    shift : bool, optional
        Shift the grid by half cell-size due to `pcolormesh()`.
    masked : bool, optional
        Mask NaN values, i.e., plot only valid entries.
    contourf : bool, optional
        Plot interpolated field instead of grid.

    Returns
    -------
    out : Basemap object
        The provided map projection object.

    """
    if np.ndim(lon) == 1:
        lon, lat = np.meshgrid(lon, lat)
    if shift:
        lon -= (lon[1] - lon[0]) / 2.
        lat -= (lat[1] - lat[0]) / 2.
    # map lon/lat into proj and fig coords.
    xx, yy = m(lon, lat)
    if masked:
        grid = np.ma.masked_invalid(grid)
    if contourf:
        m.contourf(xx, yy, grid, **kwargs)
    else:
        m.pcolormesh(xx, yy, grid, **kwargs)
    return m


def get_gtif_subreg(m, filename, res=1):
    bbox_sub = (m.llcrnrlon, m.llcrnrlat, m.urcrnrlon, m.urcrnrlat) 
    # img coords (assuming polar stere)
    x, y, data, bbox_ll = get_gtif(filename, lat_ts=-71)          
    # downsample
    x, y, data = x[::res], y[::res], data[::res,::res]
    # img coords -> fig coords
    x, y, data = align_data_with_fig(x, y, data) 
    # fig domain
    m1 = make_proj_stere(bbox_ll)
    # subregion corners in fig coords 
    x0, y0 = m1(bbox_sub[0], bbox_sub[1])
    x1, y1 = m1(bbox_sub[2], bbox_sub[3])
    # extract mask subregion
    j, = np.where((x0 < x) & (x < x1))
    i, = np.where((y0 < y) & (y < y1))
    data_sub = data[i[0]:i[-1], j[0]:j[-1]]
    x_sub = x[j[0]:j[-1]]
    y_sub = y[i[0]:i[-1]]
    # plot mask
    #data_sub = np.ma.masked_values(data_sub, 0)
    #m.imshow(data_sub, cmap=plt.cm.gray)
    return [x_sub, y_sub, data_sub]


### Matplotlib utilities

def text(ax, x, y, s, edgecolor=None, edgealpha=0.1, edgewidth=0.75, 
         npmb=16, **kwargs):
    """Matplotlib text command augmented with poor man's bold.

    See plt.text() doctring for complete documentation.

    """
    h = [ax.text(x, y, s, **kwargs)]
    h[0].zorder += 1
    if edgecolor is not None:
        if 'bbox' in kwargs:
            del(kwargs['bbox'])
        kwargs['color'] = edgecolor
        kwargs['alpha'] = edgealpha
        aspect = ax.get_aspect()
        dx, dy = ax.get_position().size * ax.figure.get_size_inches() * 72.0
        x1, x2 = ax.get_xbound()
        y1, y2 = ax.get_ybound()
        dx = edgewidth * (x2 - x1) / dx
        dy = edgewidth * (y2 - y1) / dy
        if aspect == 'equal':
            dx = dy
        m = np.sqrt(0.5)
        dx = dx / m
        dy = dy / m
        for i in range(npmb):
            phi = 2.0 * np.pi * (i + 0.5) / npmb
            x_ = x + dx * np.cos(phi)
            y_ = y + dy * np.sin(phi)
            #x_ = x + dx * np.maximum(-m, np.minimum(m, np.cos(phi)))
            #y_ = y + dy * np.maximum(-m, np.minimum(m, np.sin(phi)))
            h += [ax.text(x_, y_, s, **kwargs)]
    return h


def get_cmap(*args, **kwargs):
    """Matplotlib enhanced colormap. 
    
    See 'create_cmap' for details.
    """
    n = kwargs.pop('n', 2001)
    v, r, g, b, a = create_cmap(*args, **kwargs)
    cdict = {'red': np.c_[v, r, r],
             'green': np.c_[v, g, g],
             'blue': np.c_[v, b, b] }
    cmap = mpl.colors.LinearSegmentedColormap('cmap', cdict, n)
    return cmap


def colorbar(fig, cmap, clim, title=None, rect=None, ticks=None, 
             ticklabels=None, boxcolor='k', boxalpha=1.0, boxwidth=0.2,
             shifttitle=1.8, shiftlabel=-0.2, **kwargs):
    """Matplotlib enhanced colorbar.

    **kwargs are passed to the text() command.
    
    Original by Geoffrey Ely.
    Modified by Fernando Paolo.

    """
    if rect is None:
        rect = 0.25, 0.04, 0.5, 0.02
    axis = clim[0], clim[1], 0, 1
    ax = fig.add_axes(rect)
    x = axis[0], axis[0], axis[1], axis[1], axis[0]
    y = axis[2], axis[3], axis[3], axis[2], axis[2]
    ax.plot(x, y, '-', c=boxcolor, lw=boxwidth*2, alpha=boxalpha, clip_on=False)
    ax.imshow([np.arange(1001)], cmap=cmap, extent=axis)
    ax.axis('off')
    ax.axis('tight')
    ax.axis(axis)
    title_voffset = shifttitle
    label_voffset = shiftlabel
    if title:
        x = 0.5 * (clim[0] + clim[1])
        text(ax, x, title_voffset, title, ha='center', va='baseline', **kwargs)
    if ticks is None:
        ticks = clim[0], 0.5 * (clim[0] + clim[1]), clim[1]
    if ticklabels is None:
        ticklabels = ticks
    for i, x in enumerate(ticks):
        s = '%s' % ticklabels[i]
        text(ax, x, label_voffset, s, ha='center', va='top', **kwargs)
    return ax


def length_scale(ax, x, y, w=None, label='%s', style='k-', linewidth=1, 
                color='k', **kwargs):
    """Draw a length scale bar between the points (x[0], y[0]) and (x[1], 
    y[1]).

    Original by Geoffrey Ely.
    Modified by Fernando Paolo.
    """
    x0 = 0.5 * (x[0] + x[1])
    y0 = 0.5 * (y[0] + y[1])
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    l = np.sqrt(dx*dx + dy*dy)
    if not w:
        x = ax.get_xlim()
        y = ax.get_ylim()
        x = abs(x[1] - x[0])
        y = abs(y[1] - y[0])
        if ax.get_aspect() == 'equal':
            w = 0.005 * (y + x)
        else:
            w = 0.01 / l * (y * abs(dx) + x * abs(dy))
    try:
        label = label % l
    except TypeError:
        pass
    rot = (dx, -dy), (dy, dx)
    x = -l, l, np.nan, -l, -l, np.nan,  l, l
    y =  0, 0, np.nan, -w,  w, np.nan, -w, w
    x, y = 0.5 / l * np.dot(rot, [x, y])
    theta = np.arctan2(dy, dx) * 180.0 / np.pi
    h1 = ax.plot(x0 + x, y0 + y, style, clip_on=False, lw=linewidth, c=color)
    h2 = text(ax, x0, y0, label, ha='center', va='center', rotation=theta, 
              color=color, **kwargs)
    return h1, h2


def compass_rose(ax, x, y, r, style='k-', **kwargs):
    """
    Original by Geoffrey Ely.
    Modified by Fernando Paolo.
    """
    theta = 0.0
    if 'rotation' in kwargs:
        theta = kwargs['rotation']
    kwargs.update(rotation_mode='anchor')
    c  = np.cos(theta / 180.0 * np.pi)
    s  = np.sin(theta / 180.0 * np.pi)
    x_ = x + r * np.array([(c,  s), (-c, -s)])
    y_ = y + r * np.array([(s, -c), (-s,  c)])
    h  = [ax.plot(x_, y_, style, clip_on=False)]
    x_ = x + r * np.array([(c, -c), (s, -s)]) * 1.3
    y_ = y + r * np.array([(s, -s), (-c,  c)]) * 1.3
    h += [text(ax, x_[0,0], y_[0,0], 'E', ha='left', va='center', **kwargs),
          text(ax, x_[0,1], y_[0,1], 'W', ha='right', va='center', **kwargs),
          text(ax, x_[1,0], y_[1,0], 'S', ha='center', va='top', **kwargs),
          text(ax, x_[1,1], y_[1,1], 'N', ha='center', va='bottom', **kwargs),]
    return h


def savefig(fig, fh=None, format=None, distill=False, **kwargs):
    """
    Enhanced version of Matplotlib savefig command.

    Takes the same arguments as savefig.  Saves to disk if a filename is
    given. Otherwise return a StringIO file descriptor, or a numpy array.  PDF is
    distilled using Ghostscript to produce smaller files.
    """
    import cStringIO
    if isinstance(fh, basestring):
        if format is None:
            format = fh.split('.')[-1]
        fh = open(os.path.expanduser(fh), 'wb')
    else:
        if format is None:
            format = 'array'
    out = cStringIO.StringIO()
    if format == 'array':
        if 'dpi' not in kwargs:
            kwargs['dpi'] = fig.dpi
        dpi = kwargs['dpi']
        n = fig.get_size_inches()
        n = int(n[1] * dpi), int(n[0] * dpi), 4
        fig.savefig(out, format='raw', **kwargs)
        out = np.fromstring(out.getvalue(), 'u1').reshape(n)
    elif distill and format == 'pdf':
        fig.savefig(out, format='eps', **kwargs)
        out = distill_eps(out)
    else:
        fig.savefig(out, format=format, **kwargs)
        out.reset()
    if fh is None:
        return(out)
    else:
        with fh:
            fh.write(out.getvalue())
        return


def digitize2(img, xlim=(-1, 1), ylim=(-1, 1), color='r'):
    """
    Digitize points on an image and rectify to a rectangular coordinate system.
    """
    import coord
    fig = plt.gcf()
    fig.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img)
    ax.axis('tight')
    ax.axis('off')
    plt.draw()
    plt.show()
    ax.hold(True)
    xx, yy = [], []
    for j in 0, 1:
        for k in 0, 1:
            print('Left-click %r' % [xlim[j], ylim[k]])
            x, y = fig.ginput(1, -1)[0]
            xx += [x]
            yy += [y]
            ax.plot([x], [y], '+' + color)
            plt.draw()
    xx = xx[:2], xx[2:]
    yy = yy[:2], yy[2:]
    print("""
    Left-click, space: add point
    Right-click, delete: cancel last point
    Enter: new line segment
    Enter twice: finish
    """)
    x0 = 0.5 * (xlim[1] + xlim[0])
    y0 = 0.5 * (ylim[1] + ylim[0])
    dx = 0.5 * (xlim[1] - xlim[0])
    dy = 0.5 * (ylim[1] - ylim[0])
    xr, yr = [], []
    while 1:
        xy = fig.ginput(-1, -1)
        if len(xy) == 0:
            break
        x, y = zip(*xy)
        ax.plot(x, y, '+-'+color)
        plt.draw()
        x, y = coord.ibilinear(xx, yy, x, y)
        x = x0 + dx * x
        y = y0 + dy * y
        print(x)
        print(y)
        xr += [x]
        yr += [y]
    return xr, yr


def contour(*args, **kwargs):
    """
    Extract contour polygons using matplotlib.
    """
    concat = True
    pp = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if concat:
        for cc in ax.contour(*args, **kwargs).collections:
            p = []
            for c in cc.get_paths():
                p += c.to_polygons() + [[[np.nan, np.nan]]]
            if p:
                del p[-1]
                pp += [np.concatenate(p).T]
            else:
                pp += [None]
    else:
        for cc in ax.contour(*args, **kwargs).collections:
            p = []
            for c in cc.get_paths():
                p += c.to_polygons()
            pp += [p]
    plt.close(fig)
    return pp


def intitle(title='', loc=1, size=None, pad=0., borderpad=0.5, borderwidth=4,
            borderalpha=1, ax=None, **kwargs):
    """Add title inside the figure. Same locations as 'label'.

    Examples
    --------
    # with pyplot single
    plt.figure()
    intitle('inner title', 3)

    # with pyplot subplot
    ax = plt.subplot(211)
    intitle('inner title', 3, ax=ax)

    # with pandas
    ax = df.plot()
    intitle('inner title', 3, ax=ax)

    """
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['axes.labelsize'])
    if ax is None:
        ax = plt.subplot(111)
    at = AnchoredText(title, loc=loc, prop=size, pad=pad, 
                      borderpad=borderpad, frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w",
                                              linewidth=borderwidth,
                                              alpha=borderalpha)])
    return ax


"""
Draws Hinton diagrams using matplotlib ( http://matplotlib.sf.net/ ).
Hinton diagrams are a handy way of visualizing weight matrices, using
colour to denote sign and area to denote magnitude.
 
By David Warde-Farley -- user AT cs dot toronto dot edu (user = dwf)
  with thanks to Geoffrey Hinton for providing the MATLAB code off of 
  which this is modeled.
 
Redistributable under the terms of the 3-clause BSD license 
(see http://www.opensource.org/licenses/bsd-license.php for details)
"""
 
def _blob(x, y, area, colour):
    """
    Draws a square-shaped blob with the given area (< 1) at
    the given coordinates.

    """
    hs = np.sqrt(area) / 2
    xcorners = np.array([x - hs, x + hs, x + hs, x - hs])
    ycorners = np.array([y - hs, y - hs, y + hs, y + hs])
    plt.fill(xcorners, ycorners, colour, edgecolor=colour)
 
def hinton(m, maxweight=None):
    """Draws a Hinton diagram for visualizing a weight matrix. 

    Temporarily disables matplotlib interactive mode if it is on, 
    otherwise this takes forever.

    Note
    ----
    Because the values of the matrix are assumed to be *weights*, NaN values
    are replaced by zero (otherwise the plot doesn't come good).

    """
    W = m.copy()
    W[np.isnan(W)] = 0

    reenable = False
    if plt.isinteractive():
        plt.ioff()
    
    plt.clf()
    height, width = W.shape
    if not maxweight:
        maxweight = 2**np.ceil(np.log(np.max(np.abs(W)))/np.log(2))
        
    plt.fill(np.array([0, width, width, 0]), 
             np.array([0, 0, height, height]),
             'gray')
    
    plt.axis('off')
    plt.axis('equal')
    for x in xrange(width):
        for y in xrange(height):
            _x = x+1
            _y = y+1
            w = W[y, x]
            if w > 0:
                _blob(_x - 0.5,
                      height - _y + 0.5,
                      min(1, w/maxweight),
                      'white')
            elif w < 0:
                _blob(_x - 0.5,
                      height - _y + 0.5, 
                      min(1, -w/maxweight), 
                      'black')
    if reenable:
        plt.ion()

# tools for shaded relief

def shade(data, cmap=plt.cm.gray, azimute=165.0, altitude=45.0, scale=1.0,
           minsat=1, maxsat=0, minval=0, maxval=1):
    """Creates shaded relief (hard light approach).

    Convenient wrap around matplotlib.colors.LightSource to facilitate
    scaling and saturation settings.

    Parameters
    ----------
    scale : float, default 1.0
        Scaling factor for shading. Larger number => larger gradients.
    minsat, maxsat, minval, maxval: int, default 1, 0, 0, 1
        Min and max saturation/values (from 0 to 1).

    See docstring of function 'shade2'.

    """
    # 1) create light source object
    ls = mpl.colors.LightSource(azimute=azimute, altitude=altitude)
    ls.hsv_min_sat = minsat
    ls.hsv_max_sat = maxsat
    ls.hsv_min_val = minval
    ls.hsv_max_val = maxval
    # 2) convert data to rgb array including shading from light source
    rgb = ls.shade(data * scale, cmap)
    return rgb


def shade2(data, intensity=None, cmap=plt.cm.gray, azimute=165.0, \
           altitude=45.0, scale=1.0):
    """Creates shaded relief (soft light approach).
    
    Sets shading for data array based on intensity layer or the data's value
    itself.

    Parameters
    ----------
    data : 2d array or masked array
        The grid to shade.
    intensity : 2d array of same shape as data
        The intensity layer for shading. If None, the data itself is used after
        getting the hillshade values (see hillshade).
    cmap : colormap, default plt.cm.gray
        e.g. matplotlib.colors.LinearSegmentedColormap instance.
    azimute, altitude, scale : floats, default 10.0, 165.0, 45.0
        Parameters for hillshade function (see hillshade).

    Output
    ------
    rgb : rgb set of the 'Pegtop soft light composition' (see Notes) of the
        data and intensity. It should be plotted with imshow().

    Credits
    -------
    Ran Novitsky Nof (initial version)
    Fernando Paolo (additions/modifications)

    Based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight

    See also
    --------
    hillshade
    shade

    Notes
    -----
    Matplotlib module enables hillshade method (shade) using a LightSource
    class (v 0.99). The problem is it uses the data itself as intensity and
    data. It is very useful for viewing a DEM but sometimes you would like the
    DEM as intensity underlying some other data. Another problem is that the
    shade method is producing a very light colored image sometimes even white
    where intensity is high.

    The difference in the shading colors is derived from the method used to
    produce it. While the Matplotlib method uses "hard light" method this
    function uses a "soft light" method. Matplotlib is converting the RGB
    colors to HSV and then calculating the new saturation and value according
    to the intensity. This function uses a formula based on the description of
    ImageMagick's 'pegtop_light', which is much faster as it is a single
    formula. Another advantage is the option to use separate layers for
    intensity and colors.

    Soft light: refers to light that tends to "wrap" around objects, casting
    diffuse shadows with soft edges. Soft light is when a light source is large
    relative to the subject.

    Hard light: the shadows produced will have 'harder' edges with less
    transition between illumination and shadow.

    Functions modified from:
    http://rnovitsky.blogspot.com/2010/04/using-hillshade-image-as-intensity.html

    """
    if intensity is None:
        # hilshading the data
        intensity = hillshade(data, azimute=azimute, altitude=altitude, scale=scale)
    else:
        # or normalize the intensity
        intensity = (intensity - intensity.min()) / \
                    (intensity.max() - intensity.min())
    # get rgb of normalized data based on cmap
    rgb = cmap((data - data.min()) / float(data.max() - data.min()))[:,:,:3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2 * d * rgb + (rgb**2) * (1 - 2 * d)
    return rgb


def hillshade(data, azimute=165.0, altitude=45.0, scale=1.0):
    """Convert data to hillshade (intensity).

    Parameters
    ----------
    data : 2d array
        The grid to be used as shading.
    azimute : float, default 165.0
        Direction where the light comes from: 0=south, 90=east, 180=north,
        270=west.
    altitude : float, default 45.0
        Altitude where the light comes from: 0=horison, 90=zenith.
    scale : float, default 1.0
        Z scaling factor for shading. Larger number => larger gradients.

    Output
    ------
    intensity : 2d array
        Normalized hilshade.

    See also
    --------
    shade2

    """
    # deg -> radians
    az = np.deg2rad(azimute)
    alt = np.deg2rad(altitude)
    # gradient in x and y directions
    dx, dy = np.gradient(data * scale)
    slope = 0.5 * np.pi - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dx, dy)
    intensity = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * \
                np.cos(-az - aspect - 0.5 * np.pi)
    intensity = (intensity - intensity.min()) / \
                (intensity.max() - intensity.min())
    return intensity

#---------------------------------------------------------------

def adjust_spines(ax, spines, pad=10):
    """Adjust the axis of a plot.

    Example
    -------
    # display only 'y-' and 'x-axis'
    adjust_spines(ax, ['left', 'bottom'], 10)

    Notes
    -----
    Do not share axis in 'subplot' environment!

    From
    ----
    http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html

    """
    for loc, spine in ax.spines.items():
        if loc in spines:
            # outward by 10 points
            spine.set_position(('outward', pad))
            spine.set_smart_bounds(True)
        else:
            # don't draw spine
            spine.set_color('none')             
    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def get_limits(x, decimals=1):
    """Get the ceiling and floor limits according number of decimals.
    
    Useful to set the y-/x-limits in matplotlib.

    """
    mn0 = x.min()
    mx0 = x.max()
    mn = mn0.round(decimals)
    mx = mx0.round(decimals)
    if mn > mn0:
        mn -= 1. / 10**decimals
    if mx < mx0:
        mx += 1. / 10**decimals
    return mn, mx


def meridians(start=0, stop=360, y0=-90, yn=90, spacing=10, npts=1000):
    """Draw meridians (arcs) from y0 to yn spaced by 'spacing' from 'start' to
    'stop' longitudes.
    """
    x, y = [], []
    lats = np.linspace(y0, yn, npts)  # single meridian (an arc)
    if np.ndim(spacing) == 0:
        spacing = np.arange(start, stop+spacing, spacing)
    for lon in spacing:
        lons = np.repeat([lon], lats.shape[0])
        x = np.hstack((x, lons))
        y = np.hstack((y, lats))
    return x, y


def parallels(start=-90, stop=90, x0=0, xn=360, spacing=10, npts=1000):
    """Draw parallels (circles) from x0 to xn spaced by 'spacing' from 'start'
    to 'stop' latitudes.
    """
    x, y = [], []
    lons = np.linspace(x0, xn, npts)  # single parallel (a circle)
    if np.ndim(spacing) == 0:
        spacing = np.arange(start, stop+spacing, spacing)
    for lat in spacing:
        lats = np.repeat([lat], lons.shape[0])
        x = np.hstack((x, lons))
        y = np.hstack((y, lats))
    return x, y


def rcparams():
    """Set optimal figure layout parameters."""
    plt.rcParams['font.family'] = 'arial'
    plt.rcParams['font.size'] = 14 # 16
    plt.rcParams['axes.labelsize'] = 18 # 20
    plt.rcParams['legend.fontsize'] = 14 # 16
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = 0.2
    plt.rcParams['grid.color'] = '0.5'
    #plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['xtick.major.size'] = 0
    plt.rcParams['ytick.major.size'] = 0
    plt.rcParams['xtick.labelsize'] = 15 # 18
    plt.rcParams['ytick.labelsize'] = 15 # 18
    plt.rcParams['xtick.major.pad'] = 6
    plt.rcParams['ytick.major.pad'] = 6
