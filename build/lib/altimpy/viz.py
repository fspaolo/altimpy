"""
Module with some utility functions for visualization/plotting.

Notes
-----
Some of the functionalities were taken/modified from:
http://earth.usc.edu/~gely/coseis/www/index.html

"""
# Fernando Paolo <fpaolo@ucsd.edu>
# December 15, 2011

import os
import numpy as np
import subprocess, cStringIO
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

### Visualization utilities

def create_colormap(cmap, colorexp=1.0, nmod=0, modlim=0.5, upsample=True, 
              invert=False):
    """
    Color map creator.

    Parameters
    ----------
    cmap: Either a named colormap from colormap_lib or a 5 x N array,
        with rows specifying: (value, red, green, blue, alpha) components.
    colorexp: Exponent applied to the values to shift the colormap.
    nmod: Number of brightness modulations applied to the colormap.
    modlim: Magnitude of brightness modulations.
    upsample: Increase the number of samples if non-linear map (colorexp != 1)
    invert: Intert the order of colors.
    """
    if type(cmap) is str:
        cmap = colormap_lib[cmap]
    cmap = np.array(cmap, 'f')
    if invert:
        cmap = cmap[:,::-1]
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


def cpt(*args, **kwargs):
    """
    GMT style colormap. See `create_colormap` for details.
    """
    v, r, g, b, a = create_colormap(*args, **kwargs)
    cmap = ''
    fmt = '%-10r %3.0f %3.0f %3.0f     %-10r %3.0f %3.0f %3.0f\n'
    for i in range(len(v) - 1):
        cmap += fmt % (
            v[i],   255 * r[i],   255 * g[i],   255 * b[i],
            v[i+1], 255 * r[i+1], 255 * g[i+1], 255 * b[i+1],
        )
    return cmap


colormap_lib = {
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
}


### Map projection utilities

def get_gtif(fname, lat_ts=-71, lon_0=0, lat_0=-90):
    """
    Reads a GeoTIFF image and returns the respective 2D array.
    ==>it assumes polar stereographic proj<==

    Return
    ------
    x, y : 1D arrays containing the coordinates.
    img : 2D array containing the image.
    bbox_ll : lon/lat limits (lllon,lllat,urlon,urlat).

    Notes
    -----
    MOA parameters: http://nsidc.org/data/moa/users_guide.html
    lat_ts is standard lat (lat of true scale), or "latitude_of_origin"
    lon_0/lat_0 is proj center (NOT grid center!)

    To get GeoTIFF metadata:

        $ gdalinfo file.tif
    
    """
    try:
        from osgeo import osr, gdal
        import pyproj as pj
    except:
        msg = """some of the following modules are missing: 
        `osgeo` (GDAL) or `pyproj`"""
        raise ImportError(msg)
    ds = gdal.Open(fname) 
    img = ds.ReadAsArray()       # img -> matrix
    gt = ds.GetGeoTransform()    # coordinates
    nx = ds.RasterXSize          # number of pixels in x
    ny = ds.RasterYSize          # number of pixels in y 
    dx = gt[1]                   # pixel size in x
    dy = gt[5]                   # pixel size in y 
    xmin = gt[0]
    ymax = gt[3]
    # from: http://gdal.org/gdal_datamodel.html
    ymin = ymax + nx*gt[4] + ny*dy 
    xmax = xmin + nx*dx + ny*gt[2] 
    # Polar stereo coords x,y
    x = np.arange(xmin, xmax, dx)    
    # in reverse order -> raster origin = urcrn
    y = np.arange(ymax, ymin, dy)  
    # bbox of raster img in x,y 
    bbox_xy = (xmin, ymin, xmax, ymax)
    # bbox of raster img in lon,lat (to plot proj)
    p1 = pj.Proj(proj='stere', lat_ts=lat_ts, lon_0=lon_0, lat_0=lat_0)
    xmin, ymin = p1(bbox_xy[0], bbox_xy[1], inverse=True)
    xmax, ymax = p1(bbox_xy[2], bbox_xy[3], inverse=True)
    bbox_ll = (xmin, ymin, xmax, ymax)
    print 'image limits (left/right/bottom/top):'
    print '(x,y)', bbox_xy[0], bbox_xy[2], bbox_xy[1], bbox_xy[3]
    print '(lon,lat)', bbox_ll[0], bbox_ll[2], bbox_ll[1], bbox_ll[3]
    return [x, y, img, bbox_ll]


def align_data_with_fig(x, y, data, res=10):
    """
    Align map grid with figure frame. 
    
    Map proj origin (x0,y0) is at the center of the grid 
    (polar stereo) and y-dim increases up, figure origin is at 
    llcorner (cartesian) and y-dim increases down. 

    """
    x = np.asarray(x)
    y = np.asarray(y)
    x = x[::res] - x.min()        # shift
    y = y[::-res] - y.min()       # reverse y-dim and shift
    data = data[::-res, ::res]    # reverse y-dim
    return [x, y, data]


def plot_moa_subregion(m, moafile, res=10):
    bbox_sub = (m.llcrnrlon, m.llcrnrlat, m.urcrnrlon, m.urcrnrlat) 
    # real proj coords (stere)
    x, y, data, bbox_moa = get_gtif(moafile, lat_ts=-71)          
    # fig coords
    x, y, data = align_data_with_fig(x, y, data, res) 
    # MOA fig domain
    m1 = make_proj_stere(bbox_moa)
    # subregion corners in fig coords 
    x0, y0 = m1(bbox_sub[0], bbox_sub[1])
    x1, y1 = m1(bbox_sub[2], bbox_sub[3])
    # select MOA subregion
    j, = np.where((x0 < x) & (x < x1))
    i, = np.where((y0 < y) & (y < y1))
    data2 = data[i[0]:i[-1], j[0]:j[-1]]
    x2 = x[j[0]:j[-1]]
    y2 = y[i[0]:i[-1]]
    # plot MOA img
    data2 = np.ma.masked_values(data2, 0)
    m.imshow(data2, cmap=cm.gray)
    return [m, x2, y2, data2]


def make_proj_stere(bbox, lat_ts=-71, lon_0=0, lat_0=-90):
    """
    Make a `basemap` polar stereographic projection.

    bbox : is (lon0, lat0, lon1, lat1), the ll and ur corners.
    lat_ts : is standard lat (true scale in the projection).
    lon_0, lat_0 : is the proj center (NOT grid center!).

    Default values are the MOA parameters: 
    http://nsidc.org/data/moa/users_guide.html
    """
    # Ellipsoid: http://nsidc.org/data/polar_stereo/ps_grids.html
    a = RE
    b = a*np.sqrt(1.0 - ec2) 
    m = Basemap(
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


def plot_grid_proj(m, lon, lat, grid, cell=None, shift=True, contourf=False, **kw):
    """
    Plot a 2D rectangular grid on a given map projection `m`.

    m : Basemap projection
    lon, lat : 1D arrays of grid coordinates
    grid : 2D array (the grid)
    """
    if shift:
        # shift the grid due to `pcolormesh()`
        lon -= (lon[1] - lon[0])/2.
        lat -= (lat[1] - lat[0])/2.
    lon, lat = np.meshgrid(lon, lat)
    xx, yy = m(lon, lat)
    grid = np.ma.masked_invalid(grid)
    if contourf:
        m.contourf(xx, yy, grid, 25, **kw)
    else:
        m.pcolormesh(xx, yy, grid, **kw)
    if cell is not None:
        lon, lat = util.box(cell)
        x, y = m(lon, lat)
        m.plot(x, y, 'k', linewidth=2)
    return m


def get_gtif_subreg(m, filename, res=10):
    bbox_sub = (m.llcrnrlon, m.llcrnrlat, m.urcrnrlon, m.urcrnrlat) 
    # img coords (assuming polar stere)
    x, y, data, bbox_ll = get_gtif(filename, lat_ts=-71)          
    # img coords -> fig coords
    x, y, data = align_data_with_fig(x, y, data, res) 
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
    #m.imshow(data_sub, cmap=cm.gray)
    return [x_sub, y_sub, data_sub]


### Matplotlib utilities

def text(ax, x, y, s, edgecolor=None, edgealpha=0.1, edgewidth=0.75, 
         npmb=16, **kwargs):
    """
    Matplotlib text command augmented with poor man's bold.
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


def colormap(*args, **kwargs):
    """
    Matplotlib enhanced colormap. See `create_colormap` for details.
    """
    from matplotlib.colors import LinearSegmentedColormap
    v, r, g, b, a = create_colormap(*args, **kwargs)
    n = 2001
    cmap = { 'red':np.c_[v, r, r],
           'green':np.c_[v, g, g],
            'blue':np.c_[v, b, b] }
    cmap = LinearSegmentedColormap('cmap', cmap, n)
    return cmap


def colorbar(fig, cmap, clim, title=None, rect=None, ticks=None, 
             ticklabels=None, boxcolor='k', boxalpha=1.0, 
             boxwidth=0.2, **kwargs):
    """
    Matplotlib enhanced colorbar.
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
    if title:
        x = 0.5 * (clim[0] + clim[1])
        text(ax, x, 2, title, ha='center', va='baseline', **kwargs)
    if ticks is None:
        ticks = clim[0], 0.5 * (clim[0] + clim[1]), clim[1]
    if ticklabels is None:
        ticklabels = ticks
    for i, x in enumerate(ticks):
        s = '%s' % ticklabels[i]
        text(ax, x, -0.6, s, ha='center', va='top', **kwargs)
    return ax


def lengthscale(ax, x, y, w=None, label='%s', style='k-', **kwargs):
    """
    Draw a length scale bar between the points (x[0], y[0]) and (x[1], y[1]).
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
    h1 = ax.plot(x0 + x, y0 + y, style, clip_on=False)
    h2 = text(ax, x0, y0, label, ha='center', va='center', rotation=theta, **kwargs)
    return h1, h2


def compass_rose(ax, x, y, r, style='k-', **kwargs):
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
    h += [
        text(ax, x_[0,0], y_[0,0], 'E', ha='left', va='center', **kwargs),
        text(ax, x_[0,1], y_[0,1], 'W', ha='right', va='center', **kwargs),
        text(ax, x_[1,0], y_[1,0], 'S', ha='center', va='top', **kwargs),
        text(ax, x_[1,1], y_[1,1], 'N', ha='center', va='bottom', **kwargs),
    ]
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


def digitize(img, xlim=(-1, 1), ylim=(-1, 1), color='r'):
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


def add_inner_title(ax, title='', loc=1, size=None, **kwargs):
    """
    Add title inside the figure. Same locations as `label`.

    Example
    -------
    fig = plt.figure()
    ax = fig.add_subplot((111))
    ax = add_inner_title(ax, 'title', 3)
    """
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size, pad=0., 
        borderpad=0.5, frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=4)])
    at.patch.set_alpha(0.5)
    return at


def hinton(W, max_weight=None, ax=None):
    """
    Draws a Hinton diagram for visualizing a weight matrix.
    """
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import NullLocator
    W = W.T
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
    if not max_weight:
        max_weight = 2**np.ceil(np.log(np.abs(W).max())/np.log(2))
    ax.patch.set_facecolor('gray')
    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_locator(NullLocator())
    for (x,y),w in np.ndenumerate(W):
        if w > 0: color = 'white'
        else:     color = 'black'
        size = np.sqrt(np.abs(w))
        rect = Rectangle([x - size / 2, y - size / 2], size, size,
            facecolor=color, edgecolor=color)
        ax.add_patch(rect)
    ax.autoscale_view()
    # Reverse the yaxis limits
    ax.set_ylim(*ax.get_ylim()[::-1])


def plot_matrix(mat, title='', loc=1, plot=None, **kw):
    """
    Plot the representation of a matrix: matshow, hinton or spy.

    plot : can be 'matshow', 'hinton' or 'spy'. If `None` (default) 
    plots matshow and hinton diagrams.
    """
    from matplotlib.ticker import NullLocator
    if plot is None:
        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.matshow(mat)
        hinton(mat, ax=ax2)
        t = add_inner_title(ax1, title, loc=loc)
        t = add_inner_title(ax2, title, loc=loc)
    else:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot((111))
        if plot == 'matshow':
            ax.matshow(mat)
            ax.xaxis.set_major_locator(NullLocator())
            ax.yaxis.set_major_locator(NullLocator())
        elif plot == 'hinton':
            hinton(mat, ax=ax)
        elif plot == 'spy':
            ax.spy(mat, precision=0.1, markersize=6, **kw)
        else:
            raise ValueError('wrong argument `plot=%s`' % plot)
        t = add_inner_title(ax, title, loc=loc)
    return fig


