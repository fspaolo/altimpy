"""
Module with input/output functions.

"""

import os
import re
import numpy as np
import scipy as sp
import tables as tb
import pandas as pd
import datetime as dt
import netCDF4 as nc

try:
    from osgeo import osr, gdal
    import pyproj as pj
except:
    msg = """one of the following modules are missing: 
    `osgeo` (GDAL) or `pyproj`"""
    raise ImportError(msg)

from altimpy import is_meter


# definition of Table structures for HDF5 files

class TimeSeries(tb.IsDescription):
    sat_name = tb.StringCol(20, pos=1)
    ref_time = tb.StringCol(20, pos=2)
    date = tb.StringCol(20, pos=3)
    year = tb.Int32Col(pos=4)
    month = tb.Int32Col(pos=5)
    dh_mean = tb.Float64Col(pos=6)
    dh_error = tb.Float64Col(pos=7)
    dg_mean = tb.Float64Col(pos=8)
    dg_error = tb.Float64Col(pos=9)
    n_ad = tb.Int32Col(pos=10)
    n_da = tb.Int32Col(pos=11)


class TimeSeriesGrid(tb.IsDescription):
    sat_name = tb.StringCol(20, pos=1)
    ref_time = tb.StringCol(20, pos=2)
    date = tb.StringCol(20, pos=3)
    year = tb.Int32Col(pos=4)
    month = tb.Int32Col(pos=5)


def close_files():
    for fid in tb.file._open_files.values():
        fid.close() 


def add_cols_to_tbl(fname, tname, cols):
    """Add columns to an existing table."""
    # Open it again in append mode
    f = tb.openFile(fname, "a")
    table = f.getNode(tname) 
    # Get a description of table in dictionary format
    descr = table.description._v_colObjects
    descr2 = descr.copy()
    # Add a column to description
    for cname, cval in cols.items():
        descr2[cname] = tb.Col.from_dtype(cval.dtype)
    # Create a new table with the new description
    table2 = f.createTable('/', tname[1:]+'2', descr2, "temporary table", tb.Filters(9))
    # Copy the user attributes
    table.attrs._f_copy(table2)
    # Fill the rows of new table with default values
    for i in xrange(table.nrows):
        table2.row.append()
    # Flush the rows to disk
    table2.flush()
    # Copy the columns of source table to destination
    for col in descr:
        getattr(table2.cols, col)[:] = getattr(table.cols, col)[:]
    # Fill the new column(s)
    for cname, cval in cols.items():
        getattr(table2.cols, cname)[:] = cval[:] 
    # Remove the original table
    table.remove()
    # Move table2 to table
    table2.move('/', tname[1:])
    # Print the new table
    print "new table with added column(s):", f
    # Finally, close the file
    f.close()


def save_arr_as_tbl(fname, tname, cols):
    """
    Given 1D arrays save (or add if file exists) a Table.

    fname : name of new or existing file.
    tname : name of new table.
    cols : a dictionary {'colname': colval, ...}.
    """
    # Create column description
    descr = {}
    for i, (cname, cval) in enumerate(cols.items()):
        descr[cname] = tb.Col.from_dtype(cval.dtype, pos=i)
    f = tb.openFile(fname, 'a')  # if doesn't exist create it
    table = f.createTable('/', tname, descr, "", tb.Filters(9))
    table.append([v for k, v in cols.items()])
    table.flush()
    print "file with new table:", f


def save_arr_as_mat(fname, arrs, complib='blosc'):
    """
    Given 1D and/or 2D arrays save as a column matrix (2D array).

    fname : name of file to be saved.
    arrs : a list with 1D/2D arrays with *same first dimension*.
    """
    nrow, ncol = 0, 0
    for a in arrs:
        if a.ndim > 1:
            ncol += a.shape[1]
        else:
            ncol += 1
        nrow = a.shape[0]
    f = tb.openFile(fname, 'w')
    atom = tb.Atom.from_dtype(np.dtype('f8'))
    shape = (nrow, ncol)
    filters = tb.Filters(complib=complib, complevel=9)
    d = f.createCArray('/','data', atom=atom, shape=shape, filters=filters)
    j1, j2 = 0, 0
    for a in arrs:
        if a.ndim > 1:
            j2 += a.shape[1]
        else:
            a = a.reshape(nrow, 1)
            j2 += a.shape[1]
        d[:,j1:j2] = a
        j1 = j2
    print "file with new array:", f
    f.close()


### save to netcdf file


class NetCDF(object):
    """Quick and dirty way to save ndarrays to NetCDF4.
    
    Example
    -------
    >>> import numpy
    >>> arr = numpy.arange(100).reshape(10,10)
    >>> f = NetCDF('file.nc')
    >>> f.create_var('name_var', ('name_dim1', 'name_dim2'), arr)
    created dimension: name_dim1 (n=10)
    created dimension: name_dim2 (n=10)
    created variable: name_var ('name_dim1', 'name_dim2')
    >>> f.close()
    """
    def __init__(self, fname):
        self.filename = fname
        self.file = nc.Dataset(fname, 'w', format='NETCDF4')

    def create_var(self, varname, dimnames, arr):
        dimnames = np.atleast_1d(dimnames)
        shape = list(arr.shape)
        if len(shape) > 2:  
            # unlimited lenght for first dim (z or t)
            shape[0] = 0
        for n, dname in zip(shape, dimnames):
            # create dim if doesn't exit
            if not self.file.dimensions.has_key(dname):
                self.file.createDimension(dname, n)
                print 'created dimension: %s (n=%d)' % (dname, n)
        # create var if doesn't exist 
        if not self.file.variables.has_key(varname):
            var = self.file.createVariable(varname, 'f8', dimnames)
            var[:] = arr[:]
            print 'created variable:', varname, dimnames

    def close(self):
        self.file.close()


def read_gtif(fname, max_dim=1024.):
    """Load GeoTIFF image into array (RGB or grayscale).
    
    To get GeoTIFF metadata:

        $ gdalinfo file.tif
    
    """
    # get nx, ny, nz and geoinfo
    dataset = gdal.Open(fname) 
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    count = dataset.RasterCount
    geotrans = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    # convert format WKT -> Proj4
    spatialref = osr.SpatialReference()
    spatialref.ImportFromWkt(projection)
    spatialproj = spatialref.ExportToProj4()

    # get bands
    bands = [dataset.GetRasterBand(k) for k in range(1, count+1)]
    null = bands[0].GetNoDataValue()

    # scale to downsample
    if max_dim:
        scale_cols = cols / max_dim
        scale_rows = rows / max_dim
        scale_max = max(scale_cols, scale_rows)
        nx = round(cols / scale_max)
        ny = round(rows / scale_max)
    else:
        nx = cols
        ny = rows

    # read downsampled arrays
    image = [b.ReadAsArray(buf_xsize=nx, buf_ysize=ny) for b in bands]

    # rgb (3d) or grayscale (2d) image
    if count > 1:
        image = np.dstack(image)
    else:
        image = image[0]

    # resolution, pixel with and height
    # origin, raster upper left cornner
    # rotation, 0 if image is "north up"
    dx = geotrans[1]
    dy = geotrans[5]
    x_left = geotrans[0]
    y_top = geotrans[3]
    rot1 = geotrans[2]
    rot2 = geotrans[4]

    # lower right cornner: http://gdal.org/gdal_datamodel.html
    x_right = x_left + cols * dx + rows * rot1
    y_bottom = y_top + cols * rot2 + rows * dy 

    print "WKT format:\n", spatialref
    print "Proj4 format:\n", spatialproj
    print 'res:', dx, dy
    print 'xlim:', x_left, x_right
    print 'ylim:', y_top, y_bottom
    print 'rot:', rot1, rot2
    return image, (x_left, y_top, x_right, y_bottom)


# DEPRECATED - use read_gtif instead.
def get_gtif(fname, lat_ts=-71, lon_0=0, lat_0=-90, units='m'):
    """Reads a GeoTIFF image and returns the respective 2d array.

    It assumes polar stereographic projection.
    If units='km', converts x/y from m to km.

    Return
    ------
    x, y : 1d arrays containing the coordinates.
    img : 2d array containing the image.
    bbox_ll : lon/lat limits (lllon, lllat, urlon, urlat).

    Notes
    -----
    MOA parameters: http://nsidc.org/data/moa/users_guide.html
    lat_ts is standard lat (lat of true scale), or "latitude_of_origin".
    lon_0/lat_0 is proj center (NOT grid center!).
    The MOA image has x/y units in m.

    To get GeoTIFF metadata:

        $ gdalinfo file.tif
    
    """
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
    ymin = ymax + nx * gt[4] + ny * dy 
    xmax = xmin + nx * dx + ny * gt[2] 

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
    if units == 'km' and is_meter(x):  # if coords in m, convert to km
        x /= 1e3
        y /= 1e3
    print 'image limits (left/right/bottom/top):'
    print '(lon,lat)', bbox_ll[0], bbox_ll[2], bbox_ll[1], bbox_ll[3]
    return [x, y, img, bbox_ll]


def read_cindex(fname, from_year=1992, to_year=2012, missing_value=-999,
                comments='#', pandas=False, name=None, round_index=6):
    """Read ascii climate-index table into a time series.
    
    pandas=False - (defaul) returns a list with arrays [t,y]
    pandas=True - returns a pandas Series
    name - pandas Series name (column in DataFrame)
    round_index - round to n decimals the Series index
    """
    table = np.loadtxt(fname, comments=comments)
    table[table==missing_value] = np.nan
    t = np.arange(table[0,0], table[-1,0]+1, 1/12.) 
    y = table[:,1:].flatten()
    ind, = np.where((t >= from_year) & (t <= to_year))
    s = [t[ind], y[ind]]
    if pandas:
        s = pd.Series(s[1], index=np.round(s[0], round_index), name=name)
    return s


def write_slabs(fid, name, data, group=None):
    """Save 3d array (to HDF5) into several slabs (in the 0-axis) for XDMF."""
    if group is None:
        group = 'data'
    g = fid.create_group('/', group)
    for i, d in enumerate(data):
        fid.create_array(g, name +'_%02d' % i, d)
