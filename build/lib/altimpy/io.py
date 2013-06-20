import os
import re
import numpy as np
import scipy as sp
import tables as tb
import datetime as dt


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
    """
    Add columns to an existing table.
    """
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

