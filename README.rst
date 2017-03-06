
AltimPy 
=======

**Set of tools for processing satellite altimetry data**

A Python package with high-level functions for constructing time 
series of satellite altimetric measurements (surface heigh and
backscatter).

The package provides routines to read the raw data from the binary
format Ice Data Records (IDRs), and write/read the processed results 
to the high-performance file format HDF5; as well as processing 
algorithms for filtering, crossover finding, gridding, uncertainty
estimation, time-series construction, parallelization, visualization, 
etc.

Functionalities are divided into the following modules:

* ``io.py`` - read and write routines
* ``const.py`` - definition of constants
* ``convert.py`` - time/coordinates/etc. conversion functions
* ``tseries.py`` - crossover time series construction
* ``interp.py`` - wrapper for kriging and other interpolation routines
* ``filter.py`` - high-level time-series filtering functions
* ``bootstrap.py`` - bootstrap methods for uncertanty estimation
* ``regions.py`` - definition of "boxes" for Antarctic ice shelves
* ``util.py`` - miscellaneous utility functions
* ``viz.py`` - plotting functionalities
* ``mpi.py`` - simple parallelization functions (MPI)


Download and install the library
---------------------------------
::

    git clone https://github.com/fspaolo/altimpy 
    cd altimpy 
    python setup.py install --user 

Some day I will organize and document this library... ;)
