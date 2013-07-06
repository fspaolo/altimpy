
AltimPy - Set of tools for processing satellite altimetry data
==============================================================

A Python package with high-level functions for constructing time 
series of satellite altimetric measurements (surface elevation and
backscatter).

The package provides routines to read the raw data from the binary
format Ice Data Records (IDRs), and write/read the processed results 
to the high-performance file format HDF5; as well as processing 
algorithms for filtering, crossover finding, griding, uncertainty
estimation, parallel computing, visualization, etc.

Functions are divided into the following modules:

* `const.py` - definition of constants
* `convert.py` - time/format/etc. conversion functions
* `io.py` - read and write routines
* `util.py` - miscellaneous utility functions
* `viz.py` - visualization functionalities
* `mpi.py` - simple parallelization functions (MPI)

More coming soon...
