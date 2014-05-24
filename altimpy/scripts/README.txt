=====================================================
Processing scheme for time series of elevation change 
=====================================================

Notes
-----
keep separated files (do no merge) as long as possible.
The more (and smaller) files the better for parallel (and serial) processing!

General processing
------------------
1. read raw IDR/WDR (parallel)... ok
2. convert to HDF5 (parallel)... ok
3. apply mask (parallel)... ok
4. sepparate tracks by flags (parallel)... ok
5. separate time: years, seasons, months... ok
6. separate in sectors w/overlap in lon (fast operation)... ok

Specific processing
-------------------
9.  filter data: shelf/land/ocean/buff & GLA specific [filtflags.py, filtgla.py]... ok
10. merge times (check if it's too much data) [mpi_merge.py]... ok
11. separate asc/des FLAGS in different files for 'x2sys' (much faster!) [filttracks.py]... ok
12. HDF5 to BIN for 'x2sys' (fast operation)... ok
13. find crossovers [mpi_ntasks.py, mpi_x2sys.py]. Warning: lon -/+180 is changed to 0/360!... ok
14. predict tides: OTPSm, check struct idr/gla! (on Triton)... ok
15. merge ad/da crossovers [mpi_merge.py, mpi_merge2.py] code_ok... ok
16. average bins (and remove sector overlaps): xovers -> grid (read code header) [x2grid.py], [1] code_ok... ok
17. join all subgrids in time/space (serial, slow) [gridjoin.py] code_ok... ok
[apply ICESat biases here!]
18. average time series: grid-cells (serial, very slow; read code header) [averagets2.py] code_ok... ok
19. independent backscatter correction (serial, fast) [backscatter3.py] code_ok... ok
20. merge satellites (serial, extremely fast) [mergesats2.py] code_ok... ok
21. cross calibration (serial, slow) [crosscalib2.py] code_ok... ok
22. error/nobs propagation (serial, slow) [crosscalib3.py] code_ok... ok
23. joint backscatter correction (serial, slow) [backscatter2.py] code_ok... ok

24. [no!] interpolate empty bins [kriging.py] (but do not use interpolated values to cross-calib)!
25. post processing (reference, spherical->cartesian, filter, regrid...) [post_proc.py]

# constant backscatter correction
# tvariable backscatter correction

Visualization
-------------
1. create file with multiple layers (2d arrays, for time dimension) [post_proc.py] (for XDMF)
2. create XDMF for data/error/coord to visualize in ParaView [write_xdmf.py]


[1] do not apply 'iterative' std to ICESat (and perhaps to any sat)
