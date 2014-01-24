======================================================
 Processing scheme for time series of elevation change 
======================================================

General
-------
1. read raw IDR/WDR (parallel)... ok
2. convert to HDF5 (parallel)... ok
3. apply mask (parallel)... ok
4. sepparate tracks by flags (parallel)... ok
5. separate time: years, seasons, months... ok
6. separate in sectors w/overlap in lon (fast operation)... ok

Notes
-----
keep separated files (no merge) as long as possible.
The more (and smaller) files the better for parallel processing!

Specific
--------
9.  filter data: shelf/land/ocean/buff & GLA specific [filtflags.py, filtgla.py]... ok
10. merge times (check if it's too much data) [mpi_merge.py]... ok
11. separate asc/des FLAGS in different files for 'x2sys' (much faster!) [filttracks.py]... ok
12. HDF5 to BIN for 'x2sys' (fast operation)... ok
13. find crossovers [mpi_ntasks.py, mpi_x2sys.py]. Warning: lon -/+180 is changed to 0/360!... ok
14. predict tides: OTPSm, check struct idr/gla! (on Triton)... ok
15. merge ad/da crossovers [mpi_merge.py, mpi_merge2.py] code_ok... ok
16. average bins (and remove sector overlaps): xovers -> grid (read code header) [x2grid.py, run.py], [1] code_ok... ok
17. join all subgrids in time/space (serial) [gridjoin.py] code_ok... ok
[18. apply ICESat biases!]

18. average time series: grid-cells (serial; read code header) [averagets2.py] code_ok
19. independent backscatter correction
20. merge satellites (serial, extremely fast) [mergesats.py] code_ok
[interpolate before or after cross-calib??? -> do not use interpolated values to cross-calib!]
21. cross calibration
22. join backscatter correction

# 22. constant backscatter correction
# 23. tvariable backscatter correction

[1] do not apply `iterative` std to ICESat (and perhaps to any sat)


# NOTE: before using Xgrid check DIR and FILE permissions!


#===========================
# read IDR/WDR (bin -> txt)
#===========================
#python batch.py -j read-seasat -s -c '/Users/fpaolo/code/idr/readidr_ra1 -i 3 -v -d /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/ID04/*.ID04
#python batch.py -j read-gfo -s -c '/Users/fpaolo/code/idr/readidr_ra1 -i 2 -v -d /Volumes/LaCie1TB/ANTARCTICA/GFO/work' /Volumes/LaCie1TB/ANTARCTICA/GFO/ID04/*.ID04
#python batch.py -j read-geogm -s -c '/Users/fpaolo/code/idr/readidr_ra1 -i 2 -v -d /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/ID04/*.ID04
#python batch.py -j read-geoerm -s -c '/Users/fpaolo/code/idr/readidr_ra1 -i 3 -v -d /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/ID04/*.ID04
#python batch.py -j read-ers1 -s -c '/Users/fpaolo/code/idr/readidr_ra1 -i 3 -v -d /Volumes/LaCie1TB/ANTARCTICA/ERS1/work' /Volumes/LaCie1TB/ANTARCTICA/ERS1/ID05/*.ID05
#python batch.py -j read-ers1-2 -s -c '/data/xgrid/code/idr/readidr_ra1 -i 3 -v -d /data/xgrid/ers1/raw' /data/xgrid/ers1/raw/*.ID05
#python batch.py -e fspaolo@gmail.com -j read-ers2-2 -s -c '/data/xgrid/code/idr/readidr_ra1 -i 3 -v -d /data/xgrid/ers2/raw' /data/xgrid/ers2/raw/*.ID05
#python batch.py -e fspaolo@gmail.com -j read-envi -s -c '/data/xgrid/code/idr/readidr_ra2 -v -d /data/xgrid/envi/raw/antarctica_new' /data/xgrid/envi/raw/antarctica_new/*.ID05
###python batch.py -e fspaolo@gmail.com -j read-envi-bin -s -c '/data/xgrid/code/idr/readidr_ra2 -v -b -d /data/xgrid/envi/raw/antarctica_new' /data/xgrid/envi/raw/antarctica_new/*.ID05

# NOTE: #PBS -l walltime=0:40:00


#======================
# compress ASCII files
#======================
# TODO

#=================
# convert to HDF5
#=================
#python batch.py -j convert-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/txth5.py -v -l blosc' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*.ID04.txt
#python batch.py -j convert-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/txth5.py -v -l blosc' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*.ID04.txt
#python batch.py -j convert-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/txth5.py -v -l blosc' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*.ID04.txt
#python batch.py -j convert-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/txth5.py -v -l blosc' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*.ID04.txt
#python batch.py -j convert-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/txth5.py -v -l blosc' /data/xgrid/ers1/raw/*.ID05.txt
#python batch.py -j convert-ers2-2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/txth5.py -v -l blosc' /data/xgrid/ers2/raw/*.ID05.txt
#python batch.py -j convert-envi -e fspaolo@gmail.com -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/txth5.py -v -l blosc' /data/xgrid/envi/raw/antarctica_new/*.ID05.txt


#============
# apply mask 
#============
# note: see file extension !!!!!
#------------
#python xg-batch.py -j mask-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/mask/runmask.py -f /Users/fpaolo/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*.ID04.h5
#python xg-batch.py -j mask-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/mask/runmask.py -f /Users/fpaolo/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*.ID04.h5
#python xg-batch.py -j mask-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/mask/runmask.py -f /Users/fpaolo/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*.ID04.h5
#python xg-batch.py -j mask-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/mask/runmask.py -f /Users/fpaolo/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*.ID04.h5
#python xg-batch.py -j mask-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/mask/runmask.py -f /data/xgrid/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /data/xgrid/ers1/raw/*.ID05.h5
#python xg-batch.py -j mask-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/mask/runmask.py -f /data/xgrid/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /data/xgrid/ers2/raw/*.ID05.h5
#python xg-batch.py -j mask-envi -s -e fspaolo@gmail.com -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/mask/runmask.py -f /data/xgrid/code/mask/mask_ice_1km_2008_0410c.h5 -b 3' /data/xgrid/envi/raw/*.ID05.h5

# using Scripps Mask:
# maskrun.py -f /home/fpaolo/code/mask/scripps_antarctica_mask1km_v1.h5

# for RA
# NOTE: #PBS -l walltime=0:25:00

# for GLA
# NOTE: #PBS -l walltime=0:10:00



#=================
# separate tracks
#=================
# note: extremely fast (C module)!
#-----------------
#python xg-batch.py -j tracksep-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/tracksep.py -m' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_mask.h5
#python xg-batch.py -j tracksep-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/tracksep.py -m' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_mask.h5
#python xg-batch.py -j tracksep-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/tracksep.py -m' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_mask.h5
#python xg-batch.py -j tracksep-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/tracksep.py -m' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_mask.h5
#python xg-batch.py -j tracksep-amery-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/tracksep.py -m -f' /data/xgrid/ers1/proc/*_amery_shelf.h5
#python xg-batch.py -j tracksep-amery-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/tracksep.py -m -f' /data/xgrid/ers2/proc/*_amery_shelf.h5
#python xg-batch.py -j tracksep-amery-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/tracksep.py -m -f' /data/xgrid/envi/proc/*_amery_shelf.h5


#==========
# add bias << check why the code is so slow! --> because floating and *grounded* ice points!
#==========
# note1: see file extension!
# note2: this is the *time and spatialy invariant* bias
# note3: do not *apply* bias, add extra column!
#----------
#python xg-batch.py -j bias-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/bias/increment.py -i 0.32' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_sep.h5
#python xg-batch.py -j bias-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/bias/increment.py -i 0.11' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_sep.h5
#python xg-batch.py -j bias-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/bias/increment.py -i 0.11' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_sep.h5
#python xg-batch.py -j bias-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/bias/increment.py -i 9.99' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_sep.h5
#python xg-batch.py -j bias-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/bias/increment.py -i -0.45' /data/xgrid/ers1/raw/*_sep.h5
#python xg-batch.py -j bias-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/bias/increment.py -i -0.05' /data/xgrid/ers2/raw/*_sep.h5
#python xg-batch.py -j bias-envi -s -e fspaolo@gmail.com -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/bias/increment.py -i -0.45' /data/xgrid/envi/raw/*_sep.h5


# Obs: ATTENTION WITH PIG at this point !!!
# For PIG region see bottom of the script


#==========================================================
# separate time: YEAR (default) / MONTH (-m) / SEASON (-s)
#==========================================================
#python xg-batch.py -j seasons-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/timesep.py -s' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_bias_????.h5
#python xg-batch.py -j seasons-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/timesep.py -s' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_bias_????.h5
#python xg-batch.py -j seasons-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/timesep.py -s' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_bias_????.h5
#python xg-batch.py -j seasons-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/separate/timesep.py -s' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_bias_????.h5
#python xg-batch.py -j seasons-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/timesep.py -s' /data/xgrid/ers1/raw/*_bias.h5
#python xg-batch.py -j seasons-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/timesep.py -s' /data/xgrid/ers2/raw/*_bias.h5
#python xg-batch.py -j seasons-envi -e fspaolo@gmail.com -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/separate/timesep.py -s' /data/xgrid/envi/raw/*_bias.h5


#=====================
# separate in sectors
#=====================
#python batch.py -j region3-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_bias_????_????.h5
#python batch.py -j region3-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_bias_????_????.h5
#python batch.py -j region3-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_bias_????_????.h5
#python batch.py -j region3-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_bias_????_????.h5
#python batch.py -j region3-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/ERS1/work/*_bias_????_????.h5
#python batch.py -j region3-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/ERS2/work/*_bias_????_????.h5
#python batch.py -j region3-envi -e fspaolo@gmail.com -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0' /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/*_bias_????_????.h5


#===============
# select region
#===============
# note: the data base has already applied `mask`, `sep`, `bias`, `year`/`season`/`month`
# so select region from files: *_mask_sep_bias_??????.h5
#---------------
#python batch.py -j region-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*.ID04.h5
#python batch.py -j region-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*.ID04.h5
#python batch.py -j region-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*.ID04.h5
#python batch.py -j region-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*.ID04.h5
#python batch.py -j region-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/ERS1/work/*.ID05.h5
#python batch.py -j region-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/ERS2/work/*.ID05.h5
#python batch.py -j region-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/*.ID04.h5
#------
# FRIS
#------
#python xg-batch.py -j region-fris-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 270 340 -82 -74 -s _fris -v' /data/xgrid/ers1/raw/*_??????.h5
#python xg-batch.py -j region-fris-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 270 340 -82 -74 -s _fris -v' /data/xgrid/ers2/raw/*_??????.h5
#python xg-batch.py -e fspaolo@gmail.com -j region-fris-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 270 340 -82 -74 -s _fris -v' /data/xgrid/envi/raw/*_??????.h5
#------
# ROSS
#------
#python xg-batch.py -j region-ross-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 150 220 -82 -75 -s _ross -v' /data/xgrid/ers1/raw/*_??????.h5
#python xg-batch.py -j region-ross-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 150 220 -82 -75 -s _ross -v' /data/xgrid/ers2/raw/*_??????.h5
#python xg-batch.py -j region-ross-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 150 220 -82 -75 -s _ross -v' /data/xgrid/envi/raw/*_??????.h5
#-------
# AMERY
#-------
#python xg-batch.py -j region-amery-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 64 80 -75 -67 -s _amery -v' /data/xgrid/ers1/proc/*_bias_??????.h5
#python xg-batch.py -j region-amery-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 64 80 -75 -67 -s _amery -v' /data/xgrid/ers2/proc/*_bias_??????.h5
#python xg-batch.py -j region-amery-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/select/region.py -r 64 80 -75 -67 -s _amery -v' /data/xgrid/envi/proc/*_bias_??????.h5


# region3.py -r 0 360 -90 90 -d 10 180 -l 1 1 0 0 *_bias.h5


#=============
# merge files 
#=============
# note1: change the pattern!
# note2: very fast!
#-------------
#python xg-batch.py -j merge-season-seasat -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/merge2.py -p _\d\d\d\d_\d\d\d\d_reg\d\d -o seasat_' '/Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_mask_sep_bias_*_reg??.h5'
#python xg-batch.py -j merge-season-geogm -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/merge2.py -p _\d\d\d\d_\d\d\d\d_reg\d\d -o geogm_' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_reg??.h5
#python xg-batch.py -j merge-season-geoerm -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/merge2.py -p _\d\d\d\d_\d\d\d\d_reg\d\d -o geoerm_' '/Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_reg??.h5'
#python xg-batch.py -j merge-season-gfo -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/merge2.py -p _\d\d\d\d_\d\d\d\d_reg\d\d -o gfo_' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_reg??.h5
#python xg-batch.py -j merge-amery-ers1 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p _\d\d\d\d\d\d_amery -o ers1_' /data/xgrid/ers1/proc/*_??????_amery.h5
#python xg-batch.py -j merge-amery-ers2 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p _\d\d\d\d\d\d_amery -o ers2_' /data/xgrid/ers2/proc/*_??????_amery.h5
#python xg-batch.py -j merge-amery-envi -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p _\d\d\d\d\d\d_amery -o envi_' /data/xgrid/envi/proc/*_??????_amery.h5
#
#python xg-batch.py -j merge-ross-ers1 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o ers1_ -s _ross' /data/xgrid/ers1/ross/ers1_??????_??????_ross_*.h5
#python xg-batch.py -j merge-ross-ers2 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o ers2_ -s _ross' /data/xgrid/ers2/ross/ers2_??????_??????_ross_*.h5
#python xg-batch.py -j merge-ross-envi -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o envi_ -s _ross' /data/xgrid/envi/ross/envi_??????_??????_ross_*.h5

#python xg-batch.py -j merge-fris-ers1 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o ers1_ -s _fris' /data/xgrid/ers1/fris/ers1_??????_??????_fris_*.h5
#python xg-batch.py -j merge-fris-ers2 -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o ers2_ -s _fris' /data/xgrid/ers2/fris/ers2_??????_??????_fris_*.h5
#python xg-batch.py -j merge-fris-envi -s -1 -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/merge2.py -p \d\d\d\d\d\d_\d\d\d\d\d\d -o envi_ -s _fris' /data/xgrid/envi/fris/envi_??????_??????_fris_*.h5


#==============
# filter flags 
#==============
# note1: filter ice-shelf to reduce computing time (considerably)!
# note2: atention with GEOSAT_GM!!!
#--------------
#python xg-batch.py -j filter-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/filter/filterflags.py' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_sep.h5
#python xg-batch.py -j filter-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/filter/filterflags_geogm.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_sep.h5
#python xg-batch.py -j filter-geoerm-fris -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/filter/filterflags.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_sep.h5
#python xg-batch.py -j filter-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/filter/filterflags.py' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_sep.h5
#python xg-batch.py -j filter-amery-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/filter/filterflags.py' /data/xgrid/ers1/proc/ers1_??????_amery.h5
#python xg-batch.py -j filter-amery-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/filter/filterflags.py' /data/xgrid/ers2/proc/ers2_??????_amery.h5
#python xg-batch.py -j filter-amery-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/filter/filterflags.py' /data/xgrid/envi/proc/envi_??????_amery.h5


#=====================
# convert HDF5 -> BIN 
#=====================
#python xg-batch.py -j h5tobin-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/h5bin.py' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/seasat_*
#python xg-batch.py -j h5tobin-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/h5bin.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/geogm_*
#python xg-batch.py -j h5tobin-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/h5bin.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_*
#python xg-batch.py -j h5tobin-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/h5bin.py' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/gfo_*
#python xg-batch.py -j h5tobin-amery-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/h5bin.py' /data/xgrid/ers1/proc/ers1_*_amery_shelf_*.h5
#python xg-batch.py -j h5tobin-amery-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/h5bin.py' /data/xgrid/ers2/proc/ers2_*_amery_shelf_*.h5
#python xg-batch.py -j h5tobin-amery-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/convert/h5bin.py' /data/xgrid/envi/proc/envi_*_amery_shelf_*.h5


#==================
# rename the files 
#==================
# note: use `rename_files.py`
#------------------
#python xg-batch.py -j renameh5-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/seasat_??????_??*
#python xg-batch.py -j renameh5-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/geogm_????_????_reg*
#python xg-batch.py -j renameh5-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_????_????_reg*
#python xg-batch.py -j renameh5-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/gfo_????_????_reg*
#python xg-batch.py -j renameh5-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/ERS1/work/ers1_????_????_reg*
#python xg-batch.py -j renameh5-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/misc/rename_files.py' /Volumes/LaCie1TB/ANTARCTICA/ERS2/work/ers2_????_????_reg*
#python xg-batch.py -e fspaolo@gmail.com -j rename-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/misc/rename_files.py' /data/xgrid/envi/raw/antarctica_new/*_bias??????.h5


#=======
# x2sys
#=======
# on Xgrid
#---------
# Note1: --> use `runbatch.py` to generate *multiple* batch files 
# Note2: --> use `find_maxpts.py` to find the best reference time
#-------
#python xg-batch.py -e fspaolo@gmail.com -j x2sys-bench-all -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/x2sys/x2sys.py' /data/xgrid/test/ers1_199207_fris_shelf /data/xgrid/test/ers1_199310_fris_shelf
#python xg-batch.py -e fspaolo@gmail.com -j x2sys-bench-sep1 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/x2sys/x2sys.py' /data/xgrid/test/ers1_199207_fris_shelf_asc /data/xgrid/test/ers1_199310_fris_shelf_des
#python xg-batch.py -e fspaolo@gmail.com -j x2sys-bench-sep2 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /data/xgrid/code/x2sys/x2sys.py' /data/xgrid/test/ers1_199207_fris_shelf_des /data/xgrid/test/ers1_199310_fris_shelf_asc

# on Triton
#----------
# ???


#========================
# complib: blosc -> zlib
#========================
# Not needed if using OPTS Fortran code !!!
# note: quite fast!
#------------------------
#python xg-batch.py -j complib-season-seasat  -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/envi_2008_0608_reg00-*.h5
#python xg-batch.py -j pig-complib-season-oc-land-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /Volumes/LaCie1TB/ANTARCTICA/ERS2/work/envi_pig_2008_0608-*.h5
#python xg-batch.py -j pig-complib-season-oc-land-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/envi_pig_2008_0608-*.h5
## ice-mode
#python xg-batch.py -j complib-fris-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/ers1/fris/*_??????_??????_fris.h5
#python xg-batch.py -j complib-fris-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/ers2/fris/*_??????_??????_fris.h5
#python xg-batch.py -j complib-fris-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/envi/fris/*_??????_??????_fris.h5
#python xg-batch.py -j complib-ross-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/ers1/ross/*_??????_??????_ross.h5
#python xg-batch.py -j complib-ross-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/ers2/ross/*_??????_??????_ross.h5
#python xg-batch.py -j complib-ross-envi -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/convert/complib.py -v' /data/xgrid/envi/ross/*_??????_??????_ross.h5


#=======================
# tide/load corrections
#=======================
# Fortran: see /Users/fpaolo/code/tide/OTPSm/tide.py
# Matlab: see /Users/fpaolo/code/tide/tmd_toolbox/runtide_name.m
#-----------------------


#=============================
# average bins: xover -> grid
#=============================
# see /Users/fpaolo/code/gridding/xover2grid.py
# see /Users/fpaolo/code/gridding/xover2box.py
# --> extremely fast <--
#-----------------------------
# FRIS: -r -100 -20 -82 -75
# ROSS: -r 150 220 -82 -76
# AMERY: -r 66 76 -74 -68


#===========================
# average time series: grid
#===========================
# see /Users/fpaolo/code/tseries/averagets.py
#---------------------------

# CONTINUE EDITING THE CODE FROM HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#====================================
# independent backscatter correction
#====================================
# --> extremely fast <--
#------------------------------------


#==================
# merge satellites
#==================
# see /Users/fpaolo/code/tseries/mergesats.py
#------------------


#===================
# cross calibration
#===================
# see /Users/fpaolo/code/tseries/crosscalib.py
#-------------------


#========================
# backscatter correction
#========================
# see /Users/fpaolo/code/backscatter/backscatter.py
# --> extremely fast <--
#------------------------


#==============
# plot results
#==============
# see /Users/fpaolo/code/tseries/plot.py
#--------------


###### END OF PROCESSING ##############################################



#============
# test Xgrid
#============
#python xg-batch.py -j x2sys-test1 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198801_26 /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198804_26
#python xg-batch.py -j x2sys-test2 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198801_26 /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198807_26 
#python xg-batch.py -e fspaolo@gmail.com -j x2sys-test3 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198804_26 /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_198807_26
#python xg-batch.py -j x2sys-test1 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /data/alt/ra/antarctica/geosaterm/work/geoerm_198801_26 /data/alt/ra/antarctica/geosaterm/work/geoerm_198804_26
#python xg-batch.py -j x2sys-test2 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /data/alt/ra/antarctica/geosaterm/work/geoerm_198801_26 /data/alt/ra/antarctica/geosaterm/work/geoerm_198807_26 
#python xg-batch.py -e fspaolo@gmail.com -j x2sys-test3 -1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py' /data/alt/ra/antarctica/geosaterm/work/geoerm_198804_26 /data/alt/ra/antarctica/geosaterm/work/geoerm_198807_26


#------------------------------------------------------------------------


#python xg-batch.py -j x2sys-reg$j-seasat  -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py -r /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/envi_2008_0608_reg' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/seasat_????_????_reg$j

#python xg-batch.py -j x2sys-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py -r /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/envi_200807' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/geogm_??????
#python xg-batch.py -j x2sys-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py -r /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/envi_200807' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/geoerm_??????
#python xg-batch.py -j x2sys-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/x2sys/x2sys.py -r /Volumes/LaCie1TB/ANTARCTICA/ERS1/work/ers1_199407' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/gfo_??????



### region: PIG-shelf = -103(257)/-92(268)/-76/-74, PIG-land = -103(257)/-90(270)/-77/-74

#python batch.py -j pig-region-seasat -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 268 -76 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/SEASAT/work/*_bia.h5
##python batch.py -j pig-region-gfo -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 0 360 -82 -60 -v' /Volumes/LaCie1TB/ANTARCTICA/GFO/work/*_bias.h5
#python batch.py -j pig-region-geogm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 268 -76 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_GM/work/*_bias.h5
#python batch.py -j pig-region-geoerm -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 268 -76 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/GEOSAT_ERM/work/*_bias.h5

#python batch.py -j pig-region-land-ers1 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 270 -77 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/ERS1/work/*_bias.h5
#python batch.py -j pig-region-land-ers2 -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 270 -77 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/ERS2/work/*_bias.h5
#python batch.py -j pig-region-land-envi -e fspaolo@gmail.com -s -c '/Library/Frameworks/EPD64.framework/Versions/Current/bin/python /Users/fpaolo/code/select/region.py -r 257 270 -77 -74 -e _pig -v' /Volumes/LaCie1TB/ANTARCTICA/ENVISAT/work/*_bias.h5
