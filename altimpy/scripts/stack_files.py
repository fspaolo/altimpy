"""
Merge the output-files of ParaView's selection tool.

Stack several single-time-series csv-files into two txt-files: 
1) a file containing all time series as columns (data matrix): *_mat.txt
2) a file containing all time series as single column (data vector): *_vec.txt

Notes 
-----
- center (mean = 0) and normalize (std = 1) each time series.
- for (2), window each time series to assure continuity in concatenation. 
- (1) is suited for MSSA[1] and (2) is suited for MTM/MEM[2].

[1] multi-channel singular spectrum analysis
[2] multi-taper/maximum entropy methods

Selecting from ParaView
-----------------------
Select Cells With Polygon > Plot Selection Over Time > Save Data

Usage
-----
$ python stack_files.py *.csv

Fernando Paolo <fpaolo@ucsd.edu>
Dec 1, 2013
"""

import os
import sys
import re
import numpy as np
import pandas as pd
import altimpy as ap
import statsmodels.api as sm
from scipy import signal
from glob import glob

def window(x):
    return x * signal.windows.hanning(len(x))

def normalize(x):
    x -= x.mean()  # center (mean = 0)
    x /= x.std()   # normalize (std = 1)
    return x

def hpfilter(y, lamb=7):
    return sm.tsa.filters.hpfilter(y, lamb=lamb)[1]

def gradient(y, dt=.25):
    return np.gradient(y.values, dt)


files = sys.argv[1:]
#files = glob('/Users/fpaolo/data/shelves/all/*.csv')
files.sort(key=ap.human_order)

# take root of file name (no numbers)
root = os.path.splitext(''.join(c for c in files[0] if not c.isdigit()))[0]
fname1 = root + '_mat.txt'
fname2 = root + '_vec.txt'

# matrix file
#------------
'''
THIS PROVIDES THE BEST RESUT FOR MSSA (with full seasonal data)
Using the following parameters:
sampling: 0.25, Varimax: yes, PCA channels: 10, window: 26, components: 10
Obs: do not detrend!
'''
d = pd.concat([pd.read_csv(f, usecols=[0]) for f in files], axis=1).dropna(1, 'all')
#d = d.apply(hpfilter)
#d = d.apply(gradient)
d = d.apply(normalize)
d.to_csv(fname1, sep=' ', header=False, index=False, float_format='%.6f')
#d.to_csv(root+'_grad_mat.txt', sep=' ', header=False, index=False, float_format='%.6f')
'''
d = d.mean(axis=1)
d = normalize(d)     
d.to_csv(root+'_mean_vec.txt', sep=' ', header=False, index=False, float_format='%.6f')
'''

# vector file
#------------
d = d.apply(window)  # window each time series
d = d.unstack()      # concatenate all time series
d = normalize(d)     # normalize concatenated series
#d = pd.Series(signal.detrend(d))
d.to_csv(fname2, header=False, index=False, float_format='%.6f')

print 'stacked files ->', fname1, fname2
