"""
Merge the output-files of ParaView's selection tool.

Stack several single-time-series csv-files into two txt-files: 
1) a file containing all time series as columns (data matrix): *_mat.txt
2) a file containing all time series as single column (data vector): *_vec.txt

Notes 
-----
- center (mean = 0) and normalize (std = 1) each time series.
- for (2), window each time series to assure continuity in concatenation. 
- (1) is suited for MSSA and (2) is suited for MTM.

Usage
-----
$ python stack_files.py *.csv

"""

import os
import sys
import re
import numpy as np
import pandas as pd
import altimpy as ap
from scipy import signal

def window(x):
    return x * signal.windows.hanning(len(x))

def normalize(x):
    x -= x.mean()  # center (mean = 0)
    x /= x.std()   # normalize (std = 1)
    return x


files = sys.argv[1:]
files.sort(key=ap.human_order)

# take root of file name (no numbers)
root = os.path.splitext(''.join(c for c in files[0] if not c.isdigit()))[0]
fname1 = root + '_mat.txt'
fname2 = root + '_vec.txt'

# matrix file
d = pd.concat([pd.read_csv(f, usecols=[0]) for f in files], axis=1).dropna(1, 'all')
d = d.apply(normalize)
d.to_csv(fname1, sep=' ', header=False, index=False, float_format='%.6f')

# vector file
#d = d.apply(signal.detrend)
d = d.apply(window)
d = d.unstack()
d.to_csv(fname2, header=False, index=False, float_format='%.6f')

print 'stacked files ->', fname1, fname2
