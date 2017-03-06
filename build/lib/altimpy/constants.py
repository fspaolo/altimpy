import numpy as np

CDR = d2r = 180/np.pi         # conversion degrees to radians (180/pi)
E2 = ec2 = 6.694379852*1e-3   # eccentricity squared
E = ec = np.sqrt(E2)          # eccentricity
PI = np.pi 
RE = 6378.1370                # Earth radius, updated 2/11/08 (see email from Shad O'Neel)
#RE = 6378.273                # original value
