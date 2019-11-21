"""
    io.py
    ---------------------------------------------------------------------------
    Copyright (c) 2019 Blaine Rister et al., see LICENSE for details.
    ---------------------------------------------------------------------------
    Example of using the Python wrappers for I/O
    ---------------------------------------------------------------------------
"""

import os
from pySift3D import imutil

dataDir = 'data'
niiInName = os.path.join(dataDir, '1.nii.gz')

im, units = imutil.im_read(niiInName)
