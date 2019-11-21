"""
    sift3d.py
    ---------------------------------------------------------------------------
    Copyright (c) 2019 Blaine Rister et al., see LICENSE for details.
    ---------------------------------------------------------------------------
    Python wrappers for SIFT3D I/O routines
    ---------------------------------------------------------------------------
"""

import numpy as np
import ctypes

import pySift3D

def __im_to_np__(im):
    """
        Convert an Image class to a numpy array
    """

    # Extract the dimensions and units
    dims = np.array([int(im.dims[i]) for i in xrange(pySift3D.imNDims)]) 
    nc = int(im.nc)
    units = np.array([float(im.units[i]) for i in xrange(pySift3D.imNDims)])

    # Initialize the data numpy array
    outShape = (nc,) + tuple(dims)
    numel = nc * np.prod(dims)
    imNp = np.ctypeslib.as_array(im.data, shape=(numel,))
    imNp = imNp.reshape(outShape)

    # Reorder from (c, x, y, z) to (x, y, z, c)
    imNp = np.moveaxis(imNp, 0, -1)

    return imNp, units

def im_read(path):
    """
        Read an image given the file path
    """

    # Load the data
    im = pySift3D.imReadFun(path.encode('utf-8'))
    read_error_code = im.error_code

    # Raise error messages
    if read_error_code is not 0:
        raise ValueError(pySift3D.imReadGetErrorFun(read_error_code))

    # Convert to a Numpy array
    imNp, units = __im_to_np__(im) 

    # Copy the data and free the C memory
    imNp = np.require(imNp, 'O')
    pySift3D.imFreeFun(ctypes.byref(im))

    # Return the converted image, as well as units
    return imNp, units
