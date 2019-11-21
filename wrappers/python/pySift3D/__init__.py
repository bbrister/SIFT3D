
"""
Python module initializer for SIFT3D.
"""

import os
import ctypes

# The number of dimensions in an Image
#TODO could query this from the library, probably not subject to change
imNDims = 3

# Image struct
class Image(ctypes.Structure):

        # C structure fields
        _fields_ = [
                ('data', ctypes.POINTER(ctypes.c_float)), 
                ('dims', ctypes.c_int * imNDims),
                ('units', ctypes.c_double * imNDims),
                ('nc', ctypes.c_int),
                ('error_code', ctypes.c_int)
        ]

# Load the library
libImutilName = 'libsift3Dpython.so'
scriptDir = os.path.abspath(os.path.dirname(__file__))
prefixes = [scriptDir, '/usr/lib', '/usr/local/lib', '/usr/local/lib/sift3d']

# Try to find the library using each of the available prefixes, in order
dll = None
searched = []
for prefix in prefixes:
    searchName = os.path.join(prefix, libImutilName)
    if not os.path.exists(searchName):
        searched.append(searchName)
        continue

    dll = ctypes.cdll.LoadLibrary(searchName)
    break

if dll is None:
   raise OSError('Cannot find library ' + libImutilName + '. Searched the ' +
        'following paths: ' + '\n'.join(searched))

# Extract the im_read function
imReadFun = dll.python_im_read
imReadFun.argtypes = [
        ctypes.c_char_p # Pass an ASCII-encoded string
    ]
imReadFun.restype = Image

# Extract the im_read_get_error function
imReadGetErrorFun = dll.python_im_read_get_error
imReadGetErrorFun.argtypes = [ctypes.c_int]
imReadGetErrorFun.restype = ctypes.c_char_p

# Extract the im_free funtion
imFreeFun = dll.python_im_free
imFreeFun.argTypes = [
        ctypes.POINTER(Image)
    ]
    # Returns void

