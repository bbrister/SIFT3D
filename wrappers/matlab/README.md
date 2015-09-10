# SIFT3D for Matlab

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

## Contents

This is a Matlab toolbox wrapping SIFT3D. It contains the following functions:
- detectSift3D.m - Detect SIFT3D keypoints from a 3-dimensional array.
- extractSift3D.m - Extract SIFT3D descriptors from a 3-dimensional array.

This toobox does not include all of SIFT3D's functionality, but we are adding more soon.

## Installation instructions

If Matlab is detected in your system, a Matlab toolbox is compiled in the build/wrappers/matlab subdirectory during the normal CMake build process.

To use the toolbox, simply add it to your Matlab path. This can be accomplished by adding the following line to your startup.m file:

        run('/path/to/build/wrappers/matlab/setupSift3D')

where /path/to/build is the path to your SIFT3D build directoy. If you do not have a startup.m file, you can simply run this command prior to calling any SIFT3D functions.

### Relocating the toolbox

We do not recommend moving the toolbox libraries (.a, .so, .dylib). If you do, the operating system may not be able to find them when it tries to load the mex files. A better solution is to compile SIFT3D in the place you want the toolbox to reside. See the "Out-of-source build" section of the main README for more information.

### Known issues

CMake might fail to find Matlab on your system, especially on Mac OSX. In that case, you should see "Matlab not found" after running the cmake command. You can fix this by manually specifying the location of Matlab in the variable Matlab_ROOT_DIR. For example,

        cmake .. -DMatlab_ROOT_DIR=/path/to/MATLAB

or on Mac OSX,

        cmake .. -DMatlab_ROOT_DIR=/path/to/Matlab.app

where /path/to/MATLAB is the location of your Matlab installation.

## Usage instructions

For instructions on using the Matlab functions, use the help page, e.g.

        help detectSift3D

See /examples for sample programs.

## Advanced features

This toolbox also includes a test suite, Sift3DTest.m. It requires [xUnit](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework) to run. You can run the test suite with the following command:

        runtests
