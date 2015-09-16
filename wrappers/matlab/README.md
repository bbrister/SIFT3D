# SIFT3D for Matlab

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

## Contents

This is a Matlab toolbox wrapping SIFT3D. It contains the following functions:
- detectSift3D.m - Detect SIFT3D keypoints from a 3-dimensional array.
- extractSift3D.m - Extract SIFT3D descriptors from a 3-dimensional array.

This toobox does not include all of SIFT3D's functionality, but we are adding more soon.

## Installation instructions

If you installed SIFT3D from binaries, a Matlab toolbox is included in the lib/wrappers/matlab subdirectory of your installation. If you compiled SIFT3D from souce, the same will happen if Matlab was detected in your system.

To use the toolbox, simply add it to your Matlab path. This can be accomplished by adding the following line to your startup.m file:

        run('/path/to/sift3d/lib/wrappers/matlab/setupSift3D')

where /path/to/sift3d is the path to your SIFT3D installation. If you do not have a startup.m file, you can simply run this command prior to calling any SIFT3D functions.

### Relocating the toolbox

We do not recommend moving the toolbox libraries (.a, .so, .dylib). If you do, the operating system may not be able to find them when it tries to load the mex files. A better solution is to install SIFT3D in the place you want the toolbox to reside.

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

This toolbox also includes a test suite, Sift3DTest.m. The test suite is only found in the source distribution, not the binary distributions, and must be run from the build tree. It requires [xUnit](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework) to run. You can run the test suite with the following command:

        runtests