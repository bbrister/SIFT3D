% registerMatlab
%
% Example of image registration in Matlab
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Load the images
[src, srcUnits] = imRead3D('data/1.nii.gz');
[ref, refUnits] = imRead3D('data/2.nii.gz');

% Register
[A, matchSrc, matchRef] = registerSift3D(src, ref, 'srcUnits', ...
    srcUnits, 'refUnits', refUnits);

% Clear MEX memory
clear mex