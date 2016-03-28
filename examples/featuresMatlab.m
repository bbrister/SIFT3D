% featuresMatlab.m
%
% This script shows how to extract SIFT3D features from a volumetric image.
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Load the image
[im, units] = imRead3D('data/1.nii.gz');

% Detect keypoints
keys = detectSift3D(im, 'units', units);

% Extract descriptors
[desc, coords] = extractSift3D(keys);

% Clear MEX memory
clear mex