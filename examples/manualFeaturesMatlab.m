% manualFeaturesMatlab.m
%
% This script shows how to manually define keypoints in a volumetric image,
% and extract SIFT3D feature descriptors at those locations.
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Load the image
[im, units] = imRead3D('data/1.nii.gz');

% Define a keypoint
keys = keypoint3D([100 100 20]);

% Assign its orientation
[keys, conf] = orientation3D(keys, im, units);

% Extract a descriptor
[desc, coords] = extractSift3D(keys, im, units);
