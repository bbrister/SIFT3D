% Load the data
load data/data.mat

% Detect keypoints
keys = detectSift3D(im1);

% Extract descriptors
[desc, coords] = extractSift3D(keys);
