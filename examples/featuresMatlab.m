% Load the data
im = imRead3D('data/1.nii.gz');

% Detect keypoints
keys = detectSift3D(im);

% Extract descriptors
[desc, coords] = extractSift3D(keys);
