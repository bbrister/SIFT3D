function keys = detectSift3D(im, units)
%detectSift3D(im, units) Detect Sift3D keypoints in an image.
%  Arguments:
%    im - An [MxNxP] array, where voxels are indexed in (x, y, z) order.
%    units - (Optional) See imRead3D. If units are specified, the detected
%       keypoints are isotropic even when im is not.
%
%  Return values:
%    keys - A [Qx1] array of keypoint structs. See keypoint.m for the 
%      struct definition.
%
%  Keypoint coordinates are defined in the space of their pyramid level.
%  To convert them to the input image space, use the following 
%  transformation:
%      key.coords * pow2(-key.octave)
%
%  Examples:
%      im = rand(50, 50, 50);
%      keys = detectSift3D(im);
%
%      [im, units] = imRead3D('someFile.dcm');
%      keys = detectSift3D(im, units);
%
%  See also:
%    extractSift3D, imRead3D, imWrite3D, keypoint3D, setupSift3D
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Default parameters
if nargin < 2
    units = [];
end

% Verify inputs
if nargin < 1
        error('Not enough arguments');
end

if isempty(im)
    error('im is empty')
end

if ndims(im) ~= 3
    error(['im must have 3 dimensions, detected ' num2str(ndims(im))]);
end

units = checkUnits3D(units);

% Scale and convert the image to single precision
im = single(im);
im = im / max(im(:));

% Detect features
keys = mexDetectSift3D(im, units);

end

