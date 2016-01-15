%orientation3D(keys, im, units) Assign 3D orientations to keypoints 
%   in an image. This is not needed if your keypoints were detected with
%   detectSift3D, but may be useful if your keypoints were defined
%   manually with keypoint3D.
%
%  Arguments:
%    keys - [Mx1] array of keypoint structs. See keypoint3D.m.
%    im - An [MxNxP] array, where voxels are indexed in (x, y, z) order.
%    units - (Optional) See imRead3D. (Default: [1 1 1])
%
%  Return values:
%    keys - The same as the input, but with the "ori" fields filled with
%       rotation matrices assigned based on the data in im.
%    conf - An [Mx1] vector of confidence scores. The scores are in the 
%       interval [0, 1], where 1 is the most confident, 0 is the least. A 
%       higher score means the assigned orientation is more likely to be 
%       robust. A negative score means the orientation could not be assigned.
%
%  Note: For some data, this function cannot find a stable orientation, in
%  which case it throws a warning and returns an identity matrix R. The conf
%  score for these elements is negative.
%
%  Example:
%       im = rand(20, 20, 20);
%       keys = keypoint3D([10 10 10]);
%       [keys, conf] = orientation3D(keys, im);
%
%  See also:
%       keypoint3D, extractSift3D, imRead3D, setupSift3D
%
%  Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
function [R, conf] = orientation3D(keys, im, units)

% Verify inputs
if isempty(im)
    error('im is empty')
end

if ndims(im) ~= 3
    error(['im must have 3 dimensions, detected ' num2str(ndims(im))]);
end

if nargin < 3 || isempty(units)
        units = [];
end

units = checkUnits3D(units);

% Scale and convert the image to single precision
im = single(im);
im = im / max(im(:));

% Assign the orientations
[R, conf] = mexOrientation3D(keys, im, units);

% Check the validity
if any(conf) < 0
    warning('Failed to assign a stable orientation to some keypoints')
end

end
