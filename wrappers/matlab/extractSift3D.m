function [desc, coords] = extractSift3D(keys, im, units)
%extractSift3D(keys, im, units) Extract Sift3D descriptors from keypoints.
%  Arguments:
%    keys - An array of n keypoint structs. See keypoint3D.m for the
%      format.
%    im - (Optional) An [MxNxP] array, where voxels are indexed in
%      (x, y, z) order. If im is empty or not provided, this function uses 
%      the Gaussian scale-space pyramid from the most recent call to 
%      detectSift3D. (Note: this data is overwritten by calls to
%      registerSift3D.)
%    units - (Optional) See imRead3D. If units are specified, the extracted
%       descriptors are isotropic even when im is not.
%
%  Return values:
%    desc - An [n x 768] array of descriptors. The ith row is a descriptor
%      corresponding to keys(i).
%    coords - An [n x 3] array of descriptor coordinates, defined in the
%      space of the input image (see the description of the "im" argument).
%
%  Examples:
%    % Extract without units
%    im = rand(50, 50, 50);
%    keys = detectSift3D(im);
%    [desc, coords] = extractSift3D(keys);
%
%    % Extract with units
%    [im, units] = imRead3D('someFile.dcm');
%    keys = detectSift3D(im, units);
%    [desc, coords] = extractSift3D(keys);
%
%    % Extract from manually-defined keypoints
%    [im, units] = imRead3D('someFile.dcm');
%    keys = keypoint3D([0 0 0]);
%    keys = orientation3D(keys, im, units);
%    [desc, coords] = extractSift3D(keys, im, units);
%
%  See also:
%    detectSift3D, imRead3D, keypoint3D, orientation3D, setupSift3D
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Required field dimensions
coordsSizeReq = [1 3];
oriSizeReq = [3 3];

% The number of descriptor elements
descNumel = 768;

% Default parameters
if nargin < 2
    im = [];
end

if nargin < 3
    units = [];
end

% Verify inputs
if nargin < 1
    error('Not enough arguments')
end

if ~isa(keys, 'struct')
   error('keys must be a struct array') 
end

% Do nothing if we have no keypoints
if isempty(keys)
   warning('keys is empty')
   desc = [];
   coords = [];
   return
end

coordsSize = size(keys(1).coords);
if any(coordsSize ~= coordsSizeReq)
    error(['keys.coords must have size [' num2str(coordsSizeReq) '].' ...
        'Detected size: [' num2str(coordsSize) ']']);
end

if ~isa(keys(1).coords, 'double')
    error('keys.coords must have type double');
end

if ~isscalar(keys(1).scale)
   error('keys.scale must be a scalar'); 
end

if ~isa(keys(1).scale, 'double')
    error('keys.scale must have type double');
end

oriSize = size(keys(1).ori);
if any(oriSize ~= oriSizeReq)
    error(['keys.ori must have size [' num2str(oriSizeReq) '].' ...
        'Detected size: [' num2str(oriSize) ']']);
end

if ~isa(keys(1).ori, 'double')
    error('keys.ori must have type double');
end

if nargin < 2
    im = [];
end

% Convert to single precision
im = single(im);

% Verify and scale the input image, if any
if (~isempty(im))
    if (ndims(im) ~= 3)
        error(['im must have 3 dimensions, detected ' num2str(ndims(im))]);
    end
    
    im = im / (max(im(:)) + eps);
end

units = checkUnits3D(units);

% Detect features
ret = mexExtractSift3D(keys, im, units);

% Splice the outputs
coords = ret(:, 1 : 3);
desc = ret(:, 4 : end);

assert(size(desc, 2) == descNumel);

end

