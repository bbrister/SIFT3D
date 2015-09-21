function [desc, coords] = extractSift3D(keys, im)
%extractSift3D(keys, im) Extract Sift3D descriptors from keypoints.
%  Arguments:
%    keys - An array of n keypoint structs. See detectSift3D for the
%      format.
%    im - (Optional) An [MxNxP] array, where voxels are indexed in
%      (x, y, z) order. If im is empty or not provided, this function uses 
%      the Gaussian scale-space pyramid from the most recent call to 
%      detectSift3D.
%
%  Return values:
%    desc - An [n x 768] array of descriptors. The ith row is a descriptor
%      corresponding to keys(i).
%    coords - An [n x 3] array of descriptor coordinates, defined in the
%      space of the input image (see the description of the "im" argument).
%
%  Example:
%    im = rand(50, 50, 50);
%    keys = detectSift3D(im);
%    [desc, coords] = extractSift3D(keys);
%
%  See also:
%    detectSift3D, imRead3D, setupSift3D
%
% Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

% Required field dimensions
coordsSizeReq = [1 3];
oriSizeReq = [3 3];

% The number of descriptor elements
descNumel = 768;

% Verify inputs
if nargin < 1 || isempty(keys)
    error('Not enough arguments');
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
    
    im = im / max(im(:));
end


% Detect features
ret = mexExtractSift3D(keys, im);

% Splice the outputs
coords = ret(:, 1 : 3);
desc = ret(:, 4 : end);

assert(size(desc, 2) == descNumel);

end

