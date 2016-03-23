function [A, matchSrc, matchRef] = registerSift3D(src, ref, varargin)
%registerSift3D(src, ref, options) register a pair of
%  images using SIFT3D keypoints. Detects keypoints, extracts and matches
%  descriptors, and then fits an affine transformation using the RANSAC
%  algorithm.
%
%  Arguments:
%    src, ref - A pair of [MxNxP] arrays, where voxels are indexed in
%      (x, y, z) order. src is the "source" or "moving" image, and ref is
%      the "reference" or "fixed" image.
%
%  Options:
%    srcUnits, refUnits - The physical units for src and ref,
%      respectively. See imRead3D for the format. Units should be provided
%      when attempting to register images of different resolutions, i.e.
%      5mm slices to 1mm slices. (Default: [1 1 1])
%    nnThresh - The matching threshold, in the interval (0, 1].
%       (Default: 0.8)
%    errThresh - The RANSAC inlier threshold, in the interval (0, inf).
%       This is a threshold on the squared Euclidean distance in real-world
%       units. (Default: 5.0)
%    numIter - The number of RANSAC iterations. (Default: 500)
%    resample - If given, resamples src and ref to have the same resolution
%       prior to registration. This is slow. Use it when the inputs have
%       vastly different units.
%
%  Return values:
%    A - A [4x3] matrix giving an affine transformation from the
%      coordinates of ref to those of src. The transformation can be
%      applied as follows:
%          [xt yt zt]' = A * [x y z 1]';
%    matchSrc, matchRef - A pair of [3xQ] arrays giving the matched
%      keypoint locations in src and ref, respectively. Each row is an
%      (x, y, z) coordinate, so that the ith row of matchSrc matches the
%      ith row of matchRef. These matches may contain outliers.
%
%  Examples:
%    % Load the images
%    [src, srcUnits] = imRead3D('src.dcm');
%    [ref, refUnits] = imRead3D('ref.dcm');
%
%    % Register with units
%    [A, matchSrc, matchRef] = registerSift3D(src, ref, 'srcUnits', ...
%       [1 1 1], 'refUnits', [1 1 5], 'resample');
%
%    % Register without units
%    [A, matchSrc, matchRef] = registerSift3D(src, ref, 'resample');
%
%  See also:
%    imRead3D, imWrite3D, setupSift3D
%
%  Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Option names
srcUnitsStr = 'srcUnits';
refUnitsStr = 'refUnits';
numIterStr = 'numIter';
nnThreshStr = 'nnThresh';
errThreshStr = 'errThresh';
resampleStr = 'resample';

% Parse options
parser = inputParser;
parser.addParamValue(srcUnitsStr, [])
parser.addParamValue(refUnitsStr, [])
parser.addParamValue(numIterStr, [])
parser.addParamValue(nnThreshStr, [])
parser.addParamValue(errThreshStr, [])
parser.addOptional(resampleStr, false)
parse(parser, varargin{:})
srcUnits = parser.Results.srcUnits;
refUnits = parser.Results.refUnits;
numIter = parser.Results.numIter;
nnThresh = parser.Results.nnThresh;
errThresh = parser.Results.errThresh;
resample = parser.Results.resample;

% Verify inputs
narginchk(2, inf)
if isempty(src)
    error('src is empty')
end
if isempty(ref)
    error('ref is empty')
end
srcUnits = checkUnits3D(srcUnits, 'srcUnits');
refUnits = checkUnits3D(refUnits, 'refUnits');
validateattributes(numIter, {'numeric'}, {'real', '>', 0}, 'numIter')
validateattributes(nnThresh, {'numeric'}, {'real', '>=', 0, '<', 1}, ...
    'nnThresh')
validateattributes(errThresh, {'numeric'}, {'real', '>', 0}, 'errThresh')
validateattributes(resample, {'logical'}, {'scalar'}, 'resample')

% Convert the images to single precision and scale to [0, 1]
src = imFormat(src);
ref = imFormat(ref);

% Collect the options in a struct
optStruct = struct(numIterStr, numIter, ...
    nnThreshStr, nnThresh, ...
    errThreshStr, errThresh);

% Register
[A, matchSrc, matchRef] = mexRegisterSift3D(src, ref, srcUnits, ...
    refUnits, resample, optStruct);
end

% Helper function to convert images to single precision and scale to [0,1]
function im = imFormat(im)
im = single(im);
im = im / (max(im(:)) + eps);
end
