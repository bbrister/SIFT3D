function keys = detectSift3D(im, varargin)
%detectSift3D(im, options) Detect Sift3D keypoints in an image.
%  Arguments:
%    im - An [MxNxP] array, where voxels are indexed in (x, y, z) order.
%
%  Options:
%    units - See imRead3D. If units are specified, the detected
%       keypoints are isotropic even when im is not.
%    firstOctave - The first octave a pyramid. Must be a non-negative 
%       integer. (default: 0)
%    peakThresh - The smallest allowed absolute DoG value, as a fraction
%       of the largest. Must be in the interval (0, 1]. (default: 0.10)
%    cornerThresh - The smalled allowed corner score, on the interval
%       [0, 1]. (default: 0.5)
%    numOctaves - The number of octaves to process. Must be an integer.
%       (default: process as many as we can)
%    numKpLevels - The number of pyramid levels per octave in which
%       keypoints are found. Must be a positive integer. (default: 3)
%    sigmaN - The nominal scale parameter of the input data, on the
%       interval (0, inf). (default: 1.15)
%    sigma0 - The scale parameter of the first level of octave 0, in the
%       interval (0, inf). (default: 1.6)
%
%  Return values:
%    keys - A [Qx1] array of keypoint structs. See keypoint3D.m for the
%      struct definition.
%
%  Note: This function will reserve memory that persists after it has
%  finished. To release all SIFT3D memory, use 'clear mex'.
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
%      keys = detectSift3D(im, 'units', units);
%
%  See also:
%    extractSift3D, imRead3D, imWrite3D, keypoint3D, setupSift3D
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Option names
unitsStr = 'units';
firstOctaveStr = 'firstOctave';
peakThreshStr = 'peakThresh';
cornerThreshStr = 'cornerThresh';
numOctavesStr = 'numOctaves';
numKpLevelsStr = 'numKpLevels';
sigmaNStr = 'sigmaN';
sigma0Str = 'sigma0';

% Parse options
parser = inputParser;
parser.addParamValue(unitsStr, [])
parser.addParamValue(firstOctaveStr, [])
parser.addParamValue(peakThreshStr, [])
parser.addParamValue(cornerThreshStr, [])
parser.addParamValue(numOctavesStr, [])
parser.addParamValue(numKpLevelsStr, [])
parser.addParamValue(sigmaNStr, [])
parser.addParamValue(sigma0Str, [])
parse(parser, varargin{:})
units = parser.Results.units;
firstOctave = parser.Results.firstOctave;
peakThresh = parser.Results.peakThresh;
cornerThresh = parser.Results.cornerThresh;
numOctaves = parser.Results.numOctaves;
numKpLevels = parser.Results.numKpLevels;
sigmaN = parser.Results.sigmaN;
sigma0 = parser.Results.sigma0;

% Verify inputs
narginchk(1, inf)
if isempty(im)
    error('im is empty')
end
if ndims(im) ~= 3
    error(['im must have 3 dimensions, detected ' num2str(ndims(im))]);
end
units = checkUnits3D(units);
if ~isempty(firstOctave)
    validateattributes(firstOctave, {'numeric'}, ...
        {'real', 'integer', 'scalar', 'nonnegative'}, 'firstOctave')
end
if ~isempty(peakThresh)
    validateattributes(peakThresh, {'numeric'}, ...
        {'real', 'positive', 'scalar'}, 'peakThresh')
end
if ~isempty(cornerThresh)
    validateattributes(cornerThresh, {'numeric'}, ...
        {'real', 'nonnegative', 'scalar', '<=', 1}, 'cornerThresh')
end
if ~isempty(numOctaves)
    validateattributes(numOctaves, {'numeric'}, ...
        {'real', 'integer', 'scalar', 'positive'}, 'numOctaves')
end
if ~isempty(numKpLevels)
    validateattributes(numKpLevels, {'numeric'}, ...
        {'real', 'integer', 'scalar', 'positive'}, 'numKpLevels')
end
if ~isempty(sigmaN)
    validateattributes(sigmaN, {'numeric'}, ...
        {'real', 'positive', 'scalar'}, 'sigmaN')
end
if ~isempty(sigma0)
    validateattributes(sigma0, {'numeric'}, ...
        {'real', 'positive', 'scalar'}, 'sigma0')
end
if sigmaN >= sigma0
    error('Cannot have sigmaN >= sigma0')
end

% Scale and convert the image to single precision
im = single(im);
im = im / (max(im(:)) + eps);

% Collect the options in a struct
optStruct = struct(firstOctaveStr, firstOctave, ...
    peakThreshStr, peakThresh, ...
    cornerThreshStr, cornerThresh, ...
    numOctavesStr, numOctaves, ...
    numKpLevelsStr, numKpLevels, ...
    sigmaNStr, sigmaN, ...
    sigma0Str, sigma0);

% Detect features
keys = mexDetectSift3D(im, units, optStruct);

end

