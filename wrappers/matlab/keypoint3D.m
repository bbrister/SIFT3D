function keys = keypoint(coords, scale, ori)
%keypoint3D(coords, scales, orientations) Create an array of Keypoint 
% structs.
%
% Inputs:
%   coords - [Mx3] matrix of locations, where each row is an [x, y, z]
%     coordinate triple, 0-indexed (required)
%   scale - [Mx1] vector of nonnegative scale parameters (default: 1.6)
%   ori - [3x3xM] array, where ori(:, :, i) is the [3x3] rotation matrix 
%   of the ith keypoint (default: identity matrix)
%
% Leaving an argument empty results in setting that value to the default.
%
% Return values:
%   keys - an [Mx1] array of keypoint structs. Each struct has the following
%      fields:
%      key.coords - The [x y z] coordinates, 0-indexed.
%      key.scale - The scale coordinate.
%      key.ori - A [3x3] rotation matrix representing the 3D orientation.
%      key.octave - The pyramid octave index.
%      key.level - The pyramid level index within that octave.
%
% Example:
%   coords = [1 1 1; 2 2 2];
%   scale = [1 2];
%   ori = repmat(eye(3), [1 1 2]);
%   keys = keypoint3D(coords, scale, ori);
%
% See also:
%   orientation3D, extractSift3D, setupSift3D
%
% Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.

% Default parameters
n = 3;
scaleDefault = 1.6;
oriDefault = eye(n);

% Error tolerance
tol = 1E3 * eps;

% Verify coords
if nargin < 1 || isempty(coords)
    error('Coords argument must be specified')
elseif ~isa(coords, 'double')
    error('Coords argument must be of type double')
elseif ~isreal(coords)
    error('Coords argument must be real-valued')
end

if size(coords, 2) ~= n || ~ismatrix(coords)
   error(['coords argument has invalid dimensions [' ...
        num2str(size(coords)) '] must be [mx3]'])
end
m = size(coords, 1);

% Verify scale
if nargin < 2 || isempty(scale)
    scale = repmat(scaleDefault, [m 1]);
elseif isequal(size(scale), [1 m])
        scale = scale';
elseif ~isequal(size(scale), [m 1])
    error(['scale argument has invalid dimensions [' ...
        num2str(size(scale)) '] must be [mx1]']);
elseif ~isa(scale, 'double')
    error('scale argument must be of type double')
elseif ~isreal(scale)
    error('scale argument must be real-valued')
elseif any(scale <= 0)
    error('scale argument must be positive')
end

% Verify orientation
if nargin < 3 || isempty(ori)
    ori = repmat(oriDefault, [1 1 m]);
elseif ndims(ori) == 2 && ~isequal(size(ori), [n n]) || ...
    ndims(ori) == 3 && ~isequal(size(ori), [n n m])
    nStr = num2str(n);
    error(['ori argument has invalid dimensions [' num2str(size(ori)) ...
        '], must be [' nStr 'x' nStr 'xm]'])
elseif ~isa(ori, 'double')
    error('scale argument must be of type double')
elseif ~isreal(scale)
    error('scale argument must be real-valued')
else
    % Verify rotation matrices
    for i = 1 : m
        
        R = ori(:, :, i);
        
        % Verify determinant
        d = det(R);
        if (abs(d - 1) > tol)
            error(['det(ori(:, :, ' num2str(i) ')) = ' num2str(d) ...
              ', must be equal or near to 1'])
        end
        
        % Verify orthogonality
        RRt = R * R';
        if (abs(RRt - eye(size(R))) > tol)
            error(['ori(:, :, ' num2str(i) ')) must be orthogonal'])
        end
    end 
end

% Create default octave and level arrays
octave = 0;
level = 0;

% Create the struct
keys = struct('coords', num2cell(coords, 2), 'scale', ...
    num2cell(scale, 2), 'ori', reshape(num2cell(ori, [1 2]), [m 1]), ...
    'octave', num2cell(octave, 2), 'level', num2cell(level, 2));

end
