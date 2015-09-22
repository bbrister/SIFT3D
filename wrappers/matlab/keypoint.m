function keys = keypoint(coords, scale, ori)
%keys = keypoint(coords, scales, orientations)
% Create an array of Keypoint structs.
%
% Inputs:
%   -coords: [mx3] matrix of locations, where each row is an [x, y, z]
%   coordinate triple, 0-indexed (required)
%   -scale: [mx1] vector of nonnegative scale parameters (default: 1.6)
%   -ori [3x3xm] array, where ori(:, :, i) is the [3x3] rotation matrix 
%   of the ith keypoint (default: [3x3] identity matrix)
%
% Leaving an argument empty results in setting that value to the default.
%
% Example:
%   coords = [1 1 1; 2 2 2];
%   scale = [1 2];
%   ori = repmat(eye(3), [1 1 2]);
%   keys = keypoint(coords, scale, ori);

% Default parameters
n = 3;
scaleDefault = 1.6;
oriDefault = eye(n);

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
elseif size(scale) == [1 m]
        scale = scale';
elseif size(scale) ~= [m 1]
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
        if (abs(d - 1) > 0.1)
            error(['det(ori(:, :, ' num2str(i) ')) = ' num2str(d) ...
              ', must be equal or near to 1'])
        end
        
        % Verify orthogonality
        RRt = R * R';
        if (abs(RRt - eye(size(R))) > 0.1)
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