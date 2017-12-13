function matches = matchSift3D(desc1, coords1, desc2, coords2, nnThresh)
%matchSift3D(desc1, desc2)
% Match SIFT3D desriptors from a pair of images. 
%
%  Arguments:
%    desc1: Descriptors returned from extractSift3D, for the first image.
%    coords1: Coordinates returned from extractSift3D, for the second 
%      image. (Default: zeros)
%    desc2: Like desc1, but for the second image.
%    coords2: Like coords1, but for the second image. (Default: zeros)
%    nnThresh: The matching threshold, in the interval (0, 1]. A higher 
%      value means fewer matches will be returned. (Default: 0.8) 
%
%  Return values:
%    matches - An [m x 2] array, where matches(m, 1) indexes a row of
%      desc1, and matches(m, 2) indexes a row of desc2. The descriptors
%      given by these
%
%  If either of coords1, coords2 are empty, default values will be used.
%
%  Example:
%    % Get descriptors from the first image
%    [im1, units1] = imRead3D('someFile.dcm');
%    keys1 = detectSift3D(im1, units1);
%    [desc1, coords1] = extractSift3D(keys1);
%    
%    % Get descriptors from the second image
%    [im2, units2] = imRead3D('someOtherFile.dcm');
%    keys2 = detectSift3D(im2, units2);
%    [desc2, coords2] = extractSift3D(keys2);   
%
%    % Match the descriptors
%    matches = matchSift3D(desc1, coords1, desc2, coords2, 0.8);
%
%  See also:
%    imRead3D, imWrite3D, registerSift3D, setupSift3D

% Supply defaults
if isempty(coords1)
    coords1 = zeros(size(desc1, 1), 3);
end
if nargin < 4 || isempty(coords2)
    coords2 = zeros(size(desc2, 1), 3);
end
if nargin < 5
   nnThresh = []; 
end

% Verify inputs
narginchk(3, 4)
validateDesc('desc1', desc1, 'coords1', coords1)
validateDesc('desc2', desc2, 'coords2', coords2)
if ~isempty(nnThresh)
    validateattributes(nnThresh, {'numeric'}, ...
        {'real', 'scalar', '>=', 0, '<', 1}, 'nnThresh')
end

% Format the descriptors as a C-style matrix
desc1f = formatDesc(desc1, coords1);
desc2f = formatDesc(desc2, coords2);

% Match the descriptors
matchIdx = mexMatchSift3D(desc1f, desc2f, nnThresh);

% Convert the matches to indices in each array
validMatches = matchIdx >= 0;
matches = [find(validMatches), matchIdx(validMatches) + 1];

end

function validateDesc(descName, desc, coordsName, coords)
% Verify the descriptors and coordinates, throwing errors if invalid.

validateattributes(desc, {'numeric'}, {'real', '2d', 'ncols', 768}, ...
    descName)
validateattributes(coords, {'numeric'}, {'real', '2d', 'ncols', 3}, ...
    coordsName)

if size(desc, 1) ~= size(coords, 1)
    error([descName ' has size ' size(desc) ', but ' coordsName ' has ' ...
        'size ' size(coords) '. Number of rows must be equal.'])
end

end

function ret = formatDesc(desc, coords)
% Combine the coordinates and descriptors and convert to double
    ret = double([coords, desc]);
end
