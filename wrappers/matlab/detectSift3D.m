function keys = detectSift3D(im)
%detectSift3D(im) Detect Sift3D keypoints in an image.
%  Arguments:
%    im - An [MxNxP] array.
%
%  Return values:
%    keys - An array of n keypoint structs. Each struct has the following
%       fields:
%       key.coords - The [x y z] coordinates.
%       key.scale - The scale coordinate.
%       key.ori - A [3x3] rotation matrix representing the 3D orientation.
%
%  Example:
%     im = rand(50, 50, 50);
%     keys = detectSift3D(im);

% Verify inputs
if nargin < 1
        error('Not enough arguments');
end

if (ndims(im) ~= 3)
   error(['im must have 3 dimensions, detected ' num2str(ndims(im))]); 
end

% Scale and convert the image to single precision
im = single(im);
im = im / max(im(:));

% Detect features
[keys] = mexDetectSift3D(im);

end

