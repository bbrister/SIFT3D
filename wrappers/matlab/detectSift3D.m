function keys = detectSift3D(im, sift3d)
%detectSift3D Detect Sift3D keypoints in an image.
%   Detailed explanation goes here

% Verify inputs
if (ndims(im) != 3)
   error(['im must have 3 dimensions, detected ' num2str(ndims(im))]); 
end

if (min(size(im)) < 1)
   error(['Invalid image dimensions: ' num2str(size(im))]); 
end

% Scale and convert the image to single precision
im = single(im);
im = im / max(im(:));

% Detect features
[keys] = mexDetectSift3D(im);

end

