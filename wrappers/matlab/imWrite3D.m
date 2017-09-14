function imWrite3D(path, im, units)
%imWrite3D(im) Write a 3D image to a file.
%  Arguments:
%    path - The path to the file.
%    im - An [MxNxP] array containing the image data, where the voxels 
%       are indexed in (x, y, z) order.
%    units - (Optional) See imRead3D. Missing values default to 1.
%
%  Supported file formats:
%    NIFTI (.nii, .nii.gz)
%    Analyze (.img, .img.gz)   
%    DICOM (.dcm)
%    Directory of DICOM files (no extension)
%
%  Examples:
%    imWrite3D('image.nii.gz', im); % NIFTI
%    imWrite3D('image.dcm', im); % Multi-slice DICOM
%    imWrite3D('image', im); % Directory of DICOM slices
%    imWrite3D('image.dcm', im, [1 1 2]) % Anisotropic, multi-slice DICOM  
%
%  Notes:
%    When writing a DICOM file, the images values will be scaled and
%    rounded to 8 unsigned bits.
%
%  See also:
%    imRead3D, detectSift3D, extractSift3D, setupSift3D
%
% Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.

% Default parameters
if nargin < 3 || isempty(units)
   units = ones(3, 1); 
end

% Verify inputs
if nargin < 1 || isempty(path)
    error('path not specified')
end

if ~isa(path, 'char')
    error('path must be a string')
end

if nargin < 2 || isempty(im)
    error('im not specified')
end

if ndims(im) > 3
   error(['im has invalid dimensionality: ' num2str(ndims(im))]) 
end

units = checkUnits3D(units);

% Convert the image to single-precision
im = single(im);

% Check if this is a .gz file being written on Windows
winGz = false;
[pathstr, name, ext] = fileparts(path);
if strcmp(ext, '.gz') && ispc   
    % Strip .gz from the path
    winGz = true;
    path = fullfile(pathstr, name);
end

% Write the image
mexImWrite3D(path, im, units);

% On Windows, compress the file from within Matlab. Otherwise, it will
% crash.
if winGz
    gzip(path)
    delete(path)
end

end
