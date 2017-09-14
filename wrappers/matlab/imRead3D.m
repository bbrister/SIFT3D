function [im, units] = imRead3D(path)
%imRead3D(im) Read a 3D image from a file.
%  Arguments:
%    path - The path to the file.
%
%  Supported file formats:
%    NIFTI (.nii, .nii.gz)
%    Analyze (.img, .img.gz)   
%    DICOM or DICOM Segmentation Object (.dcm)
%    Directory of DICOM files (no extension)
%
%  Return values:
%    im - An [MxNxPxC] array containing the image data. The last dimension 
%       denotes the channels of the [MxNxP] image. The voxels are indexed 
%       in (x, y, z, c) order, where c is the channel index and (x, y, z) 
%       are the spatial coordinates.
%    units - A [3x1] vector of the (x, y, z) real-world units for the
%       voxels in im.
%
%  Examples:
%      [im, units] = imRead3D('image.nii.gz'); % NIFTI
%      [im, units] = imRead3D('image.dcm'); % Multi-slice DICOM or DSO
%      [im, units] = imRead3D('image'); % Directory of DICOM slices
%
%  Notes:
%   When possible, the image values will be in a vendor-neutral format. For 
%   example, a DICOM CT scan will be converted to Hounsfield units.
%
%   When reading a DICOM Segmentation Object (DSO), the program assumes
%   that the DICOM images to which the DSO refers are located in the same
%   directory as the DSO.
%
%  See also:
%    imWrite3D, detectSift3D, extractSift3D, setupSift3D
%
% Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.

% Verify inputs
if nargin < 1 || isempty(path)
    error('Not enough arguments')
end

if ~exist(path, 'file')
    error('File does not exist')
end

% Read the image
[im, units] = mexImRead3D(path);

end

