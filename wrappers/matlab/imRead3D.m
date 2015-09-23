function im = imRead3D(path)
%imRead3D(im) Read a 3D image from a file.
%  Arguments:
%    path - The path to the file.
%
%  Supported file formats:
%    NIFTI (.nii, .nii.gz)
%    DICOM (.dcm)
%    Directory of DICOM files (no extension)
%
%  Return values:
%    im - An [MxNxPxC] array containing the image data, scaled to the range
%    [0, 1]. The last dimension denotes the channels of the [MxNxP] image.
%    The voxels are indexed in (x, y, z, c) order, where c is the channel
%    index and (x, y, z) are the spatial coordinates.
%
%  Examples:
%      im = imRead3D('image.nii.gz'); % NIFTI
%      im = imRead3D('image.dcm'); % Multi-slice DICOM
%      im = imRead3D('image'); % Directory of DICOM slices
%
%  See also:
%    imWrite3D, detectSift3D, extractSift3D, setupSift3D
%
% Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

% Verify inputs
if nargin < 1 || isempty(path)
    error('Not enough arguments')
end

if ~exist(path, 'file')
    error('File does not exist')
end

% Read the image
im = mexImRead3D(path);

end

