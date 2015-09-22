% ioMatlab
%
% Example of using file IO functions with SIFT3D
%
% Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

% Load the from a NIFTI file
imNifti = imRead3D('data/1.nii.gz');

% Save it as a multi-slice DICOM file
imWrite3D('1.dcm', imNifti);