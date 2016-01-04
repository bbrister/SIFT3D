%registerSift3D(src, ref, srcUnits, refUnits, thresh) register a pair of 
%  images using SIFT3D keypoints. Detects keypoints, extracts and matches
%  descriptors, and then fits an affine transformation using the RANSAC 
%  algorithm.
%
%  Arguments:
%    src, ref - A pair of [MxNxP] arrays, where voxels are indexed in 
%      (x, y, z) order. src is the "source" or "moving" image, and ref is 
%      the "reference" or "fixed" image.
%    srcUnits, refUnits - (Optional) The physical units for src and ref, 
%      respectively. See imRead3D for the format. Units should be provided
%      when attempting to register images of different resolutions, i.e.
%      5mm slices to 1mm slices.
%    thresh - (Optional) The matching threshold, in the interval (0, 1).
%      (Default: 0.8)
%
%  Return values:
%    A - A [4x3] matrix giving an affine transformation from the
%      coordinates of ref to those of src. The transformation can be
%      applied as follows:
%          [xt yt zt]' = A * [x y z 1]';
%    matchSrc, matchRef - A pair of [3xQ] arrays giving the matched
%      keypoint locations in src and ref, respectively. Each row is an
%      (x, y, z) coordinate, so that the ith row of matchSrc matches the
%      jth row of matchRef. These matches may contain outliers.
%
%  Examples:
%    % Load the images
%    [src, srcUnits] = imRead3D('src.dcm');
%    [ref, refUnits] = imRead3D('ref.dcm');
%
%    % Register with the default threshold
%    [A, matchSrc, matchRef] = registerSift3D(src, ref, srcUnits, ...
%                                             refUnits);
%    
%    % Register with a custom threshold
%    [A, matchSrc, matchRef] = registerSift3D(src, ref, srcUnits, ...
%                                             refUnits, 0.8);
%
%    % Register without units
%    [A, matchSrc, matchRef] = registerSift3D(src, ref);
%
%  See also:
%    imRead3D, imWrite3D, setupSift3D
%
%  Copyright (c) 2016 Blaine Rister et al., see LICENSE for details.
function [A, matchSrc, matchRef] = registerSift3D(src, ref, srcUnits, ...
                                                  refUnits, thresh)

    % Verify inputs
    narginchk(2, inf)
    if isempty(src)
        error('src is empty')
    end
    if isempty(ref)
        error('ref is empty')
    end
    
    if nargin < 3
        srcUnits = [];
    end
    if nargin < 4
        refUnits = [];
    end
    
    srcUnits = checkUnits3D(srcUnits, 'srcUnits');
    refUnits = checkUnits3D(refUnits, 'refUnits');
    
    if nargin < 5
        thresh = [];
    else
        if ~isnumeric(thresh)
            error('thresh must be numeric')
        end
        if thresh <= 0 || thresh >= 1
            error('thresh must be in the interval (0, 1)')
        end
    end
    
    % Convert the threshold to double precision
    thresh = double(thresh);
    
    % Convert the images to single precision and scale to [0, 1]
    src = imFormat(src);
    ref = imFormat(ref);
    
    % Register
    [A, matchSrc, matchRef] = mexRegisterSift3D(src, ref, srcUnits, ...
                                                refUnits, thresh);
end

% Helper function to convert images to single precision and scale to [0,1]
function im = imFormat(im)
    im = single(im);
    im = im / max(im(:));
end