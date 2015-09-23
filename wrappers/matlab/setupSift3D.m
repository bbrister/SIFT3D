% setupSift3D
%
% Run this script to set up the paths for the SIFT3D matlab toolbox. For 
% example, add the following line to your startup.m file:
% 
%   run('/path/to/SIFT3D/build/setupSift3D.m')
%
% where /path/to/SIFT3D/build is the path to the build directory of your
% SIFT3D installation.
%
% See also:
%   imRead3D, imWrite3D, detectSift3D, extractSift3D, keypoint3D
%
% Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

% Get the full path to this file
filename = mfilename('fullpath');

% Extract the toolbox directory
toolboxname = fileparts(filename);

% Add the toolbox to the path
addpath(genpath(toolboxname));

% Clean up the workspace
clear filename toolboxname