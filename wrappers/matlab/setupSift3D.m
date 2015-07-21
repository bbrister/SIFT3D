% Get the full path to this file
filename = mfilename('fullpath');

% Extract the toolbox directory
toolboxname = fileparts(filename);

% Add the toolbox to the path
addpath(genpath(toolboxname));

% Clean up the workspace
clear filename toolboxname