% setup_paths_RFs.m
%
% adds necessary paths for running do_RFs.m on vader


function setup_paths_RFs

% if nargin < 1
% add gridfitgpu
addpath(genpath('~/Documents/MATLAB/toolboxes/gridfitgpu/'));

% add local version of vistasoft w/ gpu functionality
addpath(genpath('~/Documents/MATLAB/toolboxes_dev/vistasoft_ts/'));

% add this directory (where helper functions, etc, live) |||| NOT needed,
% put this in startup.m and it works ok; users will just need to use a
% startup that points to necessary preprocessing functions
%fprintf('%s\n',mfilename('fullpath'))
%addpath(mfilename('fullpath'));


return