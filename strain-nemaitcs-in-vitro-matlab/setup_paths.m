function setup_paths()
% SETUP_PATHS  Add every analysis folder of this package to the MATLAB path.
%   Run once per MATLAB session (or put it in startup.m):
%       >> run('/path/to/Alba_MATLAB_analysis_code/setup_paths.m')
%   The 'original_scripts' sub-folders are NOT added: they hold the
%   untouched historical scripts for reference only.
root = fileparts(mfilename('fullpath'));
p = strsplit(genpath(root), pathsep);
keep = p(~cellfun(@isempty, p));
keep = keep(cellfun(@isempty, strfind(keep, 'original_scripts')));
keep = keep(cellfun(@isempty, strfind(keep, [filesep '.'])));
addpath(strjoin(keep, pathsep));
fprintf('Added %d folders under %s to the path.\n', numel(keep), root);
end
