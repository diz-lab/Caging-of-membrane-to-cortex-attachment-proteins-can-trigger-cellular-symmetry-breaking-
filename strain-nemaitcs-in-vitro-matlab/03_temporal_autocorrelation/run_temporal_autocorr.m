%% run_temporal_autocorr.m  -  batch temporal autocorrelation of TIFF stacks
% Edit the block below, then Run. Each stack gives <name>_autocorr.mat/.csv/.pdf
% next to it, and temporal_autocorr_summary.csv is written to the input folder.
% Background-subtract stacks in ImageJ beforehand (Methods).

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
INPUT         = '';     % folder, single .tif, or cell list; '' = pick folder
P.dt          = 10.5;   % s / frame  (Alba data: 10.5)
P.frame1      = 1;      % first frame
P.frameEnd    = 300;    % last frame (0 = last in file)
P.fitFraction = 0.3;    % fit over this fraction of lags (original code: 0.3)
P.a0          = [0.9 -2000 -50];  % initial guess [A, -tau1, -tau2]
P.fitOffset   = false;  % true = add constant a4 (Methods equation)
P.makeFigures = true;
%% --------------------------------------------------

files = list_tiffs(INPUT);
out = nan(numel(files), 7); names = cell(1, numel(files));
for k = 1:numel(files)
    [~, b, e] = fileparts(files{k}); names{k} = [b e];
    try
        S = temporal_autocorr(files{k}, P);
        out(k,:) = [S.A, S.tau1_s, S.tau2_s, S.tau_long_s, S.offset, S.R2, S.n_frames];
    catch err
        warning('%s failed: %s', names{k}, err.message);
    end
end
if ~isempty(files)
    write_csv(fullfile(fileparts(files{1}), 'temporal_autocorr_summary.csv'), ...
        {'A', 'tau1_s', 'tau2_s', 'tau_long_s', 'offset', 'R2', 'n_frames'}, out, names, 'file');
end
