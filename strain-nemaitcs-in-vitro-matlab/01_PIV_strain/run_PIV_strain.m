%% run_PIV_strain.m  -  batch PIV + strain / strain rate for TIFF stacks
% Edit the block below, then Run. Every *.tif in INPUT is analysed; each
% stack gets its own <name>_PIV folder and one summary CSV is written to
% the input folder. Requires: setup_paths.m run once per session.
%
% Method (Methods text): PIV on fluorescent actin (MatPIV), window 32 px,
% overlap 0.5; strain eps = <div x> of the cumulative displacement;
% background-subtract stacks in ImageJ (rolling ball 50 px) beforehand.

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
INPUT     = '';          % folder of stacks, one .tif, or {'a.tif','b.tif'}; '' = pick folder
P.dscl    = 0.1076;      % um / pixel       (Alba data: 0.1076)
P.tscl    = 10.5;        % s / frame        (Alba data: 10.5)
P.winsize = 32;          % PIV window (px)
P.overlap = 0.5;         % PIV overlap
P.lastim  = 0;           % last frame to analyse (0 = all)
P.rateMethod  = 'pct25_75';  % 'pct25_75' | 'maxslope' | 'manual'
P.fitWindow_s = 15;          % for 'maxslope'
P.fitRange_s  = [];          % for 'manual', e.g. [100 250]
P.makeFigures = true;
%% --------------------------------------------------

files = list_tiffs(INPUT);
summary = nan(numel(files), 6); names = cell(1, numel(files));
for k = 1:numel(files)
    [~, b, e] = fileparts(files{k}); names{k} = [b e];
    try
        S = piv_strain_analysis(files{k}, P);
        summary(k,:) = [S.n_frames, S.max_strain, S.t_max_strain_s, ...
                        S.strain_rate_per_s, S.fit_start_s, S.fit_end_s];
    catch err
        warning('%s failed: %s', names{k}, err.message);
    end
end
if ~isempty(files)
    write_csv(fullfile(fileparts(files{1}), 'PIV_strain_summary.csv'), ...
        {'n_frames','max_strain','t_max_strain_s','strain_rate_per_s','fit_start_s','fit_end_s'}, ...
        summary, names, 'file');
end
