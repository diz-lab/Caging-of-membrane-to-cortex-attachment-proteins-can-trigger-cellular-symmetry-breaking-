%% run_pair_correlation_gofr.m  -  batch g(r) on binarised actin images
% Edit the block below, then Run. Input images must be binarised (e.g. in
% ImageJ) or a grey-level threshold must be given. Writes <name>_gofr.*
% per image and gofr_summary.csv in the input folder.

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
INPUT       = '';          % folder, single .tif, or cell list; '' = pick folder
P.pixPerUm  = 9.2937;      % pixels per um (= 1/0.1076 um/px)
P.frame     = 1;           % frame to analyse
P.rlim      = 35;          % max r (px)
P.dr        = 1;           % bin width (px)
P.threshold = [];          % [] = binary input; else foreground = I >= threshold
P.edgeMode  = 'original';  % 'original' (reproduces earlier g(r)) | 'interior'
P.fitBins   = [2 14];      % bins used in the exponential fit
P.makeFigures = true;
%% --------------------------------------------------

files = list_tiffs(INPUT);
out = nan(numel(files), 5); names = cell(1, numel(files));
for k = 1:numel(files)
    [~, b, e] = fileparts(files{k}); names{k} = [b e];
    try
        S = pair_correlation_gofr(files{k}, P);
        out(k,:) = [S.a(1), S.xi_um, S.a(3), S.n_particles, S.n_interior];
    catch err
        warning('%s failed: %s', names{k}, err.message);
    end
end
if ~isempty(files)
    write_csv(fullfile(fileparts(files{1}), 'gofr_summary.csv'), ...
        {'A', 'xi_um', 'c', 'n_particles', 'n_interior'}, out, names, 'file');
end
