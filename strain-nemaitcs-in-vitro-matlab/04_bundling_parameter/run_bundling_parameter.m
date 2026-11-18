%% run_bundling_parameter.m  -  batch bundling parameter lambda = nu/mu
% Edit the block below, then Run. Background-subtract images in ImageJ
% (rolling ball 50 px) and crop to the analysis size (500 x 500 px in the
% Methods) beforehand. Writes <name>_bundling.csv per image and
% bundling_summary.csv in the input folder.

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
INPUT              = '';                % folder, single .tif, or cell list; '' = pick folder
P.method           = 'normalize_total'; % 'normalize_total' (Methods) | 'subtract_mean'
P.subtractFraction = 0.9;               % for 'subtract_mean'
P.statMode         = 'original';        % 'original' (nanstd(nanstd)) | 'global'
P.frames           = [];                % [] = all frames
%% --------------------------------------------------

files = list_tiffs(INPUT);
out = nan(numel(files), 4); names = cell(1, numel(files));
for k = 1:numel(files)
    [~, b, e] = fileparts(files{k}); names{k} = [b e];
    try
        S = bundling_parameter(files{k}, P);
        out(k,:) = [S.mu(1), S.nu(1), S.lambda(1), mean(S.lambda)];
    catch err
        warning('%s failed: %s', names{k}, err.message);
    end
end
if ~isempty(files)
    write_csv(fullfile(fileparts(files{1}), 'bundling_summary.csv'), ...
        {'mu_frame1', 'nu_frame1', 'lambda_frame1', 'lambda_mean_over_frames'}, out, names, 'file');
end
