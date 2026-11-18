%% run_alignment_nematic.m  -  batch nematic order for TIFF images / stacks
% Edit the block below, then Run. Each *.tif in INPUT gets a
% <name>_alignment folder (order_t.mat, align_t.mat, alignmentParams.mat,
% nematic_order.csv, overlay PDF); a summary CSV (frame-1 value and mean
% over frames) is written to the input folder.
% Background-subtract images in ImageJ (rolling ball 50 px) beforehand,
% as in the Methods. Requires the Image Processing Toolbox (fspecial,
% improfile, bwlabel, regionprops) and nanmean (Statistics Toolbox).

clear; close all; clc;

%% ---------------- USER PARAMETERS ----------------
INPUT         = '';    % folder, single .tif, or cell list; '' = pick folder
P.winsize     = 19;    % window size (odd); Methods: 19 or 31
P.overlap     = 0.5;   % 50 % overlap
P.st          = 1;     % 1 -> 3x3 kernel of neighbouring windows
P.checkpoint  = 0;     % skip windows whose summed intensity <= this
                       % (Alba revision used 7220 for winsize 19, i.e. ~20 counts/px)
P.mask_method = 2;     % 1 = global disk, 2 = local threshold
P.frames      = [];    % [] = all frames
P.makeFigures = true;
%% --------------------------------------------------

files = list_tiffs(INPUT);
out = nan(numel(files), 3); names = cell(1, numel(files));
for k = 1:numel(files)
    [~, b, e] = fileparts(files{k}); names{k} = [b e];
    try
        S = nematic_order_stack(files{k}, P);
        out(k,:) = [S.nematic_order(1), mean(S.nematic_order), numel(S.frames)];
    catch err
        warning('%s failed: %s', names{k}, err.message);
    end
end
if ~isempty(files)
    write_csv(fullfile(fileparts(files{1}), 'nematic_order_summary.csv'), ...
        {'nematic_order_first_frame', 'nematic_order_mean_over_frames', 'n_frames'}, ...
        out, names, 'file');
end
