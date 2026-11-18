function S = nematic_order_stack(tifPath, P)
% NEMATIC_ORDER_STACK  F-actin orientation field + local nematic order, every frame.
%
%   S = nematic_order_stack('/data/actin.tif', P)
%
% Engine: FFTAlignment.m (Oakes / Gardel lab, updated by Linsmeier), called
% exactly as alignment.m (Banerjee) calls it. For each window
% (winsize x winsize, given overlap) the Gaussian-filtered window is
% Fourier transformed; the orientation is orthogonal to the axis of least
% second moment of the masked |FFT|. The local order of each window is
%   n = 2*( <cos^2(theta_ij - theta_i)> - 1/2 )
% over the (2*st+1)^2 kernel of neighbouring windows (st = 1 -> 3x3).
% Per-frame summary (as geT_Nematic_ZGS.m):
%   nematic_order = abs(nanmean(nanmean(order)))  (mean of column means)
% plus the plain mean over all valid windows for comparison.
%
% P fields (optional):
%   winsize     window size, odd, px              [19]
%   overlap     window overlap (0..1)             [0.5]
%   st          kernel half-width (windows)       [1]  -> 3x3 kernel
%   checkpoint  min. summed intensity in a window for it to be analysed
%               (absolute, in image counts summed over winsize^2 px)  [0]
%   mask_method 1 = global disk mask, 2 = local threshold mask [2]
%   frames      frames to analyse              [all]
%   outdir      [<stack folder>/<stack name>_alignment]
%   makeFigures overlay of directors on first analysed frame [true]
%
% Values used for the Alba revision (alignmentParams.mat, Aug 2025):
%   winsize 19, overlap 0.5, st 1, checkpoint 7220, mask_method 2.
%
% Outputs in outdir (same file format as alignment.m):
%   order_<t>.mat ('order'), align_<t>.mat ('pos_vec' = [col row dcol drow]),
%   alignmentParams.mat ('options'), nematic_order.csv, nematicAlignment.pdf

D = struct('winsize', 19, 'overlap', 0.5, 'st', 1, 'checkpoint', 0, ...
           'mask_method', 2, 'frames', [], 'outdir', '', 'makeFigures', true);
P = set_defaults(P, D);
if mod(P.winsize, 2) ~= 1, error('nematic_order_stack:winsize', 'winsize must be odd.'); end
[folder, base, ext] = fileparts(tifPath);
if isempty(P.outdir), P.outdir = fullfile(folder, [base '_alignment']); end
if ~exist(P.outdir, 'dir'), mkdir(P.outdir); end

nAll = numel(imfinfo(tifPath));
frames = P.frames; if isempty(frames), frames = 1:nAll; end
nemColMean = nan(numel(frames), 1); nemAll = nan(numel(frames), 1); nValid = nan(numel(frames), 1);
for k = 1:numel(frames)
    t = frames(k);
    im = imread(tifPath, t);
    if size(im, 3) > 1, im = im(:,:,1); end
    Pt = FFTAlignment(im, P.winsize, P.overlap, P.st, P.checkpoint, P.mask_method, 0);
    pos_vec = [Pt.pos, Pt.vec]; %#ok<NASGU>
    order = Pt.ordermat;
    save(fullfile(P.outdir, ['order_' num2str(t) '.mat']), 'order');
    save(fullfile(P.outdir, ['align_' num2str(t) '.mat']), 'pos_vec');
    % geT_Nematic_ZGS.m: nem = abs(nanmean(nanmean(order)))
    cm = nan(1, size(order, 2));
    for j = 1:size(order, 2)
        c = order(~isnan(order(:,j)), j);
        if ~isempty(c), cm(j) = mean(c); end
    end
    nemColMean(k) = abs(mean(cm(~isnan(cm))));
    v = order(~isnan(order));
    nemAll(k) = abs(mean(v)); nValid(k) = numel(v);
    if k == 1 && P.makeFigures
        h = figure('Visible', 'off');
        imshow(im, []); set(gca, 'YDir', 'normal'); hold on
        p = [Pt.pos, Pt.vec];
        quiver(p(:,1), p(:,2), p(:,3), p(:,4), 'ShowArrowHead', 'off', 'Color', [1 0 0]);
        hold off
        save_figure(h, fullfile(P.outdir, 'nematicAlignment')); close(h);
    end
end
options = rmfield(P, {'frames', 'outdir', 'makeFigures'}); options.figures = 0; %#ok<NASGU>
save(fullfile(P.outdir, 'alignmentParams.mat'), 'options');
write_csv(fullfile(P.outdir, 'nematic_order.csv'), ...
    {'frame', 'nematic_order', 'nematic_order_allwindows', 'n_valid_windows'}, ...
    [frames(:), nemColMean, nemAll, nValid]);
S.file = [base ext]; S.frames = frames(:); S.nematic_order = nemColMean;
S.nematic_order_allwindows = nemAll; S.outdir = P.outdir;
fprintf('%s: nematic order (frame %d) = %.4f\n', base, frames(1), nemColMean(1));
end
