function S = bundling_parameter(tifPath, P)
% BUNDLING_PARAMETER  lambda = nu / mu of the thresholded actin image, per frame.
%
%   S = bundling_parameter('/data/actin.tif', P)
%
% Two thresholding variants, both from the original scripts:
%   'normalize_total' (default; Methods text; spatial_linescan_analysis/
%       image_threshold_method/threshold_code_ZGS.m):
%       pixels == 0 -> NaN;  I <- I / sum(I);  I < 0 -> NaN
%   'subtract_mean'  (image_bundle_threshold_ZGS/threshold_code_ZGS.m):
%       I <- I - subtractFraction*mean(I);  I < 0 -> NaN   (zeros kept)
% Statistics:
%   statMode 'original' (default, reproduces earlier numbers):
%       mu = nanmean(nanmean(I))  (mean of column means)
%       nu = nanstd(nanstd(I))    (std ACROSS COLUMNS of the per-column std)
%   statMode 'global': mu, nu = mean and std over all non-NaN pixels
%       (the quantity described in the Methods text).
%   lambda = nu / mu. Note: normalising by sum(I) does not change lambda.
%
% P fields: method ['normalize_total'], subtractFraction [0.9],
%           statMode ['original'], frames [all], outdir [stack folder]
% Output: <name>_bundling.csv (frame, mu, nu, lambda); S.lambda per frame.

D = struct('method', 'normalize_total', 'subtractFraction', 0.9, ...
           'statMode', 'original', 'frames', [], 'outdir', '');
P = set_defaults(P, D);
[folder, base, ext] = fileparts(tifPath);
if isempty(P.outdir), P.outdir = folder; end
nAll = numel(imfinfo(tifPath));
frames = P.frames; if isempty(frames), frames = 1:nAll; end
res = nan(numel(frames), 4);
for k = 1:numel(frames)
    I = double(imread(tifPath, frames(k)));
    if size(I, 3) > 1, I = I(:,:,1); end
    switch lower(P.method)
        case 'normalize_total'
            I(I == 0) = NaN;
            I = I ./ sum(I(~isnan(I)));
            I(I < 0) = NaN;
        case 'subtract_mean'
            I = I - P.subtractFraction * nanmean_vec(col_nanmean(I));
            I(I < 0) = NaN;
        otherwise
            error('bundling_parameter:method', 'Unknown method %s', P.method);
    end
    switch lower(P.statMode)
        case 'original'
            mu = nanmean_vec(col_nanmean(I));
            nu = nanstd_vec(col_nanstd(I));
        case 'global'
            v = I(~isnan(I)); mu = mean(v); nu = std(v);
        otherwise
            error('bundling_parameter:statMode', 'Unknown statMode %s', P.statMode);
    end
    res(k,:) = [frames(k), mu, nu, nu/mu];
end
write_csv(fullfile(P.outdir, [base '_bundling.csv']), {'frame', 'mu', 'nu', 'lambda'}, res);
S.file = [base ext]; S.frames = res(:,1); S.mu = res(:,2); S.nu = res(:,3); S.lambda = res(:,4);
fprintf('%s: lambda (frame %d) = %.4f\n', base, frames(1), S.lambda(1));
end

function m = col_nanmean(I)   % = nanmean(I) for a matrix
m = nan(1, size(I, 2));
for j = 1:size(I, 2), c = I(~isnan(I(:,j)), j); if ~isempty(c), m(j) = mean(c); end, end
end
function s = col_nanstd(I)    % = nanstd(I) for a matrix (N-1 normalisation)
s = nan(1, size(I, 2));
for j = 1:size(I, 2), c = I(~isnan(I(:,j)), j); if numel(c) > 0, s(j) = std(c); end, end
end
function m = nanmean_vec(v), v = v(~isnan(v)); m = mean(v); end
function s = nanstd_vec(v),  v = v(~isnan(v)); s = std(v);  end
