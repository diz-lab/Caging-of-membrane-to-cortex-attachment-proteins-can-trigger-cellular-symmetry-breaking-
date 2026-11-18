function S = temporal_autocorr(tifPath, P)
% TEMPORAL_AUTOCORR  Intensity autocorrelation in time + double-exponential fit.
%
%   S = temporal_autocorr('/data/actin.tif', P)
%
% Algorithm (unchanged from ZGS_temp_autocorr.m):
%   for each lag dF = 1 .. n-1 frames:
%       C(dF) = mean_k < z_k . z_{k+dF} >,  k = 1 .. n-1-dF
%   where z_k = (I_k - mean(I_k)) / std(I_k) (whole-frame z-score).
%   C is normalised by its first point (lag = 1 frame) and fitted over the
%   first fitFraction of the lags to
%       C(t) = a1*exp(t/a2) + (1-a1)*exp(t/a3)      (a2, a3 < 0)
%   (bounded fminsearch, fminsearchbnd). tau = -a2, -a3; the longer tau is
%   reported as tau_long (network relaxation time).
%   Optional P.fitOffset = true adds a constant a4, as in the Methods
%   equation: a1*exp(t/a2) + (1-a1)*exp(t/a3) + a4.
%   R^2 is computed as in the original (MATLAB regress of C on the fitted
%   curve without intercept).
%
% P fields:
%   dt           seconds between frames                (required)
%   frame1       first frame                           [1]
%   frameEnd     last frame (0 = last in file)         [0]
%   fitFraction  fraction of lags used in the fit      [0.3]  (Methods text says ~0.4)
%   a0           initial guess [a1 a2 a3]              [0.9 -2000 -50]
%   fitOffset    add constant a4 to the model          [false]
%   outdir       [stack folder]
%   makeFigures  [true]
%
% Outputs (outdir): <name>_autocorr.mat  (param, allCors, stats - same
% variable names as actinCorrelations2Exp.mat), <name>_autocorr.csv,
% <name>_autocorr.pdf

D = struct('frame1', 1, 'frameEnd', 0, 'fitFraction', 0.3, 'a0', [0.9 -2000 -50], ...
           'fitOffset', false, 'outdir', '', 'makeFigures', true);
P = set_defaults(P, D);
if ~isfield(P, 'dt'), error('temporal_autocorr:dt', 'P.dt (s/frame) is required.'); end
[folder, base] = fileparts(tifPath);
if isempty(P.outdir), P.outdir = folder; end
if ~exist(P.outdir, 'dir'), mkdir(P.outdir); end

nAll = numel(imfinfo(tifPath));
if P.frameEnd == 0 || P.frameEnd > nAll, P.frameEnd = nAll; end
stack = double(read_tiff_stack(tifPath, P.frame1:P.frameEnd));
n_im = P.frameEnd - P.frame1 + 1;
endFit = round(n_im * P.fitFraction); startFit = 1;

totalI = squeeze(sum(sum(stack, 1), 2))';   % intensity per frame (bleaching check)
% z-score each frame once, in place (identical to doing it inside the lag loop)
for k = 1:n_im
    im = stack(:,:,k);
    stack(:,:,k) = (im - mean(im(~isnan(im)))) / std(im(:));
end
Z = reshape(stack, [], n_im); clear stack

dFrame = 1:(n_im-1);
allCors.mean = zeros(size(dFrame)); allCors.std = zeros(size(dFrame));
for ii = 1:numel(dFrame)
    dF = dFrame(ii);
    kk = 1:(n_im-1-dF);
    c = zeros(1, numel(kk));
    for q = 1:numel(kk)
        pr = Z(:,kk(q)) .* Z(:,kk(q)+dF);
        c(q) = mean(pr(~isnan(pr)));
    end
    allCors.mean(ii) = mean(c);       % NaN for the last lag (no pairs), as original
    allCors.std(ii)  = std(c);
end
allCors.time = dFrame .* P.dt;
allCors.mean = allCors.mean ./ allCors.mean(1);

%% Fit
tf = allCors.time(startFit:endFit); cf = allCors.mean(startFit:endFit);
ok = ~isnan(cf); tf = tf(ok); cf = cf(ok);
if P.fitOffset
    model = @(a,t) a(1).*exp(t./a(2)) + (1-a(1)).*exp(t./a(3)) + a(4);
    a0 = [P.a0(1:3), 0]; lb = [0 -Inf -Inf -1]; ub = [1 0 0 1];
else
    model = @(a,t) a(1).*exp(t./a(2)) + (1-a(1)).*exp(t./a(3));
    a0 = P.a0(1:3); lb = [0 -Inf -Inf]; ub = [1 0 0];
end
fun = @(a) sum((cf - model(a, tf)).^2);
a = fminsearchbnd(fun, a0, lb, ub);

% R^2 as regress(allCors.mean', X') in the original: y ~ b*X, no intercept
X = model(a, allCors.time); y = allCors.mean;
m = ~isnan(X) & ~isnan(y);
b = X(m)' \ y(m)';
R2 = 1 - sum((y(m) - b*X(m)).^2) / sum((y(m) - mean(y(m))).^2);
stats = [R2 NaN NaN NaN]; %#ok<NASGU>  (only stats(1) was used)

param = P; param.a = a; param.endFit = endFit; %#ok<NASGU>
param.totalI = totalI;
save(fullfile(P.outdir, [base '_autocorr.mat']), 'param', 'allCors', 'stats');
write_csv(fullfile(P.outdir, [base '_autocorr.csv']), {'lag_s', 'C_norm', 'C_std', 'fit'}, ...
          [allCors.time(:), allCors.mean(:), allCors.std(:), X(:)]);

taus = -a(2:3);
S.file = base; S.a = a; S.tau1_s = taus(1); S.tau2_s = taus(2);
S.tau_long_s = max(taus); S.A = a(1); S.R2 = R2; S.n_frames = n_im;
if P.fitOffset, S.offset = a(4); else, S.offset = 0; end

if P.makeFigures
    h = figure('Visible', 'off'); hold on
    scatter(allCors.time, allCors.mean, 20, 'k', 'filled');
    plot(allCors.time, X, '-r', 'LineWidth', 1.5);
    plot(allCors.time(endFit) * [1 1], [min(y(m)) 1], ':', 'Color', [.5 .5 .5]);
    hold off; box on; axis square
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('time lag \tau (s)'); ylabel('<I(t) I(t+\tau)>');
    title(sprintf('%s  \\tau_{long} = %.3g s', strrep(base,'_','\_'), S.tau_long_s), 'Interpreter', 'tex');
    save_figure(h, fullfile(P.outdir, [base '_autocorr'])); close(h);
end
fprintf('%s: tau = %.4g s, %.4g s (A = %.3f), R2 = %.3f\n', base, taus(1), taus(2), a(1), R2);
end
