function S = piv_strain_analysis(tifPath, P)
% PIV_STRAIN_ANALYSIS  PIV (MatPIV / cartPIV) + network strain on one TIFF stack.
%
%   S = piv_strain_analysis('/data/actin.tif', P)
%
% Pipeline (identical to the lab DefineParams.m + runCartPIV.m route):
%   1. cartPIV: MatPIV between frame i and i+imdist, SNR / global / local
%      median filters, NaN interpolation; displacements converted to um.
%   2. PIVstrain: cumulative strain = divergence of the cumulative
%      displacement field (cumsum over frames), one map per frame.
%   3. epsT(t) = -<div(cum. displacement)>  (spatial nanmean, sign flipped
%      so contraction is positive).  Max strain = max epsT(t).
%   4. Strain rate d(epsT)/dt from a linear fit of epsT(t) over a window:
%        'pct25_75' (default) : frames where epsT first reaches 25 % and 75 %
%                               of its maximum (runForStrainAnalysis.m)
%        'maxslope'           : steepest window of P.fitWindow_s seconds
%        'manual'             : P.fitRange_s = [t_start t_end] seconds
%
% P fields (all optional except dscl and tscl):
%   dscl        um per pixel                            (required)
%   tscl        seconds per frame                       (required)
%   winsize     PIV window, px                [32]
%   overlap     PIV window overlap            [0.5]
%   method      MatPIV method                 ['single']
%   imdist      frame spacing of PIV pairs    [1]
%   lastim      last frame to analyse, 0=all  [0]
%   snr, globfiltstd, locfiltstd, kernal, locfiltmeth, naninterp
%               filters                       [1.01, 3, 5, 3, 'median', 1]
%   rateMethod  'pct25_75' | 'maxslope' | 'manual'   ['pct25_75']
%   fitWindow_s window length for 'maxslope', s      [15]
%   fitRange_s  [t0 t1] for 'manual', s              []
%   outdir      output folder        [<stack folder>/<stack name>_PIV]
%   makeFigures true/false           [true]
%
% Outputs written to outdir:
%   ExpParams.mat       expt / pivspec / filtspec used (same format as DefineParams)
%   WS/cartPIV_WS.mat   raw PIV workspace written by cartPIV
%   PIV_strain.mat      x, y, dxt, dyt, vxt, vyt, strain maps, epsT(t), fit
%   strain_vs_time.csv  time, epsT, div mean/std, inst. strain rate, valid fraction
%   strain_vs_time.pdf, strain_heatmap.pdf
% S (returned) holds the per-stack summary.

D = struct('winsize', 32, 'overlap', 0.5, 'method', 'single', 'imdist', 1, ...
    'lastim', 0, 'snr', 1.01, 'globfiltstd', 3, 'locfiltstd', 5, 'kernal', 3, ...
    'locfiltmeth', 'median', 'naninterp', 1, 'rateMethod', 'pct25_75', ...
    'fitWindow_s', 15, 'fitRange_s', [], 'outdir', '', 'makeFigures', true);
P = set_defaults(P, D);
if ~isfield(P, 'dscl') || ~isfield(P, 'tscl')
    error('piv_strain_analysis:params', 'P.dscl (um/px) and P.tscl (s/frame) are required.');
end
[folder, base, ext] = fileparts(tifPath);
if isempty(P.outdir), P.outdir = fullfile(folder, [base '_PIV']); end
if ~exist(P.outdir, 'dir'), mkdir(P.outdir); end

%% Load stack and build the lab parameter structures (as DefineParams.m)
stack = read_tiff_stack(tifPath);
if P.lastim > 0, stack = stack(:,:,1:min(P.lastim, size(stack,3))); end
stackdim = size(stack); if numel(stackdim) == 2, stackdim(3) = 1; end
if stackdim(3) < 2, error('piv_strain_analysis:frames', 'Need >= 2 frames.'); end

expt.path = folder;  expt.filename = [base ext];
expt.dscl = P.dscl;  expt.tscl = P.tscl;
expt.timeind = linspace(0, stackdim(3)-1, stackdim(3));
expt.time = expt.timeind .* expt.tscl;
expt.lastim = 0;  expt.stackdim = stackdim;
expt.cp = [nan, nan];
expt.maxr = (min(stackdim(1:2)) - 200) .* expt.dscl;

pivspec.imdist = P.imdist;  pivspec.lastim = 0;
pivspec.winsize = P.winsize; pivspec.overlap = P.overlap; pivspec.method = P.method;
pivspec.plotoverlay = 0;
filtspec.snr = P.snr; filtspec.globfiltstd = P.globfiltstd;
filtspec.locfiltstd = P.locfiltstd; filtspec.kernal = P.kernal;
filtspec.locfiltmeth = P.locfiltmeth; filtspec.naninterp = P.naninterp;
save(fullfile(P.outdir, 'ExpParams.mat'), 'expt', 'pivspec', 'filtspec');

%% PIV
[x, y, dxt, dyt, vxt, vyt, expt, pivspec, filtspec] = ...
    cartPIV(stack, expt, pivspec, filtspec, P.outdir);
nT = size(dxt, 3);
time = (0:nT-1) * P.tscl;

%% Strain (same calls as cartPIVstats.m, without its write to pwd)
strainMap = PIVstrain(x, y, dxt, dyt);       % div of cumulative displacement
srMap     = PIVstrainrate(x, y, vxt, vyt);   % div of instantaneous velocity
nW = size(strainMap,1) * size(strainMap,2);
divMean = zeros(1, nT); divStd = zeros(1, nT); srMean = zeros(1, nT); validFrac = zeros(1, nT);
for t = 1:nT
    s = strainMap(:,:,t); s = s(~isnan(s));
    divMean(t) = mean(s); divStd(t) = std(s);
    r = srMap(:,:,t); r = r(~isnan(r)); srMean(t) = mean(r);
    validFrac(t) = sum(sum(~isnan(dxt(:,:,t)))) / nW;
end
epsT = -divMean;                                   % contraction positive

%% Strain rate
[mxEps, mxIdx] = max(epsT);
i0 = NaN; i1 = NaN;
switch lower(P.rateMethod)
    case 'pct25_75'
        [~, i0] = min(abs(epsT(1:mxIdx) - 0.25*mxEps));
        [~, i1] = min(abs(epsT(1:mxIdx) - 0.75*mxEps));
    case 'maxslope'
        w = max(2, round(P.fitWindow_s / P.tscl) + 1);   % frames in window
        best = -Inf;
        for a = 1:(nT - w + 1)
            idx = a:(a + w - 1);
            if any(isnan(epsT(idx))), continue; end
            p = polyfit(time(idx), epsT(idx), 1);
            if p(1) > best, best = p(1); i0 = idx(1); i1 = idx(end); end
        end
    case 'manual'
        [~, i0] = min(abs(time - P.fitRange_s(1)));
        [~, i1] = min(abs(time - P.fitRange_s(2)));
    otherwise
        error('piv_strain_analysis:rateMethod', 'Unknown rateMethod %s', P.rateMethod);
end
if isnan(i0) || i1 <= i0
    rate = NaN; pfit = [NaN NaN];
    warning('piv_strain_analysis:fit', '%s: fit window has < 2 frames; strain rate = NaN.', base);
else
    pfit = polyfit(time(i0:i1), epsT(i0:i1), 1);   % slope in 1/s (time in s)
    rate = pfit(1);
end

%% Save
fit.method = P.rateMethod; fit.frames = [i0 i1]; fit.t_s = [time(max(i0,1)) time(max(i1,1))];
fit.strain_rate = rate; fit.p = pfit;
params = P; %#ok<NASGU>
save(fullfile(P.outdir, 'PIV_strain.mat'), 'x', 'y', 'dxt', 'dyt', 'vxt', 'vyt', ...
     'strainMap', 'srMap', 'time', 'epsT', 'divMean', 'divStd', 'srMean', ...
     'validFrac', 'fit', 'params', 'expt', 'pivspec', 'filtspec');
write_csv(fullfile(P.outdir, 'strain_vs_time.csv'), ...
    {'time_s', 'eps_contraction', 'div_cumdisp_mean', 'div_cumdisp_std', ...
     'div_velocity_mean_per_s', 'valid_fraction'}, ...
    [time(:), epsT(:), divMean(:), divStd(:), srMean(:), validFrac(:)]);

S.file = [base ext]; S.n_frames = nT; S.max_strain = mxEps;
S.t_max_strain_s = time(mxIdx); S.strain_rate_per_s = rate;
S.fit_start_s = fit.t_s(1); S.fit_end_s = fit.t_s(2); S.outdir = P.outdir;

%% Figures
if P.makeFigures
    h = figure('Visible', 'off');
    plot(time, epsT, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 3); hold on
    if ~isnan(rate)
        tt = time(i0:i1);
        plot(tt, polyval(pfit, tt), '-r', 'LineWidth', 2);
    end
    hold off; box on
    xlabel('time (s)'); ylabel('\epsilon = -\langle\nabla\cdot x\rangle');
    title(sprintf('%s   d\\epsilon/dt = %.3g s^{-1}', strrep(base,'_','\_'), rate), 'Interpreter', 'tex');
    save_figure(h, fullfile(P.outdir, 'strain_vs_time')); close(h);

    h = figure('Visible', 'off');
    imagesc(x(1,:), y(:,1), -strainMap(:,:,mxIdx)); axis image; colorbar
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title(sprintf('-div(cum. disp.) at t = %.1f s (max \\epsilon)', time(mxIdx)));
    save_figure(h, fullfile(P.outdir, 'strain_heatmap')); close(h);
end
fprintf('%s: max strain %.4g, strain rate %.4g 1/s (fit %.1f-%.1f s)\n', ...
        base, mxEps, rate, S.fit_start_s, S.fit_end_s);
end
