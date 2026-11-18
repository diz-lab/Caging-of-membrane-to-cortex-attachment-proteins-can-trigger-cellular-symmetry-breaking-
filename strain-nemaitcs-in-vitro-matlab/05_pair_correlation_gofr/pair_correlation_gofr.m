function S = pair_correlation_gofr(tifPath, P)
% PAIR_CORRELATION_GOFR  Radial pair-correlation g(r) of foreground pixels.
%
%   S = pair_correlation_gofr('/data/actin_binary.tif', P)
%
% Algorithm (GrForXlinksAfines.m, V. Yadav; used for the Alba g(r)):
%   foreground pixels of a binarised image are treated as particles.
%   For every "interior" particle i, particles in a (2*rlim)^2 box give the
%   local number density rho_i; g_i(r) = n_i(r, r+dr) / (2*pi*r*dr-area) / rho_i
%   (self excluded at r = 0). g(r) = mean_i g_i(r); last 4 bins dropped.
%   g/g_max is fitted over bins fitBins to  A*exp(-r/xi) + c  (fminsearchbnd).
%
% Interior particles (P.edgeMode):
%   'original' (default): the original test  -dim/2+rlim < x < dim/2-rlim,
%       written for centred simulation coordinates. On images (x = 1..dim)
%       it keeps particles with row < dim_rows/2-rlim and col < dim_cols/2-rlim,
%       i.e. the top-left quarter including pixels at the image edge.
%   'interior': rlim < x <= dim-rlim on both axes (true interior; boxes
%       never leave the image). Recommended for new data.
%
% P fields:
%   pixPerUm    pixels per um (original conv_c)          (required; Alba: 9.2937)
%   frame       frame to analyse                         [1]
%   rlim        max r, px                                [35]
%   dr          bin width, px                            [1]
%   threshold   foreground = I >= threshold. [] -> image must be binary
%               (two grey values); foreground = max value (original: I/255 == 1)
%   edgeMode    'original' | 'interior'                  ['original']
%   fitBins     [first last] bin indices for the fit     [2 14]
%   a0          initial guess [A xi_um c]                [1 1 1]
%   outdir, makeFigures
% Outputs: <name>_gofr.csv (r_um, g, g_norm, fit), <name>_gofr.mat, <name>_gofr.pdf

D = struct('frame', 1, 'rlim', 35, 'dr', 1, 'threshold', [], 'edgeMode', 'original', ...
           'fitBins', [2 14], 'a0', [1 1 1], 'outdir', '', 'makeFigures', true);
P = set_defaults(P, D);
if ~isfield(P, 'pixPerUm'), error('pair_correlation_gofr:scale', 'P.pixPerUm is required.'); end
[folder, base, ext] = fileparts(tifPath);
if isempty(P.outdir), P.outdir = folder; end

im = double(imread(tifPath, P.frame)); if size(im, 3) > 1, im = im(:,:,1); end
if isempty(P.threshold)
    u = unique(im(:));
    if numel(u) > 2
        error('pair_correlation_gofr:binary', ['%s is not binary (%d grey levels). ' ...
              'Binarise it first or set P.threshold.'], base, numel(u));
    end
    fg = im == max(u);
else
    fg = im >= P.threshold;
end
[dim_x, dim_y] = size(im);
[x, y] = find(fg);
d = [x y];
rlim = P.rlim; dr = P.dr;

switch lower(P.edgeMode)
    case 'original'
        inner = d(:,1) > -dim_x/2+rlim & d(:,1) < dim_x/2-rlim & ...
                d(:,2) > -dim_y/2+rlim & d(:,2) < dim_y/2-rlim;
    case 'interior'
        inner = d(:,1) > rlim & d(:,1) <= dim_x-rlim & d(:,2) > rlim & d(:,2) <= dim_y-rlim;
    otherwise
        error('pair_correlation_gofr:edgeMode', 'Unknown edgeMode %s', P.edgeMode);
end
din = d(inner, :);
if isempty(din), error('pair_correlation_gofr:empty', 'No interior particles in %s.', base); end
rr = 0:dr:rlim; nb = numel(rr);
da = pi*((rr+dr).^2 - rr.^2);
% Box around particle i = pixels with |row-xi| < rlim and |col-yi| < rlim
% (strict, as the original find(d > xi-rlim) & find(d < xi+rlim)). Because
% particles are pixels, the box is read straight from the mask; this gives
% the same counts as the original particle-list search, much faster.
o = -(rlim-1):(rlim-1);
[OC, OR] = meshgrid(o, o);                 % OR: row offset, OC: col offset
Dbox = sqrt(OR.^2 + OC.^2);
binBox = floor(Dbox ./ dr) + 1;            % bin c <=> rr(c) <= dist < rr(c)+dr
binBox(Dbox >= rr(end) + dr) = 0;          % beyond the last bin
res = zeros(nb, size(din, 1));
for i = 1:size(din, 1)
    xi = din(i,1); yi = din(i,2);
    r1 = max(1, xi-rlim+1); r2 = min(dim_x, xi+rlim-1);
    c1 = max(1, yi-rlim+1); c2 = min(dim_y, yi+rlim-1);
    W = fg(r1:r2, c1:c2);
    B = binBox((r1:r2)-xi+rlim, (c1:c2)-yi+rlim);
    numdensity = sum(W(:)) / (2*rlim)^2;   % particles in box (incl. i) / box area
    b = B(W & B > 0);
    n = accumarray(b(:), 1, [nb 1]);
    n(1) = n(1) - 1;                       % remove self at r = 0
    res(:, i) = n(:) ./ da(:) ./ numdensity;
end
grfin = mean(res, 2)';
grfin = grfin(1:end-4);
r_um = rr(1:end-4) ./ P.pixPerUm;

st = P.fitBins(1); fi = P.fitBins(2);
gn = grfin ./ max(grfin);
fun = @(a) sum((gn(st:fi) - (a(1).*exp(-r_um(st:fi)./a(2)) + a(3))).^2);
a = fminsearchbnd(fun, P.a0, [0 0 0], [Inf Inf Inf]);
fitc = nan(size(r_um)); fitc(st:fi) = a(1).*exp(-r_um(st:fi)./a(2)) + a(3);

save(fullfile(P.outdir, [base '_gofr.mat']), 'r_um', 'grfin', 'gn', 'a', 'P');
write_csv(fullfile(P.outdir, [base '_gofr.csv']), {'r_um', 'g', 'g_over_gmax', 'fit'}, ...
          [r_um(:), grfin(:), gn(:), fitc(:)]);
S.file = [base ext]; S.r_um = r_um; S.g = grfin; S.a = a; S.xi_um = a(2);
S.n_particles = size(d, 1); S.n_interior = size(din, 1);
if P.makeFigures
    h = figure('Visible', 'off'); hold on; box on
    scatter(r_um, gn, 25, 'k', 'filled');
    plot(r_um(st:fi), fitc(st:fi), '--r', 'LineWidth', 1.5);
    hold off; xlabel('r (\mum)'); ylabel('g(r)/g(r)_{max}'); ylim([0 1]);
    title(sprintf('%s  \\xi = %.3g \\mum', strrep(base,'_','\_'), a(2)), 'Interpreter', 'tex');
    save_figure(h, fullfile(P.outdir, [base '_gofr'])); close(h);
end
fprintf('%s: %d interior particles, xi = %.4g um\n', base, S.n_interior, a(2));
end
