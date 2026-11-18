function plot_autocorr_overlay(matFiles, labels, outFile)
% PLOT_AUTOCORR_OVERLAY  Overlay several C(tau) curves + fits (generic multi_plt_ZGS).
%   plot_autocorr_overlay({'a_autocorr.mat','b_autocorr.mat'}, {'0 chol','30 chol'}, 'overlay.pdf')
%   plot_autocorr_overlay('/folder')  -> every *_autocorr.mat in the folder
if ischar(matFiles) && exist(matFiles, 'dir') == 7
    d = dir(fullfile(matFiles, '*_autocorr.mat'));
    matFiles = cellfun(@(n) fullfile(d(1).folder, n), {d.name}, 'UniformOutput', false);
end
if nargin < 2 || isempty(labels)
    labels = cellfun(@(f) strrep(regexprep(f, '.*[\\/]|_autocorr\.mat', ''), '_', ' '), matFiles, 'UniformOutput', false);
end
h = figure; hold on; cols = lines(numel(matFiles));
for k = 1:numel(matFiles)
    L = load(matFiles{k}); a = L.param.a; t = L.allCors.time;
    scatter(t, L.allCors.mean, 15, cols(k,:), 'filled', 'DisplayName', labels{k});
    f = a(1).*exp(t./a(2)) + (1-a(1)).*exp(t./a(3));
    if numel(a) > 3, f = f + a(4); end
    plot(t, f, '-', 'Color', cols(k,:), 'HandleVisibility', 'off');
end
hold off; box on; axis square; set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('time lag \tau (s)'); ylabel('<I(t) I(t+\tau)>'); legend('Location', 'southwest');
if nargin >= 3 && ~isempty(outFile), save_figure(h, outFile); end
end
