function plot_gofr_overlay(matFiles, labels, outFile)
% PLOT_GOFR_OVERLAY  Overlay g(r)/g_max curves + fits (generic gofr_plt.m).
%   plot_gofr_overlay({'a_gofr.mat','b_gofr.mat'}, {'0 chol','30 chol'}, 'gr.pdf')
%   plot_gofr_overlay('/folder')   -> every *_gofr.mat in the folder
if ischar(matFiles) && exist(matFiles, 'dir') == 7
    d = dir(fullfile(matFiles, '*_gofr.mat'));
    matFiles = cellfun(@(n) fullfile(d(1).folder, n), {d.name}, 'UniformOutput', false);
end
if nargin < 2 || isempty(labels)
    labels = cellfun(@(f) strrep(regexprep(f, '.*[\\/]|_gofr\.mat', ''), '_', ' '), matFiles, 'UniformOutput', false);
end
h = figure; hold on; cols = lines(numel(matFiles));
for k = 1:numel(matFiles)
    L = load(matFiles{k}); a = L.a; r = L.r_um; fb = L.P.fitBins;
    scatter(r, L.gn, 20, cols(k,:), 'filled', 'DisplayName', labels{k});
    rf = r(fb(1):fb(2));
    plot(rf, a(1).*exp(-rf./a(2)) + a(3), '--', 'Color', cols(k,:), 'HandleVisibility', 'off');
end
hold off; box on; xlabel('r (\mum)'); ylabel('g(r)/g(r)_{max}'); ylim([0 1]);
legend('Location', 'northeast');
if nargin >= 3 && ~isempty(outFile), save_figure(h, outFile); end
end
