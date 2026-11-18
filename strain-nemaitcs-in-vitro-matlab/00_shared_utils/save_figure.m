function save_figure(h, fname)
% SAVE_FIGURE  Save figure h as PDF (and .fig in MATLAB); never throws.
try
    [d, b] = fileparts(fname);
    print(h, fullfile(d, [b '.pdf']), '-dpdf');
    if ~exist('OCTAVE_VERSION', 'builtin')
        savefig(h, fullfile(d, [b '.fig']));
    end
catch err
    warning('save_figure:failed', 'Could not save %s: %s', fname, err.message);
end
end
