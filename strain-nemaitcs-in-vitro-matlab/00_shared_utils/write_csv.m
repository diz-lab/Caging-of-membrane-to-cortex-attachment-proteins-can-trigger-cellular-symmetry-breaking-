function write_csv(fname, headers, data, rowLabels, labelHeader)
% WRITE_CSV  Write a numeric matrix (plus optional text first column) to CSV.
%   write_csv(fname, {'a','b'}, M)
%   write_csv(fname, {'a','b'}, M, {'file1','file2'}, 'file')
% Plain fprintf, so it works in every MATLAB version and in Octave.
if nargin < 4, rowLabels = {}; end
if nargin < 5, labelHeader = 'label'; end
fid = fopen(fname, 'w');
if fid < 0, error('write_csv:open', 'Cannot write %s', fname); end
h = headers;
if ~isempty(rowLabels), h = [{labelHeader}, h]; end
fprintf(fid, '%s\n', strjoin(h, ','));
for i = 1:size(data, 1)
    if ~isempty(rowLabels)
        fprintf(fid, '"%s",', rowLabels{i});
    end
    row = arrayfun(@(v) sprintf('%.10g', v), data(i,:), 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(row, ','));
end
fclose(fid);
end
