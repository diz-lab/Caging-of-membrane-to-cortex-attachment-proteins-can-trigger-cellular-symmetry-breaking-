function files = list_tiffs(input)
% LIST_TIFFS  Resolve a file, folder, or list into full paths of TIFF stacks.
%   files = list_tiffs('')               -> folder picker (uigetdir)
%   files = list_tiffs('/data/run1.tif') -> {'/data/run1.tif'}
%   files = list_tiffs('/data/cond1')    -> every *.tif / *.tiff in the folder
%   files = list_tiffs({'a.tif','b.tif'})-> the same list, checked to exist
% Files whose name starts with '.' (macOS resource forks) are skipped.
if nargin < 1 || isempty(input)
    input = uigetdir(pwd, 'Select folder with TIFF stacks');
    if isequal(input, 0), files = {}; return; end
end
if iscell(input)
    files = input(:)';
elseif exist(input, 'dir') == 7
    d = [dir(fullfile(input, '*.tif')); dir(fullfile(input, '*.tiff'))];
    names = unique({d.name});
    names = names(~strncmp(names, '.', 1));
    files = cellfun(@(n) fullfile(input, n), names, 'UniformOutput', false);
else
    files = {input};
end
for k = 1:numel(files)
    if exist(files{k}, 'file') ~= 2
        error('list_tiffs:notFound', 'File not found: %s', files{k});
    end
end
end
