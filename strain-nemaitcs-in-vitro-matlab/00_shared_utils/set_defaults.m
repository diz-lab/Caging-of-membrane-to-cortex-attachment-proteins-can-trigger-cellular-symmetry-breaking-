function P = set_defaults(P, D)
% SET_DEFAULTS  Fill missing fields of parameter struct P from defaults D.
if nargin < 1 || isempty(P), P = struct(); end
f = fieldnames(D);
for k = 1:numel(f)
    if ~isfield(P, f{k}) || isempty(P.(f{k}))
        P.(f{k}) = D.(f{k});
    end
end
end
