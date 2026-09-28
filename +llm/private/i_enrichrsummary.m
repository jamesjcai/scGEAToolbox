function s = i_enrichrsummary(Tlist, lib_labels)
%I_ENRICHRSUMMARY Condense each Enrichr table to its top 10 terms for JSON.
%   s = i_enrichrsummary(Tlist, lib_labels) has one field per library
%   label, each a cell array of structs with term, p_adj and genes (empty
%   for a library with no results). Shared by the three +llm Enrichr
%   runners, which each carried a copy.
s = struct();
for i = 1:numel(Tlist)
    T = Tlist{i};
    fname = matlab.lang.makeValidName(lib_labels(i));
    if isempty(T) || ~istable(T) || height(T) == 0
        s.(fname) = {};
        continue;
    end
    n = min(10, height(T));
    rows = cell(1, n);
    for j = 1:n
        rows{j} = struct( ...
            'term',  char(T.TermName(j)), ...
            'p_adj', T.AdjustedP_value(j), ...
            'genes', char(T.OverlappingGenes{j}));
    end
    s.(fname) = rows;
end
end
