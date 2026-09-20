function ari = i_adjustedrandindex(a, b)
%I_ADJUSTEDRANDINDEX Adjusted Rand index between two groupings of the same items.
%
%   ari = pkg.i_adjustedrandindex(a, b)
%
%   A and B are label vectors of the same length -- string, numeric,
%   categorical, or a cell array of char -- giving each item's group under two
%   groupings. ARI is the Hubert & Arabie (1985) adjusted Rand index: 1 for
%   identical partitions, 0 for the agreement two independent groupings reach
%   by chance, and negative for agreement below chance.
%
%   The index counts PAIRS of items that the two groupings place together or
%   apart, so it never compares one label's text with the other's. That is why
%   it belongs next to the per-cell agreement rate rather than replacing it:
%   two annotation methods that call the same population "CD8 T cell" and
%   "T cells, CD8+" agree on 0% of cells and still score ARI 1.
%
%   Missing labels form their own group rather than being dropped, so every
%   item contributes and the result does not quietly change meaning with the
%   number of unlabelled cells.
%
%   See also PKG.I_CELLTYPEHISTORY, GUI.CALLBACK_COMPARECELLTYPEANNOTATIONS.

a = in_aslabels(a);
b = in_aslabels(b);
if numel(a) ~= numel(b)
    error('pkg:i_adjustedrandindex:sizeMismatch', ...
        ['A and B have %d and %d elements. Both must label the same ', ...
        'items, so pass two vectors of equal length.'], numel(a), numel(b));
end

n = numel(a);
if n < 2
    % The index counts pairs, and fewer than two items form none.
    ari = NaN;
    return;
end

ga = findgroups(a);
gb = findgroups(b);

% Sparse: the contingency matrix is numTypes(a) x numTypes(b), which stays
% small for cell type labels but is unbounded for a per-cell id used as a
% grouping, and only the nonzero cells contribute to the pair counts.
M = accumarray([ga, gb], 1, [max(ga), max(gb)], [], 0, true);

sumBoth = sum(in_pairs(nonzeros(M)));
sumA = sum(in_pairs(full(sum(M, 2))));
sumB = sum(in_pairs(full(sum(M, 1))));
total = in_pairs(n);

expected = sumA*sumB/total;
maxIndex = 0.5*(sumA + sumB);

if maxIndex == expected
    % Reachable only when both groupings are trivial in the same way -- one
    % group holding every item, or every item alone -- which makes them the
    % same partition. The ratio below is 0/0 there.
    ari = 1;
    return;
end

ari = (sumBoth - expected)/(maxIndex - expected);
end

function p = in_pairs(k)
% Number of unordered pairs within a group of K items, elementwise.
p = k.*(k - 1)/2;
end

function v = in_aslabels(v)
% One representation for every accepted label type. Missing values are named
% rather than dropped, so findgroups gives them a group instead of a NaN
% index that ACCUMARRAY would reject.
v = string(v(:));
v(ismissing(v)) = "<missing>";
end
