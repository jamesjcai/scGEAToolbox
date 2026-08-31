function M = i_pairwisemi(codes)
%I_PAIRWISEMI  All-pairs mutual information (bits, Miller-Madow corrected,
%floored at 0) among the rows of CODES, a vars-by-cells matrix of positive
%integer codes. Each row may use its own number of levels, so a gene binned
%into MIBins levels and a metadata column with a handful of categories can
%sit in the same matrix - unlike I_MIVEC's pivot-vs-all trick, which needs
%one shared bin count across every candidate column.
%
%   M = I_PAIRWISEMI(codes)
%
%   Same estimator PKG.SC_SCE2CAUSALCCC documents for gene selection, applied
%   pairwise instead of pivot-vs-all. Variable counts here are small (a
%   side's ligands/receptors + intracellular genes + metadata), so the
%   O(n^2) loop this implies is not a concern the way scoring genes
%   genome-wide would be.
n = size(codes, 1);
M = zeros(n);
for a = 1:n-1
    for b = a+1:n
        m = i_mi2(codes(a, :), codes(b, :));
        M(a, b) = m;
        M(b, a) = m;
    end
end
end


function mi = i_mi2(x, y)
n = numel(x);
nx = double(max(x));
ny = double(max(y));
if nx < 2 || ny < 2
    mi = 0;   % a constant variable shares no information with anything
    return
end

C = zeros(nx, ny);
idx = sub2ind([nx, ny], double(x(:)), double(y(:)));
for k = 1:n
    C(idx(k)) = C(idx(k)) + 1;
end

Cn = C / n;
px = sum(Cn, 2);
py = sum(Cn, 1);
denom = px .* py;
t = Cn .* log2(max(Cn, realmin) ./ max(denom, realmin));
t(C == 0) = 0;
mi = sum(t(:));

nzJoint = sum(C(:) > 0);
nzX = sum(px > 0);
nzY = sum(py > 0);
mi = mi - (nzJoint - nzX - nzY + 1) / (2 * n * log(2));
mi = max(mi, 0);
end
