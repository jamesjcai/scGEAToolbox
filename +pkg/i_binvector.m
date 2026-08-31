function b = i_binvector(x, nb)
%PKG.I_BINVECTOR  Bin 1 is "not detected"; the detected values are split
%into nb-1 equal-frequency bins. Quantile-binning the raw column instead
%would drop most cells into a single bin whenever the gene is sparse, which
%is most genes. Shared by the causalCCC gene-selection
%(PKG.I_CAUSALCCCASSEMBLE) and the local network builder (SC_CAUSALCCCNET,
%via the qualified name PKG.I_BINVECTOR), so both bin expression the same
%way. Not in a "private" folder for the same reason as PKG.I_CAUSALCCCASSEMBLE.
x = double(x(:));
b = ones(numel(x), 1, 'uint8');
nz = x > 0;
if ~any(nz), return, end

k = nb - 1;
v = x(nz);
if numel(unique(v)) <= k
    [~, ~, idx] = unique(v);
    b(nz) = uint8(idx + 1);
    return
end

edges = unique(quantile(v, linspace(0, 1, k + 1)));
if numel(edges) < 2
    b(nz) = 2;
    return
end
idx = discretize(v, edges);
idx(isnan(idx)) = numel(edges) - 1;    % the max lands on the last edge
b(nz) = uint8(idx + 1);
end
