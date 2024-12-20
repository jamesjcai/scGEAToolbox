function [Tct] = e_determinecelltype(sce, ptsSelected, wvalu, wgene, celltypev, markergenev)

T = table(celltypev);
Xk = sce.X(:, ptsSelected);
gk = upper(sce.g);

S = zeros(length(celltypev), 1);
for j = 1:length(celltypev)
    g = strsplit(markergenev(j), ',');
    % g(cellfun('isempty', g)) = [];
    g = g(~cellfun('isempty', g));
    Z = 0;
    ng = 0;
    for ix = 1:length(g)
        if any(g(ix) == wgene) && any(g(ix) == gk)
            wi = wvalu(g(ix) == wgene);
            z = median(Xk(gk == g(ix), :));
            Z = Z + z * wi;
            ng = ng + 1;
        end
    end
    if Z > 0, S(j) = Z ./ nthroot(ng, 3); end
end
if all(S(:) == 0)
    Tct = cell2table({'Unknown', 0});
else
    [~, idx] = sort(S, 'descend');
    T = [T, array2table(S)];
    Tct = T(idx, :);
end
Tct.Properties.VariableNames = {'C1_Cell_Type', 'C1_CTA_Score'};
end
