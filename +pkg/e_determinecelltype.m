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
            % SUM first, so a symbol appearing on more than one row --
            % which readers do produce -- collapses to one profile
            % instead of making MEDIAN return a row vector and the
            % scalar assignment to S(j) below throw. For the usual
            % single-row case sum(.,1) is the identity.
            z = median(sum(Xk(gk == g(ix), :), 1));
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
