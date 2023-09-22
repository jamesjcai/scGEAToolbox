function [g, optimacutoff] = expandSeedGenesNet(A, genelist, seedgenes, k)


[optimacutoff] = fminbnd(@i_totlenproj, 0.65, 1.0, [], ...
    A, genelist, seedgenes, k);
g = i_expandNetOfSeedGenes(A, genelist, seedgenes, optimacutoff);


    function d = i_totlenproj(cutoff, A, genelist, seedgenes, k)
        g = i_expandNetOfSeedGenes(A, genelist, seedgenes, cutoff);
        d = abs(length(g)-k);
end

end


    function g = i_expandNetOfSeedGenes(A, genelist, seedgenes, cutoff)
    if nargin < 4, cutoff = 0.95; end
    if size(genelist, 2) ~= 1, genelist = genelist(:); end
    if size(seedgenes, 2) ~= 1, seedgenes = seedgenes(:); end
    a = max(abs(A(:)));
    A = A ./ a;
    B = A .* (abs(A) > quantile(abs(nonzeros(A)), cutoff));
    g = [];
    for k = 1:length(seedgenes)
        sgene = seedgenes(k);
        [y, i] = ismember(sgene, genelist);
        if y
            idx = (B(:, i) ~= 0 | B(i, :)' ~= 0);
            g = [g; unique(genelist(idx), 'stable')];
        end
    end
    g = unique([seedgenes; g], 'stable');
end
