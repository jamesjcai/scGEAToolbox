function [M] = e_cellscorecorrmat(X, g, gsets, methodid)

%   [posg, ctselected] = gui.i_selectMSigDBGeneSets(speciestag, false, FigureHandle);

C = zeros(size(X, 2), length(gsets));

for k = 1:length(gsets)
    if methodid == 1
        [cs] = sc_cellscore_ucell(X, g, gsets(k));
    elseif methodid == 2
        [cs] = sc_cellscore_admdl(X, g, gsets(k));
    end
    C(:, k) = cs(:);
end
M = corr(C);



for k = 1:n
    [y{k}, ~, posg] = pkg.e_cellscores(sce.X, sce.g, ...
        indx2(k), methodid, false);
    ttxt{k} = T.ScoreType(indx2(k));
end
    