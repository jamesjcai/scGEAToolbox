function [M, C] = e_cellscorecorrmat(X, g, gsets, methodid, ~)

if nargin<4, methodid = 2; end

C = zeros(size(X, 2), length(gsets));

for k = 1:length(gsets)
    tgsPos = unique(strsplit(string(gsets{k}),','));
    [cs] = sc_cellscore(X, g, tgsPos, [], methodid);
    C(:, k) = cs(:);
end
M = corr(C,'Type','Spearman');

%
%
% for k = 1:n
%     [y{k}, ~, posg] = pkg.e_cellscores(sce.X, sce.g, ...
%         indx2(k), methodid, false);
%     ttxt{k} = T.ScoreType(indx2(k));
% end
%
