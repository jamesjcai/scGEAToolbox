function [T] = sc_genestats(X, g)

if isa(X, 'SingleCellExperiment')
    g = X.g;
    X = X.X;
else
    if nargin < 2 || isempty(g), g = "Gene"+string(1:size(X, 1)).'; end
end
if issparse(X)
    try
        X=full(X);
    catch
    end
end

dropr = 1 - sum(X > 0, 2) ./ size(X, 2);
u = mean(X, 2, 'omitnan');
cv = std(X, [], 2, 'omitnan') ./ u;

T = table(g(:), u, cv, dropr);
T.Properties.VariableNames = {'Gene', 'Mean', 'CV', 'Dropout_rate'};
% gui.i_viewtable(T);



