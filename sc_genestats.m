function T = sc_genestats(X, g)
% SC_GENESTATS  Compute per-gene statistics into a tidy table
%
%   T = sc_genestats(X) accepts raw counts matrix X (genes×cells) and
%   assigns default gene names.
%
%   T = sc_genestats(X, g) uses provided gene labels (string or cellstr).
%
%   T = sc_genestats(SingleCellExperiment) extracts data and gene names
%   automatically.

    % --- Input parsing ---
    p = inputParser;
    addRequired(p, 'X', @(x) isnumeric(x) || issparse(x) || isa(x, 'SingleCellExperiment'));
    addOptional(p, 'g', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
    parse(p, X, g);

    X = p.Results.X;
    g = p.Results.g;

    % --- Handle SingleCellExperiment input ---
    if isa(X, 'SingleCellExperiment')
        g = X.g;
        X = X.X;
    end

    % --- Gene names default ---
    numGenes = size(X,1);
    if isempty(g)
        g = "Gene" + string((1:numGenes).');
    elseif iscellstr(g)
        g = string(g);
    end
    g = g(:);  % ensure column

    % --- Convert sparse to dense if needed ---
    if issparse(X)
        try
            X = full(X);
        catch
            warning('Could not convert sparse matrix to full; computations may be slow.');
        end
    end

    % --- Compute statistics ---
    dropr = 1 - sum(X > 0, 2) ./ size(X, 2);
    u     = mean(X, 2, 'omitnan');
    cv    = std(X, 0, 2, 'omitnan') ./ u;

    % --- Assemble result table ---
    T = table(g, u, cv, dropr, ...
        'VariableNames', {'Gene', 'Mean', 'CV', 'Dropout_rate'});
end

%{
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
%}


