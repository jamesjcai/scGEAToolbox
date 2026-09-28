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

if nargin == 1
    % --- Handle SingleCellExperiment input ---
    if isa(X, 'SingleCellExperiment')
        g = X.g;
        X = X.X;
    else
        error('Invalid input(s).')
    end

else

    % --- Input parsing ---
    p = inputParser;
    addRequired(p, 'X', @(x) isnumeric(x) || issparse(x) || isa(x, 'SingleCellExperiment'));
    addOptional(p, 'g', [], @(x) isempty(x) || isstring(x) || iscellstr(x));
    parse(p, X, g);

    X = p.Results.X;
    g = p.Results.g;

end

% --- Gene names default ---
numGenes = size(X,1);
if isempty(g)
    g = pkg.i_defaultgenenames(numGenes);
elseif iscellstr(g)
    g = string(g);
end
g = g(:);  % ensure column

% --- Compute statistics ---
% On X as given. This used to FULL() a sparse X first -- 8 GB at 20000
% genes x 50000 cells -- only to take row statistics, which the sparse
% forms give directly; the per-gene vectors are made full instead.
% PKG.E_ROWVAR, not STD/VAR along dim 2: on a sparse matrix those walk
% every zero (7.0 s against 0.23 s at 20000 genes x 30000 cells); the
% values agree to ~1e-12 relative, and a constant row still gives 0.
dropr = full(1 - sum(X > 0, 2) ./ size(X, 2));
u     = full(mean(X, 2, 'omitnan'));
cv    = sqrt(pkg.e_rowvar(X, "omitnan")) ./ u;

% --- Assemble result table ---
T = table(g, u, cv, dropr, ...
    'VariableNames', {'Gene', 'Mean', 'CV', 'Dropout_rate'});
end
