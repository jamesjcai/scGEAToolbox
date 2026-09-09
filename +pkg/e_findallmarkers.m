function [T] = e_findallmarkers(X, g, c, cL, logfc, minpct, showwaitbar, ...
                maxnummarkers)

% https://satijalab.org/seurat/reference/findallmarkers
if nargin < 4, cL = []; end
if nargin < 5 || isempty(logfc), logfc = 0.5; end % logfc.threshold
if nargin < 6 || isempty(minpct), minpct = 0.1; end % min pct
if nargin < 7, showwaitbar = false; end
if nargin < 8, maxnummarkers = 100; end

if isempty(cL)
    [c, cL] = findgroups(string(c));
end
if issparse(X), X = full(X); end
X = log1p(sc_norm(X));
if showwaitbar
    fw = gui.gui_waitbar_adv;
end
mC = max(c);
Tcell = cell(mC, 1);
for kc = 1:mC
    if showwaitbar
        if kc ~= mC
            gui.gui_waitbar_adv(fw, kc/mC, sprintf('Processing %s', cL{kc}));
        else
            gui.gui_waitbar_adv(fw, (kc-1)/mC, sprintf('Processing %s', cL{kc}));
        end
    end
    [t] = in_findmarkers(X(:, c == kc), X(:, c ~= kc), cL(kc));
    % IN_FINDMARKERS returns its rows already ranked. It has to: truncating
    % first and sorting afterwards -- which is what this did -- keeps the
    % first MAXNUMMARKERS genes in GENELIST order and throws the real
    % markers away. On a 3000-gene fixture with the 100 strongest markers
    % placed late in the list and 300 weak ones early, the returned "top
    % 100" contained none of the strong ones and all 100 of the weak.
    Tcell{kc} = t(1:min([maxnummarkers, size(t,1)]), :);
end
T = vertcat(Tcell{:});
[~, idx] = sort(T.p_val_adj);
T = T(idx, :);
[~, idx] = natsort(T.grp);
T = T(idx, :);
if showwaitbar
    gui.gui_waitbar_adv(fw);
end


function [t] = in_findmarkers(x, y, ctxt)
ng = size(x, 1);
p_val = ones(ng, 1);
avg_log2FC = ones(ng, 1);
avg_1 = zeros(ng, 1);
avg_2 = zeros(ng, 1);
pct_1 = ones(ng, 1);
pct_2 = ones(ng, 1);
nx = size(x, 2);
ny = size(y, 2);
for k = 1:ng
    xk = x(k, :);
    yk = y(k, :);
    p_val(k) = ranksum(xk, yk);
    % X arrives log1p-transformed, so the means have to be taken back to the
    % normalised count scale before they can be divided. Taking log2 of a
    % ratio of log-means -- which is what this computed -- is not a fold
    % change on any scale: a marker at 200 against 70 counts has a true
    % log2FC of 1.50 and came out as 0.315, below the default LOGFC
    % threshold, so the strongest markers were the ones being dropped.
    % EXPM1 inverts the LOG1P exactly and allocates only a row, which
    % matters on a matrix this size. Pseudocount and form follow SC_DEG so
    % that the two functions report the same number for the same contrast.
    avg_1(k) = mean(expm1(xk));
    avg_2(k) = mean(expm1(yk));
    avg_log2FC(k) = log2(avg_1(k) + 1) - log2(avg_2(k) + 1);
    pct_1(k) = nnz(xk > 0) ./ nx;
    pct_2(k) = nnz(yk > 0) ./ ny;
end
% PKG.E_FDR replaces the two-branch block that used to sit here. Its MAFDR
% branch and its PKG.E_FDR_BH branch are the same algorithm on clean input
% -- they agree to 2e-16 -- but they disagree whenever the p-values carry
% NaN, which is what RANKSUM returns for a gene with no counts in either
% group and what a per-cell-type run produces in bulk. One drops those from
% the family, the other counts them, so the same command gave a different
% answer depending on whether the Bioinformatics Toolbox was installed.
p_val_adj = pkg.e_fdr(p_val);
grp = repmat(ctxt, ng, 1);
t = table(grp, g, p_val, avg_log2FC, avg_1, avg_2, ...
    pct_1, pct_2, p_val_adj);
t = t(p_val_adj < 0.05 & avg_log2FC > logfc & ...
    (pct_1 > minpct | pct_2 > minpct), :);
% Rank here, before the caller truncates to MAXNUMMARKERS. Ties on the
% adjusted p-value are common -- the rank-sum test has a floor set by the
% group sizes -- so the fold change breaks them, as in SC_DEGMAST.
t = sortrows(t, ["p_val_adj", "avg_log2FC"], ["ascend", "descend"]);
end

end
