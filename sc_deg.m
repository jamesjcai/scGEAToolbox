function [T, Tup, Tdn] = sc_deg(X, Y, genelist, methodid, guiwaitbar)
%DEG analysis using Mann–Whitney U test
% https://satijalab.org/seurat/v3.1/de_vignette.html
% p_val : p_val (unadjusted)
% avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
% abs_logFC
% pct.1 : The percentage of cells where the feature is detected in the first group
% pct.2 : The percentage of cells where the feature is detected in the second group
% p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.
%
% SEE ALSO: [T]=run.r_MAST(X,Y,genelist);

if nargin < 2, error("USAGE: sc_deg(X, Y)\n"); end
if nargin < 3, genelist = string(1:size(X, 1))'; end
if nargin < 4, methodid = 1; end
if nargin < 5, guiwaitbar = false; end

ng = size(X, 1);
assert(isequal(ng, size(Y, 1)));

p_val = ones(ng, 1);
avg_log2FC = ones(ng, 1);
avg_1 = zeros(ng, 1);
avg_2 = zeros(ng, 1);
pct_1 = ones(ng, 1);
pct_2 = ones(ng, 1);

nx = size(X, 2);
ny = size(Y, 2);

Z = log(1+sc_norm([X, Y]));
X = Z(:, 1:nx);
Y = Z(:, nx+1:end);


if guiwaitbar
    fw = gui.gui_waitbar_adv;
end
for k = 1:ng
    if guiwaitbar
        gui.gui_waitbar_adv(fw, k/ng);
    end
    x = X(k, :);
    y = Y(k, :);
    switch methodid
        case 1
            % “wilcox” : Wilcoxon rank sum test (default)
            % i.e., Mann–Whitney U test
            p_val(k) = ranksum(x, y);
        case 2
            [~, p] = ttest2(x, y);
            p_val(k) = p;
        otherwise
            error('Unknown option');
    end
    avg_1(k) = mean(x);
    avg_2(k) = mean(y);
    avg_log2FC(k) = log2(avg_1(k)./avg_2(k));
    pct_1(k) = sum(x > 0) ./ nx;
    pct_2(k) = sum(y > 0) ./ ny;
end

if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
end

if guiwaitbar, gui.gui_waitbar_adv(fw); end
% sortid=(1:length(genelist))';
if size(genelist, 2) > 1
    gene = genelist';
else
    gene = genelist;
end
abs_log2FC = abs(avg_log2FC);
T = table(gene, p_val, avg_log2FC, abs_log2FC, avg_1, avg_2, ...
    pct_1, pct_2, p_val_adj);
if nargout > 1
    [Tup, Tdn] = pkg.e_processDETable(T);
end
end
