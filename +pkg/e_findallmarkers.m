function [T] = e_findallmarkers(X, g, c, cL, logfc, minpct, guiwaitbar)

% https://satijalab.org/seurat/reference/findallmarkers
if nargin < 4, cL = []; end
if nargin < 5 || isempty(logfc), logfc = 0.5; end % logfc.threshold
if nargin < 6 || isempty(minpct), minpct = 0.1; end % min pct
if nargin < 7, guiwaitbar = false; end

if isempty(cL)
    [c, cL] = grp2idx(c);
end

X = log(1+sc_norm(X));
T = table();
if guiwaitbar
    fw = gui.gui_waitbar_adv;
end
mC = max(c);
for kc = 1:mC
    if guiwaitbar
        if kc ~= mC
            gui.gui_waitbar_adv(fw, kc/mC, sprintf('Processing %s', cL{kc}));
        else
            gui.gui_waitbar_adv(fw, kc-1/mC, sprintf('Processing %s', cL{kc}));
        end
    end
    [t] = findmarkers(X(:, c == kc), X(:, c ~= kc), cL(kc));
    T = [T; t];
end
[~, idx] = sort(T.p_val_adj);
T = T(idx, :);
[~, idx] = natsort(T.grp);
T = T(idx, :);
if guiwaitbar
    gui.gui_waitbar_adv(fw);
end


    function [t] = findmarkers(x, y, ctxt)
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
            p_val(k) = ranksum(x(k, :), y(k, :));
            avg_1(k) = mean(x(k, :));
            avg_2(k) = mean(y(k, :));
            avg_log2FC(k) = log2(avg_1(k)./avg_2(k));
            pct_1(k) = sum(x(k, :) > 0) ./ nx;
            pct_2(k) = sum(y(k, :) > 0) ./ ny;
        end
        if exist('mafdr.m', 'file')
            p_val_adj = mafdr(p_val, 'BHFDR', true);
        else
            [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
        end
        grp = repmat(ctxt, ng, 1);
        t = table(grp, g, p_val, avg_log2FC, avg_1, avg_2, ...
            pct_1, pct_2, p_val_adj);
        t = t(p_val_adj < 0.05 & avg_log2FC > logfc & ...
            (pct_1 > minpct | pct_2 > minpct), :);
end
end
