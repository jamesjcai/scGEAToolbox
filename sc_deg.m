function [T, Tup, Tdn] = sc_deg(X, Y, genelist, methodid, ...
                                guiwaitbar, parentfig)
%SC_DEG - DEG analysis using Mann–Whitney U test or t-test
%
% Inputs:
%   X: matrix for group 1
%   Y: matrix for group 2
%   genelist: list of gene names (optional)
%   methodid: 1 for Mann–Whitney U test (default), 2 for t-test
%   guiwaitbar: true to show progress, false otherwise (default: false)
%
% Outputs:
%   T: table with DEG analysis results
%   Tup, Tdn: processed tables for up- and down-regulated genes
%
% https://satijalab.org/seurat/v3.1/de_vignette.html

    if nargin < 4, methodid = 1; end
    if nargin < 5, guiwaitbar = false; end
    if nargin < 6, parentfig = []; end

    ng = size(X, 1);
    assert(isequal(ng, size(Y, 1)), 'X and Y must have the same number of rows (genes)');
    
    p_val  = zeros(ng, 1);
    stats  = zeros(ng, 1);
        
    nx = size(X, 2);
    ny = size(Y, 2);
    
    if guiwaitbar
        fw = gui.myWaitbar(parentfig);
    end

    % --- Normalize counts (on raw scale, before log transform) ---
    Znorm = sc_norm([X, Y]);
    Xnorm = Znorm(:, 1:nx);
    Ynorm = Znorm(:, nx+1:end);

    % --- Compute avg_1, avg_2, and log2FC on the natural (non-log) scale ---
    %     with pseudocount = 1, matching Seurat convention
    pseudocount = 1;
    avg_1      = mean(Xnorm, 2);
    avg_2      = mean(Ynorm, 2);
    avg_log2FC = log2(avg_1 + pseudocount) - log2(avg_2 + pseudocount);

    % --- Log1p transform for statistical testing ---
    X = log1p(Xnorm);
    Y = log1p(Ynorm);

    % Vectorized pct expressed (>0 on log1p-transformed data)
    pct_1 = sum(X > 0, 2) / nx;
    pct_2 = sum(Y > 0, 2) / ny;

    % Loop through genes
    for k = 1:ng
        if guiwaitbar
            if k / ng > 0.618
                gui.myWaitbar(parentfig, fw, false, '', '', k / ng);
            end
        end
        
        x = X(k, :);
        y = Y(k, :);
        
        switch methodid
            case 1  % Mann–Whitney U test
                [px, ~, tx] = ranksum(x, y);
                p_val(k)  = px;
                stats(k)  = tx.ranksum;
            case 2  % Two-sample t-test
                [~, px, ~, tx] = ttest2(x, y);
                p_val(k)  = px;
                stats(k)  = tx.tstat;
        end
        
    end
    
    % Adjust p-values for multiple comparisons
    if exist('mafdr.m', 'file')
        p_val_adj = mafdr(p_val, 'BHFDR', true);
    else
        [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
    end
    
    if guiwaitbar
        gui.myWaitbar(parentfig, fw); 
    end

    gene       = genelist(:);
    abs_log2FC = abs(avg_log2FC);

    T = table(gene, p_val, avg_log2FC, abs_log2FC, avg_1, ...
              avg_2, pct_1, pct_2, p_val_adj, stats);

    if nargout > 1
        [paramset] = gui.i_degparamset(true, parentfig);
        [Tup, Tdn] = pkg.e_processdetable(T, paramset, parentfig);
    end
end
