function [markerlist] = sc_pickmarkers(X, genelist, c, topn, methodid)
    if nargin < 5, methodid = 1; end
    if nargin < 4, topn = 10; end
    assert(isequal(grp2idx(c), c));
    markerlist = cell(max(c), 1);
    
    switch methodid
        case 1 % Fast method
            [idxv] = run.ml_PickMarkers(X, genelist, c, topn);
            for k = 1:max(c)
                idx = idxv(1+(k - 1)*topn:k*topn);
                markerlist{k} = genelist(idx(~isnan(idx)));
            end
        case 2
            num_markers = topn * numel(unique(c));
            markerlist = run.ml_scGeneFit(X, genelist, c, num_markers);
        case 3 % LASSO (slower method)
            for k = 1:max(c)
                fprintf('Processing cell group ... %d of %d\n', k, max(c));
                markerlist{k} = i_pickmarkerslasso(X, genelist, c, k, topn);
            end
        case 4 % Slowest method
            for k = 1:max(c)
                a = i_pickmarkers(X, genelist, c, k);
                markerlist{k} = a(1:topn);
            end
    end
end


function [markerlist] = i_pickmarkerslasso(X, genelist, idv, id, topn)
    idx = idv == id;
    y = double(idx);
    if issparse(X)
        X = full(X);
    end
    [B] = lasso(X', y, 'DFmax', topn*3, 'MaxIter', 1e3);
    [~, ix] = min(abs(sum(B > 0)-topn));
    b = B(:, ix);
    idx = b > 0;
    if ~any(idx)
        warning('No marker gene found')
        markerlist = [];
        return;
    else
        markerlist = genelist(idx);
        [~, jx] = sort(b(idx), 'descend');
        markerlist = markerlist(jx);
    end
end


function [markerlist, A] = i_pickmarkers(X, genelist, idv, id)
    % IDV - cluster ids of cells
    % ID  - the id of the cluster, for which marker genes are being identified.
    % see also: run.celltypeassignation
    % Demo:
    %gx=sc_pickmarkers(X,genelist,cluster_id,2);
    %run.celltypeassignation(gx)
    X = sc_transform(X);
    K = max(idv);
    idx = idv == id;
    
    x1 = X(:, idx);
    x0 = X(:, ~idx);
    % A=[];
    T = i_sc_deg(x0, x1, genelist);
    A = T.z_val;
    totn = sum(~idx);
    for k = 1:K
        if k ~= id
            fprintf('Comparing group #%d with group #%d (out of %d)\n', ...
                id, k, K-1);
            x0 = X(:, idv == k);
            T = i_sc_deg(x0, x1, genelist);
            %a=-log(T.p_val).*sign(T.avg_logFC);
    
            w = sum(idv == k) ./ totn; % weight by number of cells
            a = w * T.z_val;
            % a=T.z_val;
    
            A = [A, a];
        end
    end
    % [~,idx]=sort(sum(A,2));
    %A(isnan(A))=0;
    %[~,idx]=sort(vecnorm(A,2,2),'descend');  % NaN messed up
    %[~,idx]=sort(-vecnorm(A,2,2));  % NaN is ignored
    [~, idx] = sort(sum(A, 2, 'omitnan'));
    markerlist = genelist(idx);
end

function [T] = i_sc_deg(X, Y, genelist)
    ng = size(X, 1);
    assert(isequal(ng, size(Y, 1)));
    
    p_val = ones(ng, 1);
    avg_logFC = ones(ng, 1);
    pct_1 = ones(ng, 1);
    pct_2 = ones(ng, 1);
    z_val = ones(ng, 1);
    parfor k = 1:ng
        x = X(k, :);
        y = Y(k, :);
        [xp, ~, xt] = ranksum(x, y, 'method', 'approximate');
        p_val(k) = xp;
        z_val(k) = xt.zval;
        avg_logFC(k) = log2(mean(x)./mean(y));
        pct_1(k) = sum(x > 0) ./ length(x);
        pct_2(k) = sum(y > 0) ./ length(y);
    end
    
    if exist('mafdr.m', 'file')
        p_val_adj = mafdr(p_val, 'BHFDR', true);
    else
        [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
    end
    sortid = (1:length(genelist))';
    if size(genelist, 2) > 1
        gene = genelist';
    else
        gene = genelist;
    end
    T = table(sortid, gene, p_val, avg_logFC, ...
        pct_1, pct_2, p_val_adj, z_val);
    % T=sortrows(T,'p_val_adj','ascend');
end
