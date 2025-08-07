function [gene_idxv, T] = ml_PickMarkers(X, genelist, c, topn, plotit)
%PICKMARKERS - adapted from GC_htmp_DE.m

% This is a support function for SC_PICKMARKERS

if nargin < 5, plotit = false; end
if nargin < 4, topn = 10; end

%% new data visualization
numC = length(unique(c));
%gene_idxv = [];
cluster_order = [];

[No_gene] = size(X, 1);
% calculat mean of gene expression
%gene_mean = zeros(No_gene,numC);
gene_DE_score = zeros(No_gene, numC);

% gene_value_idx = zeros(No_gene,1);
for i = 1:numC
    % gene_mean(:,i) = mean(X(:,c==i),2);
    cluster_order = [cluster_order; find(c == i)];
end
if issparse(X)
    X = full(X);
end
% gene_mean = grpstats(X', c, @mean)';
gene_mean = splitapply(@mean, X', c)';

[~, gene_value_idx] = max(gene_mean, [], 2);

% compute DE-score for each gene
for i = 1:numC
    zz = abs(gene_mean(:, i)-gene_mean);
    gene_DE_score(:, i) = sum(zz, 2); % sum
end

%%
% topn markers for each cluster based on DE score
gclusters = [];
gscore = [];
gene_idxv = nan(numC*topn, 1);
gclusters = nan(numC*topn, 1);
gscore = nan(numC*topn, 1);
for i = 1:numC
    zz_idx = find(gene_value_idx == i);
    zz_DEscore = gene_DE_score(zz_idx, i);
    [zzvalue, zz1] = sort(zz_DEscore, 'descend');

    mtopn = min([length(zz1), topn]);
    ts = (i - 1) * topn + 1;
    idx = ts:(ts + mtopn - 1);
    gene_idxv(idx) = zz_idx(zz1(1:mtopn));
    gclusters(idx) = i .* ones(mtopn, 1);
    gscore(idx) = zzvalue(1:mtopn);
end

%%
if any(isnan(gene_idxv))
    genelistplus = [genelist; "---"];
    gene_idxvplus = gene_idxv;
    gene_idxvplus(isnan(gene_idxv)) = 1 + length(genelist);
else
    genelistplus = genelist;
    gene_idxvplus = gene_idxv;
end


GL500 = [gene_idxv, gclusters, gscore];
T = array2table(GL500, 'VariableNames', ...
    {'Gene_indx', 'Cluster', 'DE_Score'});
T = addvars(T, string(genelistplus(gene_idxvplus)), ...
    'NewVariableNames', 'Gene_Name', 'Before', 'Gene_indx');

if plotit
    figure;
    idata = X(gene_idxv, cluster_order);
    kk = 2;
    center = mean(idata, kk);
    scale = std(idata, 0, kk);
    tscale = scale;
    %=Check for zeros and set them to 1 so not to scale them.
    scale(tscale == 0) = 1;
    %== Center and scale the data
    idata = bsxfun(@minus, idata, center);
    sdata = bsxfun(@rdivide, idata, scale);
    thresh = 3;
    colormap redbluecmap;
    clims = [-thresh, thresh];
    imagesc(sdata, clims);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);

    lgd = cell(1, numC);
    for i = 1:numC
        if i < 10
            vv = 'CC';
            vv(2:2) = num2str(i);
            lgd{i} = vv;
        else
            vv = 'CCC';
            vv(2:3) = num2str(i);
            lgd{i} = vv;
        end
    end
    No_cells_inC = [];
    for i = 1:numC
        No_cells_inC = [No_cells_inC; length(find(c == i))];
    end
    xtkval = cumsum(No_cells_inC);
    xtkval1 = zeros(size(xtkval));
    for i = 1:numC
        if i == 1
            xtkval1(i) = 0.5 .* No_cells_inC(i);
        else
            xtkval1(i) = xtkval(i-1) + 0.5 .* No_cells_inC(i);
        end
    end

    if size(idata, 1) < 200
        yticks(1:size(idata, 1));
        yticklabels(genelist(gene_idxv));
    end

    xticks(xtkval1);
    xticklabels(lgd);
    cb = colorbar;
    ax = gca;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5 * cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
end
end

%print([folder '\HeatMap_Top' num2str(topn)],'-dpdf','-r300','-fillpage'); %'-dpdf',
