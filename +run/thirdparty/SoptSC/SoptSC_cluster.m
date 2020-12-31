function [W,No_cluster,cluster_label,H,eigenvalues] = SoptSC_cluster(data,NC,No_exc_cell,No_features,folder)
% G_filter_idx
% SOptSC identifies clusters from single cell data
%
% Input
%   -data:      An m*n single cell data matrix with m rows(genes) and n columns(cells)
%   -NC:        Number of Cluster specified by user, if NC = [] (default), then the 
%               algorithm will compute the number.
%
%   -No_exc_cell: Gene selection parameter range from [0,No_cells] (No_cells represents
%                 the number of cells), we remove genes that are expressed less than 
%                 No_exc_cell cells and more than (No_cells - No_exc_cell)
%                 cells (Default value: No_exc_cell = 6)
%   -No_features: Gene selection parameter, which represent the number of highly expressed
%                 genes selected for clustering (Default: No_features = 2000)
%                 
%                 
%
% Output
%   --  W: Cell-to-cell similarity matrix.
%   --  No_cluster: Number of clusters
%   --  cluster_label: cluster labels for all cells.
%   --  latent: An nx2 matrix representing low dimensional latent space for all cells
%   --  H: Non-negative matrix such that W = H*H^T
%   --  eigenvalues: eigenvalues of the graph Laplacian of the truncated
%       consensus matrix which is used to estimate the number of clusters 


if nargin==1
    NC = [];
    No_exc_cell = 6;
    No_features = 2000;
elseif nargin == 2
    No_exc_cell = 6;
    No_features = 2000;
elseif nargin == 3
    No_features = 2000;
end
[No_genes,No_cells] = size(data);

% Data preprocess: Gene filtering
if No_exc_cell > 0
    gene_nnz = zeros(No_genes,1);
    alpha_filter = No_exc_cell./No_cells;
    for i = 1:No_genes
        gene_nnz(i) = nnz(data(i,:))./No_cells;
    end
    G_filter_idx1 = union(find(gene_nnz<=alpha_filter),find(gene_nnz>=1-alpha_filter));
    G_filter_idx = setdiff(1:No_genes,G_filter_idx1);
else
    G_filter_idx = 1:size(data,1);
end

data1 = data(G_filter_idx,:);

[coeff, score, pca_eigvalue] = pca(data1');
[~,No_Comps] = max(abs(pca_eigvalue(2:end-1) - pca_eigvalue(3:end)));

aa = max(coeff(:,1:No_Comps+1)');
bb = sort(aa,'descend');

if size(data1,1) <=1000
    No_sel_genes = size(data1,1);
else
    No_sel_genes = min([No_features round(size(data1,1))]);
end

gene_selection = find(aa>=bb(No_sel_genes));
input_data = data1(gene_selection,:);


eigenvalues = [];
nC = NC;
[No_cluster,W,idx,eigenvalues,H] = Main(nC,input_data);

% save eigenvalues of the graph Laplacian of the truncated consensus matrix
if isempty(NC)
    T = table(eigenvalues(1:min([No_cells 100])));
    writetable(T, [folder '/EigenValue.txt'], 'WriteVariableNames',false);
end

cluster_label = idx;
