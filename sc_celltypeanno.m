function [c_cell_type_tx, cL_cell_type_tx, T] = sc_celltypeanno(X, g, c_cluster_id, species)
    if nargin<4 || isempty(species), species='human'; end
    [T] = run.ml_alona_new(X, g, c_cluster_id,'species', species,'bestonly', true);
    cL_cell_type_tx = table2array(T(:,1:2:end))';
    c_cell_type_tx = cL_cell_type_tx(c_cluster_id);
end