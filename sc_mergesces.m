function [sce]=sc_mergesces(sces,method)

if nargin<2, method='intersect'; end
if ~iscell(sces), error('SCES is not a cell array'); end
sce=sces{1};
c=ones(sce.NumCells,1);
for k=2:length(sces)
    c=[c; k*ones(sces{k}.NumCells,1)];
    [sce]=i_merge2sces(sce,sces{k},method);
end
sce.c_batch_id=c;
end


function [sce]=i_merge2sces(sce1,sce2,method)

if nargin<3, method='intersect'; end
[X,g,c]=sc_mergedata(sce1.X,sce2.X,...
            sce1.g,sce2.g,method);
sce=SingleCellExperiment(X,g);
sce.c=[sce1.c; sce2.c];
sce.s=[sce1.s; sce2.s];
sce.c_batch_id=c;
if ~isempty(sce1.c_cell_cycle_tx) && ~isempty(sce2.c_cell_cycle_tx)
    sce.c_cell_cycle_tx=[sce1.c_cell_cycle_tx; sce2.c_cell_cycle_tx];
end
if ~isempty(sce1.c_cell_type_tx) && ~isempty(sce2.c_cell_type_tx)
    sce.c_cell_type_tx=[sce1.c_cell_type_tx; sce2.c_cell_type_tx];
end
if ~isempty(sce1.c_cluster_id) && ~isempty(sce2.c_cluster_id)    
    sce.c_cluster_id=[sce1.c_cluster_id; sce2.c_cluster_id];
end
if ~isempty(sce1.c_cell_id) && ~isempty(sce2.c_cell_id)
    sce.c_cell_id=[sce1.c_cell_id; sce2.c_cell_id];
end
if ~isempty(sce1.list_cell_attributes) && ~isempty(sce2.list_cell_attributes)
    for k=1:2:min([length(sce1.list_cell_attributes) length(sce2.list_cell_attributes)]) 
    	if strcmp(sce1.list_cell_attributes{k},sce2.list_cell_attributes{k})
            sce.list_cell_attributes{k}=sce1.list_cell_attributes{k};
            sce.list_cell_attributes{k+1}=[sce1.list_cell_attributes{k+1}; sce2.list_cell_attributes{k+1}];
        end
    end
end
end

