function obj = sortcells(obj,idx)
    obj.X=obj.X(:,idx);
    obj.s=obj.s(idx,:);
    obj.c=obj.c(idx);
    if ~isempty(obj.c_cell_cycle_tx)
        obj.c_cell_cycle_tx=obj.c_cell_cycle_tx(idx);
    end
    if ~isempty(obj.c_cell_type_tx)
        obj.c_cell_type_tx=obj.c_cell_type_tx(idx);
    end
    if ~isempty(obj.c_cluster_id)
        obj.c_cluster_id=obj.c_cluster_id(idx);
    end
    if ~isempty(obj.c_batch_id)
        obj.c_batch_id=obj.c_batch_id(idx);
    end
    if ~isempty(obj.c_cell_id)
        obj.c_cell_id=obj.c_cell_id(idx);
    end
    for k=2:2:length(obj.list_cell_attributes)
        obj.list_cell_attributes{k}=obj.list_cell_attributes{k}(idx);
    end
    
    a=fields(obj.struct_cell_embeddings);
    for k=1:length(a)
        if ~isempty(obj.struct_cell_embeddings.(a{k}))
            obj.struct_cell_embeddings.(a{k})=obj.struct_cell_embeddings.(a{k})(idx,:);
        end
    end
    
    a=fields(obj.struct_cell_clusterings);
    for k=1:length(a)
        if ~isempty(obj.struct_cell_clusterings.(a{k}))
             obj.struct_cell_clusterings.(a{k})=obj.struct_cell_clusterings.(a{k})(idx);
        end
    end
end
