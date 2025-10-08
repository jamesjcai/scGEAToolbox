function [sce] = sc_mergesces(sces, method, keepbatchid)
%Merges two SCE objects
%Usage: [sce]=sc_mergesces({sce1,sce2},'intersect');
%See also: SC_MERGEDATA

if nargin < 3, keepbatchid = true; end
if nargin < 2 || isempty(method), method = 'intersect'; end
validMethods = ["intersect", "union"];
method = validatestring(method, validMethods);
if ~iscell(sces), error('SCES is not a cell array.'); end
if length(sces) < 2, error('At least two SCE required.'); end
for k = 1:length(sces)
    if ~isa(sces{k}, 'SingleCellExperiment')
        error('sces{%d} is not a SingleCellExperiment object.', k);
    end
end

needappendix=false;
sce = sces{1};
c = ones(sce.NumCells, 1);
for k = 2:length(sces)
    c = [c; k * ones(sces{k}.NumCells, 1)];
    [sce, hasidoverlapx] = i_merge2sces(sce, sces{k}, method);
    if hasidoverlapx
        needappendix=true;
    end
end
if ~keepbatchid || isscalar(unique(sce.c_batch_id)) 
    sce.c_batch_id = c; 
end
if needappendix
    sce.c_batch_id = strcat(string(sce.c_batch_id), "_", string(c));
    disp('A suffix is added to SCE.C_BATCH_ID to distinguish cells'' original batch IDs.');
end

end


function [sce, hasidoverlap] = i_merge2sces(sce1, sce2, method)
    hasidoverlap = false;
    if nargin < 3, method = 'intersect'; end
    [X, g, ~] = sc_mergedata(sce1.X, sce2.X, ...
        sce1.g, sce2.g, method);
    sce = SingleCellExperiment(X, g);
    sce.metadata = [sce1.metadata; sce2.metadata];
    sce.c = [sce1.c; sce2.c];
    try
        sce.s = [sce1.s; sce2.s];
    catch
        sce.s = randn(size(X, 2), 3);
    end
    % sce.c_batch_id=c;
    
    
    if ~isempty(sce1.c_batch_id) && ~isempty(sce2.c_batch_id)
        if ~isstring(sce1.c_batch_id)
            sce1.c_batch_id = string(sce1.c_batch_id);
        end
        if ~isstring(sce2.c_batch_id)
            sce2.c_batch_id = string(sce2.c_batch_id);
        end
        sce.c_batch_id = [sce1.c_batch_id; sce2.c_batch_id];
        % intersect(sce1.c_batch_id, sce2.c_batch_id)
        % pause
        hasidoverlap = ~isempty(intersect(sce1.c_batch_id, sce2.c_batch_id));
    end
    
    if ~isempty(sce1.c_cell_cycle_tx) && ~isempty(sce2.c_cell_cycle_tx)
        sce.c_cell_cycle_tx = [sce1.c_cell_cycle_tx; sce2.c_cell_cycle_tx];
    end
    
    if ~isempty(sce1.c_cell_type_tx) && ~isempty(sce2.c_cell_type_tx)
        if ~isstring(sce1.c_cell_type_tx)
            sce1.c_cell_type_tx = string(sce1.c_cell_type_tx);
        end
        if ~isstring(sce2.c_cell_type_tx)
            sce2.c_cell_type_tx = string(sce2.c_cell_type_tx);
        end
        sce.c_cell_type_tx = [i_remove_affix(sce1.c_cell_type_tx); ...
            i_remove_affix(sce2.c_cell_type_tx)];
    end
    
    if ~isempty(sce1.c_cluster_id) && ~isempty(sce2.c_cluster_id)
        sce.c_cluster_id = [sce1.c_cluster_id; sce2.c_cluster_id];
    end
    
    if ~isempty(sce1.c_cell_id) && ~isempty(sce2.c_cell_id)
        sce.c_cell_id = [sce1.c_cell_id; sce2.c_cell_id];
    end
    
    try
        if ~isempty(sce1.list_cell_attributes) && ~isempty(sce2.list_cell_attributes)
            for k = 1:2:min([length(sce1.list_cell_attributes), length(sce2.list_cell_attributes)])
                if strcmp(sce1.list_cell_attributes{k}, ...
                        sce2.list_cell_attributes{k})
                    sce.list_cell_attributes{k} = sce1.list_cell_attributes{k};
                    sce.list_cell_attributes{k+1} = [sce1.list_cell_attributes{k+1}; ...
                        sce2.list_cell_attributes{k+1}];
                end
            end
        end
    catch
        warning('SCE.LIST_CELL_ATTRIBUTES not merged.')
    end
    
    try
        a = fields(sce.struct_cell_clusterings);
        for k = 1:length(a)
            c1 = sce1.struct_cell_clusterings.(a{k});
            c2 = sce2.struct_cell_clusterings.(a{k});
            if ~isempty(c1) && ~isempty(c2)
                sce.struct_cell_clusterings.(a{k}) = [c1; c2];
            end
        end
    catch
        warning('SCE.STRUCT_CELL_CLUSTERINGS not merged.')
    end
    
    try
        a = fields(sce.struct_cell_embeddings);
        for k = 1:length(a)
            c1 = sce1.struct_cell_embeddings.(a{k});
            c2 = sce2.struct_cell_embeddings.(a{k});
            if ~isempty(c1) && ~isempty(c2)
                sce.struct_cell_embeddings.(a{k}) = [c1; c2];
            end
        end
    catch
        % warning('SCE.STRUCT_CELL_EMBEDDINGS not merged.')
    end

end


function b = i_remove_affix(a)
    b = a;
    for k = 1:length(a)
        idx = strfind(a(k), '_{');
        if idx > 0
            b(k) = extractBefore(a(k), idx);
        end
    end
end
