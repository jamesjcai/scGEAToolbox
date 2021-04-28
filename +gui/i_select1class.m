function [thisc,clable]=i_select1class(sce)

thisc=[];
clable='';

    listitems={'Current Class (C)'};
    if ~isempty(sce.c_cluster_id)
        listitems=[listitems,'Cluster ID'];
    end
    if ~isempty(sce.c_cell_type_tx)
        listitems=[listitems,'Cell Type'];
    end
    if ~isempty(sce.c_cell_cycle_tx)
        listitems=[listitems,'Cell Cycle Phase'];
    end
    if ~isempty(sce.c_batch_id)
        listitems=[listitems,'Batch ID'];
    end
    



% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select Class','',''},...
     'SelectionMode','single','ListString',listitems);
if tf2==1
    clable=listitems{indx2};
    switch clable
        case 'Current Class (C)'
            thisc=sce.c;            
        case 'Cluster ID' % cluster id
            thisc=sce.c_cluster_id;
        case 'Batch ID' % batch id
            thisc=sce.c_batch_id;
        case 'Cell Type' % cell type
            thisc=sce.c_cell_type_tx;                
        case 'Cell Cycle Phase' % cell cycle
            thisc=sce.c_cell_cycle_tx;
    end    
end

% if isempty(thisc)
%     errordlg('Undefined');    
%     return;
% end
% if numel(unique(thisc))==1
%     warndlg("Cannot compare with an unique group");
%     return;
% end
end