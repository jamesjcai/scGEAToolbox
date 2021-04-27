function [thisc,clable]=i_select1class(sce)

thisc=[];
clable='';
listitems={'Current Class (C)','Cluster ID','Batch ID',...
           'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select Class','',''},...
     'SelectionMode','single','ListString',listitems);
if tf2==1
    switch indx2
        case 1
            thisc=sce.c;            
        case 2 % cluster id
            thisc=sce.c_cluster_id;
        case 3 % batch id
            thisc=sce.c_batch_id;
        case 4 % cell type
            thisc=sce.c_cell_type_tx;                
        case 5 % cell cycle
            thisc=sce.c_cell_cycle_tx;
    end
    clable=listitems{indx2};
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