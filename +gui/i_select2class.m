function [thisc1,clable1,thisc2,clable2]=i_select2class(sce)

thisc1=[]; clable1='';
thisc2=[]; clable2='';

    listitems={'Current Class (C)'};
    i_additem(sce.c_cluster_id, 'Cluster ID');
    i_additem(sce.c_cell_cycle_tx, 'Cell Cycle Phase');
    i_additem(sce.c_cell_type_tx, 'Cell Type');
    i_additem(sce.c_batch_id, 'Batch ID');
    
    function i_additem(itemv,itemn)
    if ~isempty(itemv)&&length(unique(itemv))>1
        listitems=[listitems,itemn];
    end
    end
        
        
    
    
n=length(listitems);
if n<3
    errordlg('Need at least two grouping variables (e.g., BATCH_ID, CLUSTER_ID, or CELL_TYPE_TXT)');
    return;
end

% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select two grouping varible:'},...
     'SelectionMode','multiple',...
     'ListString',listitems,...
     'InitialValue',[n-1 n]);
if tf2==1
    if length(indx2)~=2
        warndlg('Please select 2 grouping variables.');
        return;
    end
    [thisc1,clable1]=i_getidx(indx2(1));
    [thisc2,clable2]=i_getidx(indx2(2));
end

% if isempty(thisc)
%     errordlg('Undefined');    
%     return;
% end
% if numel(unique(thisc))==1
%     warndlg("Cannot compare with an unique group");
%     return;
% end



function [thisc,clable]=i_getidx(indx)
    clable=listitems{indx};
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

end