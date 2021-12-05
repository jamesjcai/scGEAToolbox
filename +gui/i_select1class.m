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
    
    a=evalin('base','whos');
    b=struct2cell(a);
    v=false(length(a),1);
    for k=1:length(a)
        if max(a(k).size)==sce.NumCells && min(a(k).size)==1
            v(k)=true;
        end
    end
    if any(v)
        a=a(v);
        b=b(:,v);
        listitems=[listitems,'Customized C...'];
    end


% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select one grouping variable:'},...
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
        case 'Customized C...'
            thisc=i_pickvariable;
    end
end


function [c]=i_pickvariable
    c=[];
%     a=evalin('base','whos');
%     b=struct2cell(a);
%     v=false(length(a),1);
%     for k=1:length(a)
%         if max(a(k).size)==sce.NumCells && min(a(k).size)==1
%             v(k)=true;
%         end
%     end
%     if any(v)
        %valididx=ismember(b(4,:),'double');
        %a=a(valididx);
        [indx,tf]=listdlg('PromptString',{'Select network variable:'},...
            'liststring',b(1,:),'SelectionMode','single');
        if tf==1
            c = evalin('base',a(indx).name);
        end
%    end
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

