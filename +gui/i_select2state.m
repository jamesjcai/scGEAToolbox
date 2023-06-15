function [thisc1,clable1,thisc2,clable2]=i_select2state(sce)

thisc1=[]; clable1='';
thisc2=[]; clable2='';

%function [thisc,clable,listitems,newpickclable]=i_select2state(sce)

%    thisc=[];
%    clable='';
%    newpickclable='';

    baselistitems={'Current Class (C)'};
    i_additem(sce.c_cluster_id, 'Cluster ID');
    i_additem(sce.c_cell_cycle_tx, 'Cell Cycle Phase');
    i_additem(sce.c_cell_type_tx, 'Cell Type');
    i_additem(sce.c_batch_id, 'Batch ID');
    
    function i_additem(itemv,itemn)
        if ~isempty(itemv) && length(unique(itemv))>=1
            baselistitems=[baselistitems,itemn];
        end
    end

    listitems=[baselistitems,...
        sce.list_cell_attributes(1:2:end)];
        nx=length(baselistitems);

%     a=evalin('base','whos');
%     b=struct2cell(a);
%     v=false(length(a),1);
%     for k=1:length(a)
%         if max(a(k).size)==sce.NumCells && min(a(k).size)==1
%             v(k)=true;
%         end
%     end
%     if any(v)
%         a=a(v);
%         b=b(:,v);
%         listitems=[listitems,'Customized C...'];
%     end

n=length(listitems);
if n<2
    warndlg(['This function requires at least two ' ...
        'grouping variables (e.g., BATCH_ID, ' ...
        'CLUSTER_ID, or CELL_TYPE_TXT).']);
    return;
end

% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select cell state/grouping variable:'},...
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



function [thisc,clable]=i_getidx(indx)
    clable=listitems{indx};
    switch clable
        case 'Library Size'
            thisc=sum(sce.X).';
        case 'Mt-reads Ratio'
            i = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
            lbsz = sum(sce.X, 1);
            lbsz_mt = sum(sce.X(i, :), 1);
            thisc = (lbsz_mt ./ lbsz).';        
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
        otherwise   % other properties
            nx=length(baselistitems);
            clable = sce.list_cell_attributes{2 * (indx - nx) - 1};
            thisc = sce.list_cell_attributes{2 * (indx - nx)};
            
    end
end


end