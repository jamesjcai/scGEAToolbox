function [thisc,clable,listitems,newpickclable]=i_select1state(sce,nobaseline)

    if nargin<2, nobaseline=false; end
    thisc=[];
    clable='';
    newpickclable='';

    if nobaseline
        listitems=sce.list_cell_attributes(1:2:end);
    else
        baselistitems = {'Library Size', 'Mt-reads Ratio', ...
            'Cell Cycle Phase', ...
            'Cell Type', 'Cluster ID', 'Batch ID'};
        listitems=[baselistitems,...
            sce.list_cell_attributes(1:2:end)];
    end
    listitems=listitems(~cellfun(@isempty,listitems));
    
    if ~nobaseline
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
    end

if isempty(listitems), return; end

% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};
[indx2,tf2] = listdlg('PromptString',...
    {'Select cell state/grouping variable:'},...
     'SelectionMode','single','ListString',listitems);
if tf2==1
    clable=listitems{indx2};
    switch clable
        case 'Library Size'
            thisc=sum(sce.X).';
        case 'Mt-reads Ratio'
            i = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
            lbsz = sum(sce.X, 1);
            lbsz_mt = sum(sce.X(i, :), 1);
            thisc = (lbsz_mt ./ lbsz).';
        case 'Cluster ID' % cluster id            
            thisc=sce.c_cluster_id;
        case 'Batch ID' % batch id
            thisc=sce.c_batch_id;
        case 'Cell Type' % cell type
            thisc=sce.c_cell_type_tx;
        case 'Cell Cycle Phase' % cell cycle
            if isempty(sce.c_cell_cycle_tx) || all(strcmp(unique(sce.c_cell_cycle_tx),"undetermined"))
                fw = gui.gui_waitbar;
                sce = sce.estimatecellcycle(true,1);
                gui.gui_waitbar(fw);
            end
%             [~, tx] = grp2idx(sce.c_cell_cycle_tx);
%             ttxt = sprintf('%s|', string(tx));
%             clable = sprintf('%s (%s)',clable,ttxt);
            thisc=sce.c_cell_cycle_tx;
        case 'Customized C...'
            [thisc,newpickclable]=i_pickvariable;
        otherwise   % other properties
            [~,idx]=ismember(clable,sce.list_cell_attributes(1:2:end));
            thisc = sce.list_cell_attributes{idx+1};
            %nx=length(baselistitems);
            %clable = sce.list_cell_attributes{2 * (indx2 - nx) - 1};
            %thisc = sce.list_cell_attributes{2 * (indx2 - nx)};
    end
end


function [c,x]=i_pickvariable
    c=[]; x='';
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
            x = a(indx).name;
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

