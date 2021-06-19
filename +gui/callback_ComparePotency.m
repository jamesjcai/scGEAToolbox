function callback_ComparePotency(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [a]=contains(sce.list_cell_attributes(1:2:end),'cell_potency');
    if any(a)
                answer2 = questdlg(sprintf('CELL_POTENCY has been computed.\nCompare it across cell classes?'));
                switch answer2
                    case 'Yes'
                    otherwise
                        return;
                end
        
    else
        answer = questdlg('Compute cell differentiation potency (CELL_POTENCY))?');
        switch answer
            case 'Yes'
                answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                [y,specisid]=ismember(lower(answer2),{'human','mouse'});
                if ~y, return; end
                fw=gui.gui_waitbar;
                try
                sce=sce.estimatepotency(specisid);
                guidata(FigureHandle,sce);
                gui.gui_waitbar(fw);
                catch ME
                    gui.gui_waitbar(fw);
                    errordlg(ME.message)
                    %rethrow(ME)
                end
                answer2 = questdlg(sprintf('CELL_POTENCY is computed.\nContinue to compare it across cell classes?'));
                switch answer2
                    case 'Yes'
                    otherwise
                        return;
                end
            otherwise
                return;
        end
    end
    
    [yes,idx]=ismember('cell_potency',sce.list_cell_attributes(1:2:end));
    if ~yes
        warndlg('Check SCE.LIST_CELL_ATTRIBUTES(1:2)')
        return;
    end
    [thisc,clabel]=gui.i_select1class(sce);
    if isempty(thisc)   % || numel(unique(thisc))==1
        errordlg('Undefined');    
        return;
    end
    x=sce.list_cell_attributes{idx+1};    
    figure;
    pkg.i_violinplot_groupordered(x,thisc);
    xlabel(clabel)
    
    
    % answer = questdlg('Compute each cell''s differentiation potency?');
    %{
    if ~isempty(sce.list_cell_attributes{idx+1}) && ~isempty(sce.c_cell_type_tx)
        x=sce.list_cell_attributes{idx+1};
        y=sce.c_cell_type_tx;
        figure;
        pkg.i_violinplot_groupordered(x,y);
        xlabel('Cell Type')
    end
    
    if ~isempty(sce.list_cell_attributes{idx+1}) && ~isempty(sce.c_batch_id) && length(unique(sce.c_batch_id))>1
        x=sce.list_cell_attributes{idx+1};
        y=sce.c_batch_id;
        figure;
        pkg.i_violinplot_groupordered(x,y);
        xlabel('Batch ID')
    end
    %}
        
%  %       try
%             x=sce.list_cell_attributes{idx+1};
%             f=figure('Visible','Off');
%             h1=subplot(1,2,1);
%             i=sce.c_batch_id==1;
%             pkg.i_violinplot_groupordered(x(i),y(i));
%             [aabb]=ylim();
%             h2=subplot(1,2,2);
%             i=sce.c_batch_id==2;
%             pkg.i_violinplot_groupordered(x(i),y(i));
%             [ccdd]=ylim();
%             ylim(h1,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
%             ylim(h2,[min(aabb(1),ccdd(1)) max(aabb(2),ccdd(2))]);
%             set(f,'Visible','On');
% %         catch
% %             disp('Check SCE.C_BATCH_ID');
% %         end
    
end
