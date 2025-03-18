function [requirerefresh, renamedwhat] = callback_RenameOthers(src, ~)

requirerefresh = false;
renamedwhat = [];

[FigureHandle, sce] = gui.gui_getfigsce(src);

listitems = {'Gene name','Cluster ID'};
    
       if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, listitems, ...
                'Select varible to be renamed:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select varible to be renamed:'}, ...
                'SelectionMode', 'single', ...
                'ListString', listitems, ...
                'ListSize',[220 300]);
        end

    if tf2 ~= 1, return; end    

    switch listitems{indx2}
        case 'Gene name'
            [requirerefresh] = gui.callback_RenameGenes(src);
            if requirerefresh, renamedwhat = listitems{indx2}; end

        case 'Cluster ID'
            if isempty(sce.c_cluster_id)
                sce.c_cluster_id = string(ones(sce.NumCells, 1));
            end
            if ~isstring(sce.c_cluster_id)
                sce.c_cluster_id = string(sce.c_cluster_id);
            end
            [ci, cLi] = pkg.i_grp2idxsorted(sce.c_cluster_id);
            
           if gui.i_isuifig(FigureHandle)
                [indxx, tfx] = gui.myListdlg(FigureHandle, string(cLi), ...
                    'Select cluster ID');
            else
                [indxx, tfx] = listdlg('PromptString', ...
                    {'Select cluster ID'}, ...
                    'SelectionMode', 'single', ...
                    'ListString', string(cLi), 'ListSize', [220, 300]);
            end

            if tfx == 1
                i = ismember(ci, indxx);
                
                if gui.i_isuifig(FigureHandle)
                    newctype = gui.myInputdlg({'New cluster ID name'}, 'Rename', ...
                        {cLi(ci(i))}, FigureHandle);
                else
                    newctype = inputdlg('New cluster ID name', 'Rename', [1, 50], cLi(ci(i)));
                end

                if ~isempty(newctype)
                    cLi(ci(i)) = newctype;
                    sce.c_cluster_id = string(cLi(ci));
                    requirerefresh = true;                    
                    renamedwhat = listitems{indx2};
                    guidata(FigureHandle, sce);
                end
            end        
    end

end
