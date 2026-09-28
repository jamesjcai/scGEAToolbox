function callback_ShowHgBGeneExpression(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);

species = gui.i_selectspecies(2, false, FigureHandle);
if isempty(species), return; end

switch species
    case 'mouse'
        idx1 = startsWith(sce.g, 'Hba-', 'IgnoreCase', true);
        idx2 = startsWith(sce.g, 'Hbb-', 'IgnoreCase', true);
        idx3 = false;
    case 'human'
        idx1 = strcmpi(sce.g, "HBA1");
        idx2 = strcmpi(sce.g, "HBA2");
        idx3 = strcmpi(sce.g, "HBB");
end
idx = idx1 | idx2 | idx3;

if any(idx)
    ttxt = sprintf("%s+", sce.g(idx));
    ci = full(sum(sce.X(idx, :), 1));

    hx = gui.myFigure(FigureHandle);
    hFig = hx.FigHandle;

    cm = uicontextmenu(hFig);
    m1 = uimenu(cm, 'Text', 'Save HgBGeneExpression...', "MenuSelectedFcn", {@in_callback_saveM, ci});
    hFig.ContextMenu = cm;

    gui.i_stemscatter(sce.s, ci, hx.AxHandle);

    title(hx.AxHandle, ttxt);
    hx.addCustomButton('off', {@in_callback_saveM, ci}, 'floppy-disk-arrow-in.jpg', ...
        'Save marker gene map...');
    hx.show(FigureHandle);
else
    gui.myWarndlg(FigureHandle, 'No Hgb-genes found');
end


function in_callback_saveM(~, ~, ~)
        if ~(ismcc || isdeployed)
            labels = {'Save C_CELL_ID to variable named:', ...
                'Save HgBGeneExpression to variable named:'};
            vars = {'cell_id', 'c'};
            values = {sce.c_cell_id, ci(:)};
            waitfor(export2wsdlg(labels, vars, values));   % one at a time
        else
            gui.myErrordlg(hx.FigHandle, ['This function is not available for standalone application.' ...
                ' Run scgeatool in MATLAB to use this function.']);
        end
    end
end
