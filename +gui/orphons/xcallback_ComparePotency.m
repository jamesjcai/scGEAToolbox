function xcallback_ComparePotency(src, ~)
% UNUSED function will be removed.

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
[a] = contains(sce.list_cell_attributes(1:2:end), 'cell_potency');
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
            answer2 = questdlg('Which species?', 'Select Species', 'Mouse', 'Human', 'Mouse');
            [y, specisid] = ismember(lower(answer2), {'human', 'mouse'});
            if ~y, return; end
            fw = gui.gui_waitbar;
            try
                sce = sce.estimatepotency(specisid);

                [y, idx] = ismember({'cell_potency'}, sce.list_cell_attributes(1:2:end));
                if y
                    c = sce.list_cell_attributes{idx+1};
                    figure;
                    gui.i_gscatter3(sce.s, c);
                end
                guidata(FigureHandle, sce);
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

[yes, idx] = ismember('cell_potency', sce.list_cell_attributes(1:2:end));
if ~yes
    warndlg('Check SCE.LIST_CELL_ATTRIBUTES(1:2)')
    return;
end
[thisc, clabel] = gui.i_select1class(sce);
if isempty(thisc) % || numel(unique(thisc))==1
    errordlg('Undefined');
    return;
end
x = sce.list_cell_attributes{idx+1};
figure;
pkg.i_violinplot_groupordered(x, thisc);
xlabel(clabel)
end
