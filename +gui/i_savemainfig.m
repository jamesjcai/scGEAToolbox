function i_savemainfig(src, ~, tag)
    
    [FigureHandle, ~, isui] = gui.gui_getfigsce(src);
    
    axesHandles = findall(FigureHandle, 'Type', 'axes');
    if isempty(axesHandles)    
        gui.myHelpdlg(FigureHandle, 'No figures available in the current window. Unable to save figure.', '');
        return;
    end
    if tag == 1
        filter = {'*.svg'};
        [filename, filepath] = uiputfile(filter);
        if ischar(filename)
            if gui.i_isuifig(FigureHandle)
                exportapp(FigureHandle, [filepath filename]);
            else
                saveas(FigureHandle, [filepath, filename], 'svg');
            end
        end
    elseif tag == 2
        % axx=gca;
        filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
        [filename, filepath] = uiputfile(filter);
        if ischar(filename)

            if gui.i_isuifig(FigureHandle)
                exportapp(FigureHandle, [filepath filename]);
            else
                exportgraphics(FigureHandle, [filepath, filename]);
            end
            
        end
    elseif tag == 3
        gui.i_export2pptx({FigureHandle}, {''}, FigureHandle);
    end
end

