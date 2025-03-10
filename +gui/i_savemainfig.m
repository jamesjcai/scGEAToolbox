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
            saveas(FigureHandle, [filepath, filename], 'svg');
        end
    elseif tag == 2
        % axx=gca;
        filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
        [filename, filepath] = uiputfile(filter);
        if ischar(filename)
            exportgraphics(FigureHandle, [filepath, filename]);
        end
    elseif tag == 3
        gui.i_export2pptx({FigureHandle}, {''});
    end
end

