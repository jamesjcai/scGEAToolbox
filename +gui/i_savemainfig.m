function i_savemainfig(src, ~, tag)
    
    
    if isa(src, 'matlab.apps.AppBase')    
        [parentfig] = xui.gui_getfigsce(src);
    else
        
        [parentfig] = gui.gui_getfigsce(src);
    end    
    axesHandles = findall(parentfig, 'Type', 'axes');
    if isempty(axesHandles)    
        gui.myHelpdlg(parentfig, ...
            'No figures available in the current window. Unable to save figure.');
        return;
    end
    if tag == 1
        filter = {'*.svg'};
        [filename, filepath] = uiputfile(filter);
        if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
        if ischar(filename)
            if gui.i_isuifig(parentfig)
                exportapp(parentfig, [filepath filename]);
            else
                saveas(parentfig, [filepath, filename], 'svg');
            end
        end
    elseif tag == 2
        % axx=gca;
        filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
        [filename, filepath] = uiputfile(filter);
        if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
        if ischar(filename)

            if gui.i_isuifig(parentfig)
                exportapp(parentfig, [filepath filename]);
            else
                exportgraphics(parentfig, [filepath, filename]);
            end
            
        end
    elseif tag == 3
        gui.i_export2pptx({parentfig}, {''}, parentfig);
    end
end

