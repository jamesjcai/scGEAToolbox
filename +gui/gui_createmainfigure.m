function [FigureHandle] = gui_createmainfigure(v1)
    if ~isempty(v1)
        figname = sprintf('SCGEATOOL v%s', v1);
    else
        figname = 'SCGEATOOL';
    end        
    defaultPosition = get(groot, 'DefaultFigurePosition');
    defaultWidth = defaultPosition(3);
    defaultHeight = defaultPosition(4);        
    if defaultWidth==560 && defaultHeight==420
        FigureHandle = figure('Name', figname, ...
            'position', round(1.2*[0, 0, 560, 420]), ...
            'visible', 'off', 'NumberTitle', 'off', ...
            'DockControls','off','MenuBar','none','ToolBar','Figure');
    else
        FigureHandle = figure('Name', figname, ...    
            'visible', 'off', 'NumberTitle', 'off', ...
            'DockControls','off','MenuBar','none','ToolBar','Figure');
    end
    delete(findall(FigureHandle, 'Tag', 'FigureToolBar'));
    movegui(FigureHandle, 'center');
    dt = datacursormode(FigureHandle);
    dt.Enable = 'off';
    dt.UpdateFcn = {@i_myupdatefcnx};

    function [txt] = i_myupdatefcnx(pdt, ~)
        % pos = event_obj.Position;
        % idx = event_obj.DataIndex;
        % txt = cL(c(idx));
        % https://www.mathworks.com/matlabcentral/answers/549567-disable-datatips-on-click
        pdt.Visible = 'off';
        txt = '';
    end    
end