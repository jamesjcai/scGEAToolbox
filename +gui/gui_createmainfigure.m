function [parentfig,hAx] = gui_createmainfigure(v1,useuifigure)
    if nargin<2, useuifigure = false; end

    if ~isempty(v1)
        figname = sprintf('SCGEATOOL v%s', v1);
    else
        figname = 'SCGEATOOL';
    end

    if useuifigure
        a=round(1.2*[560, 420]);
        parentfig = uifigure('Name', figname, ...
                'position', [0, 0, a], ...
                'visible', 'off');
        %p = uipanel(parentfig,'Position',[0 0 a-5]);
        hAx = uiaxes(parentfig,'Position',[65 50 a(1)-120 a(2)-60], ...
            'Visible','off');
    else
        defaultPosition = get(groot, 'DefaultFigurePosition');
        defaultWidth = defaultPosition(3);
        defaultHeight = defaultPosition(4);
        
        if defaultWidth==560 && defaultHeight==420
            parentfig = figure('Name', figname, ...
                'position', round(1.2*[0, 0, 560, 420]), ...
                'visible', 'off', 'NumberTitle', 'off', ...
                'DockControls','off','MenuBar','none','ToolBar','Figure');
        else
            parentfig = figure('Name', figname, ...    
                'visible', 'off', 'NumberTitle', 'off', ...
                'DockControls','off','MenuBar','none','ToolBar','Figure');
        end
        delete(findall(parentfig, 'Tag', 'FigureToolBar'));        
        dt = datacursormode(parentfig);
        dt.Enable = 'off';
        dt.UpdateFcn = {@i_myupdatefcnx};
        hAx = axes('Parent', parentfig, 'Visible', 'off');
    end
    movegui(parentfig, 'center');

    function [txt] = i_myupdatefcnx(pdt, ~)
        % pos = event_obj.Position;
        % idx = event_obj.DataIndex;
        % txt = cL(c(idx));
        % https://www.mathworks.com/matlabcentral/answers/549567-disable-datatips-on-click
        pdt.Visible = 'off';
        txt = '';
    end
end