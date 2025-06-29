function [para] = i_getoldsettings(src, parentfig)

    if isa(src, 'matlab.apps.AppBase')  
        para.SizeData = 10;
        ha1 = src.h;
        if ~isempty(ha1)
            
            oldMarker = ha1.Marker;
            oldSizeData = ha1.SizeData;
            oldColorMap = colormap(src.UIFigure);
            para.oldMarker = oldMarker;
            para.oldSizeData = oldSizeData;
            para.oldColorMap = oldColorMap;
        end
    
    else
    
        para.SizeData = 10;
        % if ~isprop(src, 'Parent') || ~isprop(src.Parent, 'Parent')
        %    error('Invalid source object: missing parent properties.');
        % end
        % src - is a PushTool (button handle)
        % src.Parent - is a Toolbar
        % src.Parent.Parent - is the figure.
    
        if nargin<2, parentfig = []; end
        if isempty(parentfig)
            parentfig = src.Parent.Parent;
        end
        ah = findobj(parentfig, 'type', 'Axes');
    
        % assignin('base',"ah",ah)
        ha = findobj(ah.Children, 'type', 'Scatter');
        if ~isempty(ha)
            ha1 = ha(1);
            oldMarker = ha1.Marker;
            oldSizeData = ha1.SizeData;
            oldColorMap = colormap(ah);
            para.oldMarker = oldMarker;
            para.oldSizeData = oldSizeData;
            para.oldColorMap = oldColorMap;
        end
    end

end