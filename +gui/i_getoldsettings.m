function [para] = i_getoldsettings(src, parentfig)

if isa(src, 'matlab.apps.AppBase')
    para.SizeData = 10;
    ha1 = src.h;
    if ~isempty(ha1)

        oldMarker = ha1.Marker;
        oldSizeData = ha1.SizeData;
        % The axes, not the figure: i_gscatter3 sets the map on the axes,
        % and that does not change the figure's, so reading the figure
        % handed back the default map and a 3D-to-2D projection redrew the
        % groups in parula.
        oldColorMap = colormap(src.UIAxes);
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
