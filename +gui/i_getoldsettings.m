function [para] = i_getoldsettings(src)
    para.SizeData = 10;
    if ~isprop(src, 'Parent') || ~isprop(src.Parent, 'Parent')
       error('Invalid source object: missing parent properties.');
    end
    ah = findobj(src.Parent.Parent, 'type', 'Axes');
    ha = findobj(ah.Children, 'type', 'Scatter');
    if ~isempty(ha)
        ha1 = ha(1);
        oldMarker = ha1.Marker;
        oldSizeData = ha1.SizeData;
        oldColorMap = colormap;
        para.oldMarker = oldMarker;
        para.oldSizeData = oldSizeData;
        para.oldColorMap = oldColorMap;
    end
end