function [para] = i_getoldsettings(src)
ah = findobj(src.Parent.Parent, 'type', 'Axes');
ha = findobj(ah.Children, 'type', 'Scatter');
ha1 = ha(1);
oldMarker = ha1.Marker;
oldSizeData = ha1.SizeData;
oldColorMap = colormap;
para.oldMarker = oldMarker;
para.oldSizeData = oldSizeData;
para.oldColorMap = oldColorMap;
end