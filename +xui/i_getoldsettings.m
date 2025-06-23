function [para] = i_getoldsettings(app)

    para.SizeData = 10;
    ha1 = app.h;
    if ~isempty(ha1)
        
        oldMarker = ha1.Marker;
        oldSizeData = ha1.SizeData;
        oldColorMap = colormap(app.UIFigure);
        para.oldMarker = oldMarker;
        para.oldSizeData = oldSizeData;
        para.oldColorMap = oldColorMap;
    end
end