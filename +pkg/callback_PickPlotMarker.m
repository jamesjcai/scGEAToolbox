function callback_PickPlotMarker(~,~,h)
    if h.Marker=='.'
        h.Marker='o';
        h.SizeData=10;
    else
        h.Marker='.';
        h.SizeData=20;
    end
end