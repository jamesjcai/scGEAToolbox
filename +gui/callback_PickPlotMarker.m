function callback_PickPlotMarker(src,~)    
    ah=findobj(src.Parent.Parent,'type','Axes');
    h=findobj(ah.Children,'type','Scatter');
    if h.Marker=='.'
        h.Marker='o';
        h.SizeData=10;
    else
        h.Marker='.';
        h.SizeData=20;
    end
end