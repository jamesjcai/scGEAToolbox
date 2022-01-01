function callback_PickPlotMarker(src,~)    
    ah=findobj(src.Parent.Parent,'type','Axes');
    ha=findobj(ah.Children,'type','Scatter');
    
    for k=1:length(ha)
    h=ha(k);
    if h.Marker=='.'
        h.Marker='o';
    else
        h.Marker='.';
    end
    h.SizeData=10*randi(10);
    end
end