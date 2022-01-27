function callback_PickPlotMarker(src,~)    
    ah=findobj(src.Parent.Parent,'type','Axes');
    ha=findobj(ah.Children,'type','Scatter');
    
    for k=1:length(ha)
    h=ha(k);
    if h.Marker=='.'
        h.Marker='o';
        h.SizeData=10*randi(10);
    else
        h.Marker='.';
        h.SizeData=50*randi(10);
    end    
    end
end