function callback_PickPlotMarker(src,~)    
    ah=findobj(src.Parent.Parent,'type','Axes');
    h=findobj(ah.Children,'type','Scatter');
    if h.Marker=='.'
        h.Marker='o';
        if rand>0.5
            h.SizeData=20;
        else
            h.SizeData=10;
        end
    else
        h.Marker='.';
        if rand>0.5
            h.SizeData=20;
        else
            h.SizeData=10;
        end
    end
end