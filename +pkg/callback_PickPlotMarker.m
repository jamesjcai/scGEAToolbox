function callback_PickPlotMarker(src,~,h)    
    % h=gca;
    if h.Marker=='.'
        h.Marker='o';
        h.SizeData=10;
    else
        h.Marker='.';
        h.SizeData=20;
    end
end