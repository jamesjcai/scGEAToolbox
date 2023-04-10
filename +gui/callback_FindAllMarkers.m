function callback_FindAllMarkers(src,~)
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

answer = questdlg('Select Method',...
    '','Find All Markers','Marker Gene Heatmap',...
    'Find All Markers');
switch answer
    case 'Find All Markers'

    case 'Marker Gene Heatmap'
        gui.callback_MarkerGeneHeatmap(src);
        return;
    otherwise
        return;
end

    [thisc,~]=gui.i_select1class(sce);
    if isempty(thisc), return; end

    [T]=pkg.e_findallmarkers(sce.X,sce.g,thisc,[],[],[],true);
    if ~isempty(T)
        needwait=true;
        [answer,filename]=gui.i_exporttable(T,needwait);
        if ~isempty(answer)
            disp(filename);
            helpdlg(sprintf('All Markers Table saved.'),'');
        end
    else
            helpdlg('No results.','');
    end

