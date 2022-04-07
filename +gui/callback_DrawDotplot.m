function callback_DrawDotplot(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [thisc,~]=gui.i_select1class(sce);
    if ~isempty(thisc)
        [c,cL] = grp2idx(thisc);
    else
        return;
    end
    [glist]=gui.i_selectngenes(sce);
    if isempty(glist)
        helpdlg('No gene selected.','');
        return;
    end
    try        
        f=gui.i_dotplot(sce.X,sce.g,c,cL,glist);
    catch ME
        if exist('f','var') && ishandle(f)
            close(f);
        end
        errordlg(ME.message);
    end
end