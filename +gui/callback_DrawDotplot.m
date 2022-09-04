function callback_DrawDotplot(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [thisc,~]=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [c,cL] = grp2idx(thisc);
    [answer]=questdlg('Manually order groups?','', ...
        'Yes','No','Cancel','No');
    if isempty(answer), return; end
    switch answer
        case 'Yes'
            [newidx]=gui.i_selmultidlg(cL);
            if length(newidx)~=length(cL)
                warndlg('Please select all group items.','');
                return;
            end
            cx=c;
            for k=1:length(newidx)
                c(cx==newidx(k))=k;
            end
            cL=cL(newidx);
        case 'No'
            
        case 'Cancel'
            return;
        otherwise
            return;
    end

    [glist]=gui.i_selectngenes(sce);
    if isempty(glist)
        helpdlg('No gene selected.','');
        return;
    end
    [Xt]=gui.i_transformx(sce.X);
    try
        f=gui.i_dotplot(Xt,sce.g,c,cL,glist,true);
        % f=gui.i_violinplot(sce.X,sce.g,c,cL,glist);
    catch ME
        if exist('f','var') && ishandle(f)
            close(f);
        end
        errordlg(ME.message);
    end
end