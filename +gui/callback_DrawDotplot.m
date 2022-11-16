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
            [newidx]=gui.i_selmultidlg(cL,sort(cL));
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
    glist=glist(end:-1:1);

if length(glist)>75
    for k=1:50:length(glist)
        k2=min([length(glist),k+50-1]);
        f=gui.i_dotplot(Xt,sce.g,c,cL,glist(k:k2),true);
    end
else

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

end