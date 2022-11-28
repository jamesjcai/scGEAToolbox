function callback_MultiGroupingViewer(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
        
    [thisc1,clable1,thisc2,clable2]=gui.i_select2state(sce);
     if isempty(thisc1) || isempty(thisc2)
         return;
     end

    [c,cL] = grp2idx(thisc1);
    cx1.c=c; cx1.cL=strrep(cL,'_','\_');
    [c,cL] = grp2idx(thisc2);
    cx2.c=c; cx2.cL=strrep(cL,'_','\_');
    gui.sc_multigroupings(sce,cx1,cx2,clable1,clable2);
end
