function callback_GeneHeatMap(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [thisc,~]=gui.i_select1class(sce);
    if isempty(thisc), return; end
    [c,cL] = grp2idx(thisc);
    [answer]=questdlg('Manually order groups?','');
    switch answer
        case 'Yes'
            [newidx]=gui.i_selmultidlg(cL);
            if length(newidx)~=length(cL)
                return;
            end
            cx=c;
            for k=1:length(newidx)
                c(cx==newidx(k))=k;
            end
            cL=cL(newidx);
        otherwise
    end

    [glist]=gui.i_selectngenes(sce);
    if isempty(glist)
        helpdlg('No gene selected.','');
        return;
    end


[~,gidx]=ismember(glist,sce.g);

Y=sce.X(gidx,:);
[~,cidx]=sort(c);
Y=Y(:,cidx);

szgn=grpstats(c,c,@numel);


a=zeros(1,max(c)); 
b=zeros(1,max(c));
for k=1:max(c)
    a(k)=sum(c<=k);
    b(k)=round(sum(c==k)./2);
end

figure;
imagesc(Y);
set(gca,'XTick',a-b);
set(gca,'XTickLabel',cL);
set(gca,'XTickLabelRotation',45);
set(gca,'YTick',1:length(glist));
set(gca,'YTickLabel',glist);
set(gca,'TickLength',[0 0]);
szc=cumsum(szgn);
for k=1:length(szc)
    xline(szc(k)+0.5,'y-');
end

end
