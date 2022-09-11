function callback_GeneHeatMap(src,~)
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


[~,gidx]=ismember(glist,sce.g);

%[Xt]=gui.i_transformx(sce.X);
Xt=sc_norm(sce.X);
Xt=log(Xt+1);

Y=Xt(gidx,:);
[~,cidx]=sort(c);
Y=Y(:,cidx);

methodid=1;
switch methodid
    case 1
        Y=zscore(Y,0,2);
        qx=quantile(Y(:),0.90);
        Y(Y>qx)=qx;
        qx=quantile(Y(:),0.10);
        Y(Y<qx)=qx;
    case 2
        Y=zscore(Y,0,2);
        Y = reshape( zscore(Y(:)), size(Y) );
        Y=Y./(max(abs(Y(:))));
end

szgn=grpstats(c,c,@numel);
a=zeros(1,max(c)); 
b=zeros(1,max(c));
for k=1:max(c)
    a(k)=sum(c<=k);
    b(k)=round(sum(c==k)./2);
end

% figure;
% heatmap(Y)
% assignin('base','Y',Y);
% assignin('base','g',glist) ;
% heatmap(Y,'YDisplayLabels',glist, ...
%     'XDisplayLabels',strings(size(Y,2),1), ...
%     'GridVisible',false,'ColorScaling','scaled',...
%     'ColorbarVisible',false)

hFig=figure;
imagesc(Y);
% hFig.Colormap = repmat(linspace(0, 1, 25).', 1, 3);
set(gca,'XTick',a-b);
set(gca,'XTickLabel',cL);
set(gca,'XTickLabelRotation',45);
set(gca,'YTick',1:length(glist));
set(gca,'YTickLabel',glist);
set(gca,'TickLength',[0 0]);
colormap(flipud(bone));
box on

szc=cumsum(szgn);
for k=1:length(szc)
    xline(szc(k)+0.5,'y-');
end
tb = uitoolbar('Parent', hFig);
pkg.i_addbutton2fig(tb,'on',{@gui.i_pickcolormap,c},'plotpicker-compass.gif','Pick new color map...');
pkg.i_addbutton2fig(tb,'off',@gui.i_changefontsize,'noun_font_size_591141.gif','ChangeFontSize');
pkg.i_addbutton2fig(tb,'on',@i_renamecat,'guideicon.gif','Rename groups...');
pkg.i_addbutton2fig(tb,'on',{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb,'on',@i_resetcolor,'plotpicker-geobubble2.gif','Reset color map');

    function i_renamecat(~,~)
        tg=gui.i_inputgenelist(string(cL),true);
        if length(tg)==length(cL)
            set(gca,'XTick',a-b);
            set(gca,'XTickLabel',tg(:))
            cL=tg;
        else
            errordlg('Wrong input.');
        end
    end

    function i_resetcolor(~,~)
        set(gca,'FontSize',10);
        if rand>0.5
            colormap(flipud(bone));
        else
            colormap(bone);
        end
    end

    end
