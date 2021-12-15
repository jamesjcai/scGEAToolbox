function callback_ShowClustersPop(src, ~)
    answer = questdlg('Show clusters in new figures?');
    if ~strcmp(answer, 'Yes')
        return
    end
FigureHandle=src.Parent.Parent;
sce=guidata(FigureHandle);

[c,cL]=grp2idx(sce.c);
cLa=getappdata(FigureHandle,'cL');
if ~isempty(cLa)
    cL=cLa;
end

    cmv = 1:max(c);
    idxx = cmv;
    [cmx] = countmember(cmv, c);
    answer = questdlg('Sort by size of cell groups?');
    if strcmpi(answer, 'Yes')
        [~, idxx] = sort(cmx, 'descend');
    end
    sces = sce.s;
    h=findall(FigureHandle,'type','scatter');
    if isempty(h.ZData)
        sces = sce.s(:, 1:2);
    end

    [para] = gui.i_getoldsettings(src);
    f1=figure;
    for k = 1:9
        if k <= max(c)
            subplot(3, 3, k);
            gui.i_gscatter3(sces, c, 3, cmv(idxx(k)));
            title(sprintf('%s\n%d cells (%.2f%%)', ...
                cL{idxx(k)}, cmx(idxx(k)), ...
                100 * cmx(idxx(k)) / length(c)));
        end
        colormap(para.oldColorMap);
    end
        P = get(f1,'Position');
        
    if ceil(max(c) / 9) == 2
        f2=figure;
        for k = 1:9
            kk = k + 9;
            if kk <= max(c)
                subplot(3, 3, k);
                gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
                title(sprintf('%s\n%d cells (%.2f%%)', ...
                    cL{idxx(kk)}, cmx(idxx(kk)), ...
                    100 * cmx(idxx(kk)) / length(c)));
            end
        end
        colormap(para.oldColorMap);
        set(f2,'Position',[P(1)-10*k P(2)-10*k P(3) P(4)]);                
    end
    
    if ceil(max(c) / 9) > 2
        warndlg('Group(s) #18 and above are not displayed');
    end
end


%     function [para] = i_getoldsettings(src)
%         ah = findobj(src.Parent.Parent, 'type', 'Axes');
%         ha = findobj(ah.Children, 'type', 'Scatter');
%         ha1 = ha(1);
%         oldMarker = ha1.Marker;
%         oldSizeData = ha1.SizeData;
%         oldColorMap = colormap;
%         para.oldMarker = oldMarker;
%         para.oldSizeData = oldSizeData;
%         para.oldColorMap = oldColorMap;
%     end