function callback_ShowClustersPop(src, ~)
    answer = questdlg('Select a grouping variable and show cell groups in new figures individually?');
    if ~strcmp(answer, 'Yes'), return; end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    [thisc,~]=gui.i_select1class(sce);
    if isempty(thisc), return; end

    try

    [c,cL]=grp2idx(thisc);
    cLa=getappdata(FigureHandle,'cL');
    if ~isempty(cLa) && length(cL)==length(cLa)
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

    totaln=max(c);
    numfig=ceil(totaln/9);
    for nf=1:numfig
        f=figure('visible','off');
        for k = 1:9        
            kk=(nf-1)*9+k;
            if kk<=totaln
                subplot(3, 3, k);
                gui.i_gscatter3(sces, c, 3, cmv(idxx(kk)));
                title(sprintf('%s\n%d cells (%.2f%%)', ...
                    cL{idxx(kk)}, cmx(idxx(kk)), ...
                    100 * cmx(idxx(kk)) / length(c)));
            end
            colormap(para.oldColorMap);
        end
        P = get(f,'Position');
        set(f,'Position',[P(1)-20*nf P(2)-20*nf P(3) P(4)]);
        set(f,'visible','on');
        drawnow;
    end
    catch ME
        errordlg(ME.message);
    end
end
