function uicallback_Violinplot(src, ~)

if isa(src,"SingleCellExperiment")
    FigureHandle = gcf;
    sce = src;
else
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
end





[thisc, ~] = gui.i_select1class(sce);
if isempty(thisc), return; end


[c, cL]=grp2idx(thisc);
cL = strrep(cL, '_', '\_');
% [c, cL, noanswer] = gui.i_reordergroups(thisc);
% if noanswer, return; end

[glist] = gui.i_selectngenes(sce);
if isempty(glist)
    helpdlg('No gene selected.', '');
    return;
end
[Xt] = gui.i_transformx(sce.X);
% glist=glist(end:-1:1);


tabgp = uitabgroup();
idx = 1;
focalg = glist(idx);
thisc = strrep(string(cL(c)), '_', '\_');
for k = 1:length(glist)
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',glist(k)));
    y = full(Xt(upper(sce.g) == upper(glist(k)), :));
    pkg.i_violinplot(y, thisc, true, cL);
    title(strrep(glist(k), '_', '\_'));
end
if isempty(px_new)
    movegui(hFig,'center');
else    
    movegui(hFig, px_new);
end
drawnow;
hFig.Visible=true;
return;

    try
        %f=gui.i_violinmatrix(Xt,sce.g,c,cL,"PRLR");
        for k = 1:length(glist)
            y = Xt(upper(sce.g) == upper(glist(k)), :);
            [f] = gui.i_violinplot(full(y), cL(c), glist(k), true, cL);
            % f=gui.i_violinplot(sce.X,sce.g,c,cL,glist);
            p = f.Position;
            %p(2)=0;
            p(4) = p(4) * 1.5;
            f.Position = p;
            f.Visible = "on";
            % pause(1);
            drawnow;
        end    
    catch ME
        if exist('f','var') && ishandle(f)
            close(f);
        end
        errordlg(ME.message);
    end
end
