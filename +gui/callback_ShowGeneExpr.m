function callback_ShowGeneExpr(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end

    % answer = questdlg("Select the type of expression values","",...
    %     "Raw UMI Counts","Library Size-Normalized",)

    [Xt] = gui.i_transformx(sce.X);
    if isempty(Xt), return; end

    n = length(glist);
    an1 = "Yes, same figure";
    an2 = "No, different tabs";
    if n > 1
        answer = questdlg("Plot on all genes in the same figure?", "", an1, ...
            an2, "Cancel", an1);
        if isempty(answer), return; end
        if strcmp(answer, "Cancel"), return; end
    else
        answer = an2;
    end
    
    switch answer
        case an1
            answer2 = questdlg("Type of plot:","", "stem plot", "feature plot", "stem plot");
            if isempty(answer2), return; end
            fw = gui.gui_waitbar; 
            hx = gui.myFigure;
            maxy = 0;
            if ~ispref('scgeatoolbox', 'prefcolormapname')
                setpref('scgeatoolbox', 'prefcolormapname', 'autumn');
            end
            a = getpref('scgeatoolbox', 'prefcolormapname');
            for k = 1:n
                nexttile
                switch answer2
                    case "feature plot"
                        sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
                        c = Xt(sce.g == glist(k), :);
                        gui.i_setautumncolor(c, a, true, any(c==0));
                        colorbar;
                    case "stem plot"
                        sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);                        
                end
                maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
            end
            gui.gui_waitbar(fw);
            hx.show(FigureHandle);
        case an2
            fw = gui.gui_waitbar;
            y = cell(n,1);
            for k = 1:n
                y{k} = Xt(sce.g == glist(k), :);
            end
            gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
            gui.gui_waitbar(fw);            
    end
end
