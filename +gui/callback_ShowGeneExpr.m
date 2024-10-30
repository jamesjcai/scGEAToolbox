function callback_ShowGeneExpr(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end

    % answer = questdlg("Select the type of expression values","",...
    %     "Raw UMI Counts","Library Size-Normalized",)

    [Xt] = gui.i_transformx(sce.X);
    if isempty(Xt), return; end

    n = length(glist);
    if n > 1
        answer = questdlg("Plot on the same figure?","", "Yes, same figure", ...
            "No, different figures", "Cancel", "Yes, same figure");
        if isempty(answer), return; end
        if strcmp(answer, "Cancel"), return; end
    else
        answer = "No, different figures";
    end
    
    switch answer
        case "No, different figures"
            y = cell(n,1);
            for k = 1:n
                y{k} = Xt(sce.g == glist(k), :);
            end
            gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
        case "Yes, same figure"
            answer2 = questdlg("Type of plot:","", "feature plot", "stem plot", "feature plot");
            if isempty(answer2), return; end
            hFig = figure(Visible="off");
            maxy = 0;
            for k = 1:n
                nexttile
                switch answer2
                    case "feature plot"
                        sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
                        colorbar;
                    case "stem plot"
                        sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);
                end
                maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
            end
            gui.i_movegui2parent(hFig, FigureHandle);
            
            evalin('base', 'linkprop(findobj(gcf,''type'',''axes''), {''CameraPosition'',''CameraUpVector''});');

            % pkg.i_addbutton2fig(tb, 'off', {@in_rescale, maxy}, 'networkcomp.gif', 'Normalize scales...');
            % tb = uitoolbar('Parent', hFig);
            tb = findall(hFig, 'Type', 'uitoolbar', 'Tag', 'FigureToolBar');
            if ~isempty(tb)
                uipushtool(tb, 'Separator', 'off');
                pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');            
                pkg.i_addbutton2fig(tb, 'off', @gui.i_setboxon, 'RectGate.gif', 'Box on/off');            
                pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
                pkg.i_addbutton2fig(tb, 'off', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');            
                pkg.i_addbutton2fig(tb, 'off', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');
                gui.gui_3dcamera(tb, 'AllCells');
            end
            hFig.Visible = "on";
    end

    %    gui.gui_waitbar(fw);

    
    % function in_rescale(~, ~, maxy)
    %        clim([0 maxy]);
    % end

end


