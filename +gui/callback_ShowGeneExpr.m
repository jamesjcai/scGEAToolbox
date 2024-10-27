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

    % answer = requestdlg("Plot on the same figure?","", "Yes, same figure", ...
    %     "No, different figures", "Cancel", "Yes, same figure");
    % if strcmp(answer, "Cancel"), return; end

    fw=gui.gui_waitbar;
        n = length(glist);
        y = cell(n,1);
        for k=1:n
            y{k} = Xt(sce.g == glist(k), :);
        end
        gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);

        % switch answer
        %     case "No, different figures"
        %         gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
        %     case "Yes, same figure"
        %         figure;
        %         nexttile
        % 
        % 
        % end
    gui.gui_waitbar(fw);
end