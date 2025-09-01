function callback_ShowGeneExpr(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist) 
        % gui.myHelpdlg(FigureHandle, 'No genes is selected or found.');
        return;
    end

    % answer = gui.myQuestdlg(FigureHandle, "Select the type of expression values","",...
    %     "Raw UMI Counts","Library Size-Normalized",)

    [Xt] = gui.i_transformx(sce.X,[],[],FigureHandle);
    if isempty(Xt), return; end
    n = length(glist);
    if ~ispref('scgeatoolbox', 'prefcolormapname')
        setpref('scgeatoolbox', 'prefcolormapname', 'autumn');
    end
    %a = getpref('scgeatoolbox', 'prefcolormapname');

    % if ~isMATLABReleaseOlderThan('R2025a')
    %     for k = 1:n
    %         figure;
    %         sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
    %         c = Xt(sce.g == glist(k), :);
    %         ax = gca;
    %         gui.i_setautumncolor(c, a, true, any(c==0), ax);
    %         colorbar(ax);
    %     end
    % else

    
        % an1 = 'Yes, same figure';
        % an2 = 'No, different tabs';
%{
        if n > 1
            answer = gui.myQuestdlg(FigureHandle, "Plot on all genes in the same figure?", "",...
                {an1, an2, 'Cancel'}, an1);
            if isempty(answer), return; end
            if strcmp(answer, "Cancel"), return; end
        else
            answer = an2;
        end
    %}  

        % answer = an2;
        % switch answer
        %     case an1    % 'Yes, same figure';
        %         fw = gui.myWaitbar(FigureHandle); 
        % 
        %         hx = gui.myFigure(FigureHandle);
        %         tabgp = uitabgroup(hx.FigHandle);
        %         nf = 1;
        %         tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
        %         axes('parent', tab{nf});
        % 
        %         % hx.addCustomButton('off', @in_showgenename, 'HDF_point.gif', 'Show gene names');
        % 
        %         maxy = 0;
        %         for k = 1:n
        %             nexttile
        %             sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
        %             c = Xt(sce.g == glist(k), :);
        %             gui.i_setautumncolor(c, a, true, any(c==0));
        %             colorbar;
        %             maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
        %         end    
        %         nf = 2;
        %         tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
        %         axes('parent', tab{nf});                
        %         for k = 1:n
        %             nexttile
        %             sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);                        
        %         end
        %         gui.myWaitbar(FigureHandle, fw);
        %         hx.show(FigureHandle);
        %     case an2    % 'No, different tabs';
                fw = gui.myWaitbar(FigureHandle);
                y = cell(n, 1);
                for k = 1:n
                    y{k} = Xt(sce.g == glist(k), :);
                end
                gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
                gui.myWaitbar(FigureHandle, fw);
        % end
    % end
end