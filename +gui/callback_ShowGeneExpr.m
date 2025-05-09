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
    an1 = 'Yes, same figure';
    an2 = 'No, different tabs';
    if n > 1
        answer = gui.myQuestdlg(FigureHandle, "Plot on all genes in the same figure?", "",...
            {an1, an2, 'Cancel'}, an1);
        if isempty(answer), return; end
        if strcmp(answer, "Cancel"), return; end
    else
        answer = an2;
    end
    
    switch answer
        case an1
            fw = gui.myWaitbar(FigureHandle); 
            hx = gui.myFigure;
            tabgp = uitabgroup();
            nf = 1;
            tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
            axes('parent', tab{nf});
            
            % hx.addCustomButton('off', @in_showgenename, 'HDF_point.gif', 'Show gene names');

            if ~ispref('scgeatoolbox', 'prefcolormapname')
                setpref('scgeatoolbox', 'prefcolormapname', 'autumn');
            end
            a = getpref('scgeatoolbox', 'prefcolormapname');
            maxy = 0;            
            for k = 1:n
                nexttile
                sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
                c = Xt(sce.g == glist(k), :);
                gui.i_setautumncolor(c, a, true, any(c==0));
                colorbar;
                maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
            end

            nf = 2;
            tab{nf} = uitab(tabgp, 'Title', sprintf('Tab%d',nf));
            axes('parent', tab{nf});
            
            for k = 1:n
                nexttile
                sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);                        
            end
            gui.myWaitbar(FigureHandle, fw);
            hx.show(FigureHandle);            

        % case an1    % same figure;
        %     answer2 = gui.myQuestdlg(FigureHandle, "Type of plot:","", "stem plot", "feature plot", "stem plot");
        %     if isempty(answer2), return; end
        %     fw = gui.myWaitbar(FigureHandle); 
        %     hx = gui.myFigure;
        %     maxy = 0;
        %     if ~ispref('scgeatoolbox', 'prefcolormapname')
        %         setpref('scgeatoolbox', 'prefcolormapname', 'autumn');
        %     end
        %     a = getpref('scgeatoolbox', 'prefcolormapname');
        %     for k = 1:n
        %         nexttile
        %         switch answer2
        %             case "feature plot"
        %                 sc_scattermarker(Xt, sce.g, sce.s, glist(k), 2, 5, false);
        %                 c = Xt(sce.g == glist(k), :);
        %                 gui.i_setautumncolor(c, a, true, any(c==0));
        %                 colorbar;
        %             case "stem plot"
        %                 sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);                        
        %         end
        %         maxy = max([maxy, max(Xt(sce.g == glist(k)))]);
        %     end
        %     gui.myWaitbar(FigureHandle, fw);
        %     hx.show(FigureHandle);
        case an2
            fw = gui.myWaitbar(FigureHandle);
            y = cell(n,1);
            for k = 1:n
                y{k} = Xt(sce.g == glist(k), :);
            end
            gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
            gui.myWaitbar(FigureHandle, fw);            
    end

 % function in_showgenename(~,~)
 %     % ancestor(ax1, 'figure')
 %     for l = 1:2
 %     ax = findall(tab{l}, 'Type', 'axes');
 %         for kk = 1:length(ax)
 %             % gui.i_export2pptx({ancestor(ax, 'figure')}, {''});
 %             colorbar(ax(kk),"off");
 %             box(ax(kk),"on");
 %             grid(ax(kk),"off");
 %             subtitle(ax(kk),'');
 %         end
 %     end
 % end

end