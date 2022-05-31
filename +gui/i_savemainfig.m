function i_savemainfig(src,~,tag)
    FigureHandle=src.Parent.Parent;
    if tag==1
        filter = {'*.svg'};
        [filename,filepath] = uiputfile(filter);
        if ischar(filename)
            saveas(FigureHandle,[filepath filename],'svg');
        end
    elseif tag==2
        % axx=gca;
        filter = {'*.jpg';'*.png';'*.tif';'*.pdf';'*.eps'};
        [filename,filepath] = uiputfile(filter);
        if ischar(filename)
            exportgraphics(FigureHandle,[filepath filename]);
        end
    elseif tag==3
        gui.i_export2pptx({FigureHandle},{'STGEATOOL'});
    end
end
