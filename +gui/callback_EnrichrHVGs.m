function callback_EnrichrHVGs(src, ~, sce)

FigureHandle = [];
if nargin < 3 || isempty(sce)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

answer = gui.myQuestdlg(FigureHandle, 'Which HVG detecting method to use?', '', ...
{'Splinefit Method [PMID:31697351]', ...
'Brennecke et al. (2013) [PMID:24056876]'}, ...
'Splinefit Method [PMID:31697351]');

switch answer
        case 'Splinefit Method [PMID:31697351]'
            fw = gui.myWaitbar(FigureHandle);
            try
                [T] = sc_splinefit(sce.X, sce.g, true, false);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message);
                return;
            end
            gui.myWaitbar(FigureHandle, fw);
            figtab = gui.TableViewerApp(T, FigureHandle);
            if strcmp(gui.myQuestdlg(figtab, ...
                    'Explore HVG expression profile of genes?'), 'Yes')
                try
                    gui.i_hvgsplinefitplot(sce.X, sce.g, true, true, figtab);
                catch ME
                    gui.myErrordlg(figtab, ME.message);
                end
            end

        case 'Brennecke et al. (2013) [PMID:24056876]'
            fw = gui.myWaitbar(FigureHandle);
            try
                [T] = sc_hvg(sce.X, sce.g, true, false);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message);
                return;
            end
            gui.myWaitbar(FigureHandle, fw);
            figtab = gui.TableViewerApp(T, FigureHandle);
            if strcmp(gui.myQuestdlg(figtab, ...
                    'Explore HVG expression profile of genes?'), 'Yes')
                try
                    sc_hvg(sce.X, sce.g, true, true);
                catch ME
                    gui.myErrordlg(figtab, ME.message);
                end
            end

        otherwise
            return;
    end

end
