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
        case 'Brennecke et al. (2013) [PMID:24056876]'
            fw = gui.myWaitbar(FigureHandle);
            try
                sc_hvg(sce.X, sce.g, true, true);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message);
            end
            gui.myWaitbar(FigureHandle, fw, true);
        case 'Splinefit Method [PMID:31697351]'            
            try
               gui.sc_scatter3genes(sce.X, sce.g, true, true, FigureHandle);
            catch ME
               gui.myErrordlg(FigureHandle, ME.message);
            end            
        otherwise
            return;
    end

end
