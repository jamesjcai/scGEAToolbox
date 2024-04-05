function callback_EnrichrHVGs(src, ~, sce)

FigureHandle = [];
if nargin < 3 || isempty(sce)
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);    
end

answer = questdlg('Which HVG detecting method to use?', '', ...
    'Splinefit Method [PMID:31697351]', ...
    'Brennecke et al. (2013) [PMID:24056876]', ...
    'Splinefit Method [PMID:31697351]');

    switch answer
        case 'Brennecke et al. (2013) [PMID:24056876]'
            fw = gui.gui_waitbar;
            try
                sc_hvg(sce.X, sce.g, true, true);
            catch ME
                gui.gui_waitbar(fw, true);
                uiwait(errordlg(ME.message,''));
            end
            gui.gui_waitbar(fw, true);
        case 'Splinefit Method [PMID:31697351]'            
            try
               gui.sc_scatter3genes_new(sce.X, sce.g, true, true, FigureHandle);
            catch ME
               uiwait(errordlg(ME.message,''));
            end            
        otherwise
            return;
    end

end
