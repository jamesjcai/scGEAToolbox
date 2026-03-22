function [requirerefresh] = callback_WorkonSelectedGenes(src, ~, type)

requirerefresh = false;
if nargin < 3, type = 'name'; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

switch type
    case 'name'
        [glist] = gui.i_selectngenes(sce, [], FigureHandle);
        if isempty(glist), return; end
        [y, idx] = ismember(glist, sce.g);
        if ~all(y)
            gui.myErrordlg(FigureHandle, 'Runtime error.','');
            return;
        end
    case 'hvg'
        k = gui.i_inputnumk(2000, 1, sce.NumGenes, ...
            'Enter the number of HVGs', FigureHandle);
        if isempty(k), return; end
        answer = gui.myQuestdlg(FigureHandle, 'Which HVG detecting method to use?', '', ...
            {'Splinefit Method [PMID:31697351]', ...
            'Brennecke et al. (2013) [PMID:24056876]'}, ...
            'Splinefit Method [PMID:31697351]');
        switch answer
            case 'Splinefit Method [PMID:31697351]'
                fw = gui.myWaitbar(FigureHandle);
                T = sc_splinefit(sce.X, sce.g);
            case 'Brennecke et al. (2013) [PMID:24056876]'
                fw = gui.myWaitbar(FigureHandle);
                T = sc_hvg(sce.X, sce.g);
            otherwise
                return;
        end
        glist = T.genes(1:min([k, sce.NumGenes]));
        [y, idx] = ismember(glist, sce.g);
        if ~all(y)
            gui.myErrordlg(FigureHandle, 'Runtime error.','');
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
    case 'ligandreceptor'
        mfolder = fileparts(mfilename('fullpath'));
        fw = gui.myWaitbar(FigureHandle);
        load(fullfile(mfolder, '..', 'assets', 'Ligand_Receptor', ...
             'Ligand_Receptor_more.mat'), 'ligand','receptor');
        idx = ismember(upper(sce.g), unique([ligand; receptor]));
        if ~any(idx)
            gui.myErrordlg(FigureHandle, 'Runtime error: No gene left after selection.','');
            return;
        end
        if sum(idx) < 50
            if ~strcmp(gui.myQuestdlg(FigureHandle, 'Few genes (n < 50) selected. Continue?',''), 'Yes'), return; end
        end
        gui.myWaitbar(FigureHandle, fw);
end

sce.g = sce.g(idx);
sce.X = sce.X(idx, :);
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end
