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
        % Same three rankers, named as GUI.CALLBACK_DVGENE2GROUPS names
        % them. The analytic curve leads, matching what
        % SINGLECELLEXPERIMENT.EMBEDCELLS selects genes with.
        optSpline = 'Splinefit Method [PMID:31697351]';
        optAnalytic = 'Analytic Curve (closed-form Spline-DV)';
        optBrennecke = 'Brennecke et al. (2013) [PMID:24056876]';
        answer = gui.myQuestdlg(FigureHandle, 'Which HVG detecting method to use?', '', ...
            {optAnalytic, optSpline, optBrennecke}, optAnalytic);
        if ~ismember(answer, {optAnalytic, optSpline, optBrennecke}), return; end
        % The ranking and the error below both used to leave the bar open.
        fw = gui.myWaitbar(FigureHandle);
        try
            switch answer
                case optAnalytic
                    T = sc_analyticfit(sce.X, sce.g);
                case optSpline
                    T = sc_splinefit(sce.X, sce.g);
                otherwise   % optBrennecke; the others returned above
                    T = sc_hvg(sce.X, sce.g);
            end
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
        % HEIGHT(T) rather than SCE.NUMGENES: SC_ANALYTICFIT and
        % SC_SPLINEFIT both drop genes that are zero in every cell, so
        % the table can be shorter than the SCE and indexing out to
        % NUMGENES throws. SC_HVG keeps them, so only two of the three
        % branches could reach it.
        glist = T.genes(1:min([k, height(T)]));
        [y, idx] = ismember(glist, sce.g);
        if ~all(y)
            gui.myErrordlg(FigureHandle, 'Runtime error.','');
            return;
        end
    case 'ligandreceptor'
        % No progress bar: loading this small table is instant, and the
        % error and question below used to appear over an open one.
        mfolder = fileparts(mfilename('fullpath'));
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
    otherwise
        return;   % no other TYPE is passed
end

sce = sce.selectgenesbyindex(idx);
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end
