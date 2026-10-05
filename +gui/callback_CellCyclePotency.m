function [needupdate] = callback_CellCyclePotency(src, ~, typeid)

if nargin < 3, typeid = 1; end

needupdate = false;
needestimt = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

switch typeid
    case 1
        answer = gui.myQuestdlg(FigureHandle, 'This function assigns cell cycle phase to each cell. Continue?');
    case 2
        answer = gui.myQuestdlg(FigureHandle, 'This function assigns differentiation potency [PMID:33244588] to each cell. Continue?');
    case 3
        answer = gui.myQuestdlg(FigureHandle, 'This function calculates stemness index [PMID:29625051] for each cell. Continue?');
    case 4
        answer = gui.myQuestdlg(FigureHandle, 'This function calculates the expression ratio of dissociation-associated genes [PMID:34020534] for each cell. Continue?');
    case 5
        answer = gui.myQuestdlg(FigureHandle, 'This function predicts tumor (aneuploid) and normal (diploid) cells using copykat [PMID:33462507]. Continue?');
        extprogname = 'R_copykat';
    case 6
        answer = gui.myQuestdlg(FigureHandle, 'SCEVAN [PMID:36841879] is a fast variational algorithm for the deconvolution of the clonal substructure of tumors from single-cell RNA-seq data. Continue?');
        extprogname = 'R_SCEVAN';
end
if ~strcmp(answer, 'Yes'), return; end

switch typeid
    case {2, 6}
        speciestag = gui.i_selectspecies(2, false, FigureHandle);
        if isempty(speciestag), return; end
end

switch typeid
    case {5, 6}
        preftagname = 'externalwrkpath';
        [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
        if isempty(wkdir), return; end
end


switch typeid
    case 1
        if isempty(sce.c_cell_cycle_tx) || all(strcmp(unique(sce.c_cell_cycle_tx), "undetermined"))
            needestimt = true;
        else
            answer1 = gui.myQuestdlg(FigureHandle, 'Use existing cell cycle estimation or re-compute new estimation?', ...
                '', {'Use existing', 'Re-compute', 'Cancel'}, 'Use existing');
            switch answer1
                case 'Re-compute'
                    needestimt = true;
                case 'Use existing'
                    % nothing to compute
                otherwise
                    return;   % Cancel, or the dialog closed
            end
        end
        if needestimt
            fw = gui.myWaitbar(FigureHandle);
            try
                sce = sce.estimatecellcycle(true, 1);
            catch ME
                gui.myWaitbar(FigureHandle, fw, true);
                gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                return;
            end
            needupdate = true;
            gui.myWaitbar(FigureHandle, fw);
            gui.myGuidata(FigureHandle, sce, src);
        end
        % Show the phases in the main plot, as Assign Cell Type does.
        if isprop(src, 'sce')
            [src.c, src.cL] = findgroups(string(sce.c_cell_cycle_tx));
            src.sce.c = src.c;
            src.in_RefreshAll(true, false);
            src.ix_labelclusters(true);
        end
        if needestimt
            gui.myHelpdlg(FigureHandle, ['Cell cycle phase (c_cell_cycle_tx) added. ' ...
                'Cells are now colored by phase.']);
        end
        return;
    case 2
        attribtag = "cell_potency";
        in_aaa(attribtag);
    case 3
        attribtag = "stemness_index";
        in_aaa(attribtag);
    case 4
        attribtag = 'dissociation_ratio';
        in_aaa(attribtag);
    case 5
        attribtag = 'copykat_prediction';
        in_aaa(attribtag);
    case 6
        attribtag = 'scevan_prediction';
        in_aaa(attribtag);
end


function [s] = in_aaa(attribtag)
s = [];
if ~sce.hasCellAttribute(attribtag)
    needestimt = true;
else
    answer1 = gui.myQuestdlg(FigureHandle, sprintf('Use existing %s estimation or re-compute new estimation?', ...
        attribtag), '', {'Use existing', 'Re-compute', 'Cancel'}, 'Use existing');
    switch answer1
        case 'Re-compute'
            needestimt = true;
        case 'Use existing'
            % nothing to compute
        otherwise
            return;   % Cancel, or the dialog closed
    end
end
if needestimt
    % Asked before the progress bar opens: asking inside it put the dialog
    % over an open bar, and a cancel returned with the bar still up.
    if strcmp(attribtag, 'copykat_prediction')
        speciesid = gui.i_selectspecies(2, false, FigureHandle);
        if strlength(speciesid) == 0, return; end
    end
    fw = gui.myWaitbar(FigureHandle);
    try
        switch attribtag
            case 'cell_potency'
                s = sc_potency(sce.X, sce.g, speciestag);
            case 'stemness_index'
                s = sc_stemness(sce.X, sce.g);
            case 'dissociation_ratio'
                s = pkg.sc_dissratio(sce.X, sce.g, true);
            case 'copykat_prediction'
                s = run.r_copykat(sce, wkdir, speciesid);
            case 'scevan_prediction'
                s = run.r_SCEVAN(sce, wkdir, false, speciestag);
            otherwise
                error('Invalid attribtag');
        end
    catch ME
        % The estimators and the R runners can throw; the bar used to stay up.
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    if isempty(s)
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, sprintf("%s runtime error.", attribtag),"")
        return;
    end
    sce.setCellAttribute(attribtag, s);
    gui.myWaitbar(FigureHandle, fw);
    gui.myGuidata(FigureHandle, sce, src);
    needupdate = true;
    gui.myHelpdlg(FigureHandle, sprintf(['%s added. To see the result, ' ...
        'use View -> Cell State (Ctrl + T). Then select "%s"'], ...
        attribtag, attribtag));
else
    gui.myHelpdlg(FigureHandle, sprintf(['To see the result, use View -> ' ...
        'Cell State (Ctrl + T). Then select "%s"'], ...
        attribtag));
end


end

end
