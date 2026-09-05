function [needupdatesce, T] = callback_RunPanhumanpy(src, ~)

needupdatesce = false;
[y, prepare_input_only] = gui.i_memorychecked([], []);
if ~y, return; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

% Preparing input files does not touch sce.c_cell_type_tx; only ask about
% overwriting labels when this run will actually assign new ones.
if ~prepare_input_only && ~gui.i_confirmoverwritecelltype(FigureHandle, sce)
    return;
end

% https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
% SCimilarity trained model. Download SCimilarity models.
% Note, this is a large tarball - downloading and uncompressing can take a several minutes.

% [modeldir] = gui.i_setscimilaritymodelpath;
% if isempty(modeldir), return; end

% label_ints_file = fullfile(modeldir, 'label_ints.csv');
% if exist(label_ints_file, "file")
%    answer = gui.myQuestdlg(FigureHandle, 'Unconstrained or constrained annotation','', ...
%        {'Unconstrained','Constrained'},'Unconstrained');
%    switch answer
%        case 'Unconstrained'
%            target_celltypes = '';
%        case 'Constrained'
%            T = readtable(label_ints_file, 'ReadVariableNames',true, ...
%                'VariableNamingRule', 'modify');
%            allcelltypes = natsort(string(T.x0));
            % scimilmodelpath
            % scimiltargetcel
%            if ispref('scgeatoolbox', 'scimiltargetcel')
%                preselected_celltypes = getpref('scgeatoolbox', 'scimiltargetcel');
%            else
%                preselected_celltypes = '';
%            end
%            [idx] = gui.i_selmultidialog(allcelltypes, preselected_celltypes, FigureHandle);
%            if isempty(idx), return; end
%            if idx == 0, return; end
%            target_celltypes = allcelltypes(idx);
%            setpref('scgeatoolbox', 'scimiltargetcel', target_celltypes);
%        otherwise
%            return;
%    end
% else
%    gui.myWarndlg(FigureHandle, "Missing label_ints.csv. Scimilarity model path is invalid.");
%    return;
% end

extprogname = 'py_panhumanpy';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
sce.g = upper(sce.g);


if prepare_input_only
    try
        fw = gui.myWaitbar(FigureHandle);
        run.py_panhumanpy(sce, wkdir, true, prepare_input_only);
        gui.myWaitbar(FigureHandle, fw);
        if strcmp(gui.myQuestdlg(FigureHandle, 'Input files prepared. Open the working folder?'),'Yes')
            winopen(wkdir);
        end
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    needupdatesce = false;
else
    fw = gui.myWaitbar(FigureHandle);
    try
        [c, T] = run.py_panhumanpy(sce, wkdir, true);
        assert(sce.NumCells==numel(c));
        stashname = pkg.i_stashcelltypehistory(sce);
        sce.c_cell_type_tx = c;
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myGuidata(FigureHandle, sce, src);
    needupdatesce = true;
    gui.myWaitbar(FigureHandle, fw);

    % Say what was assigned and where the labels it replaced went; the
    % pre-flight warning could only name the attribute generically.
    msg = sprintf('%d cell type(s) assigned to %d cells by panhumanpy.', ...
        numel(unique(string(c))), sce.NumCells);
    gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));
end
end
