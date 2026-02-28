function callback_scTenifoldXct2(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

% ── 1. Method selection ───────────────────────────────────────────────────
has_dlt = license('test', 'Neural_Network_Toolbox') && ...
          (exist('dlarray', 'builtin') || exist('dlarray', 'file'));

methodlist = {
    'Neural Network       (MATLAB Path A — Deep Learning Toolbox)',
    'Spectral + PCNet     (MATLAB Path B)',
    'Spectral + Pearson   (MATLAB Path B Lite — no extra toolbox)',
    'Python               (scTenifoldXct2 package, original)'
};
MTHD_NN     = 1;
MTHD_SPEC   = 2;
MTHD_LITE   = 3;
MTHD_PYTHON = 4;

def_midx = MTHD_NN;
if ~has_dlt, def_midx = MTHD_SPEC; end

[midx, tf] = gui.myListdlg(FigureHandle, methodlist, ...
    'Select scTenifoldXct2 implementation:', methodlist(def_midx));
if ~tf || isempty(midx), return; end

% Block Path A if Deep Learning Toolbox is absent
if midx == MTHD_NN && ~has_dlt
    gui.myErrordlg(FigureHandle, ...
        ['Deep Learning Toolbox is required for the Neural Network ' ...
         '(Path A) implementation but was not found on this system. ' ...
         'Please select a different method.'], ...
        'scTenifoldXct2:noDLT');
    return;
end

use_python = (midx == MTHD_PYTHON);

% ── 2. Python-only setup ──────────────────────────────────────────────────
prepare_input_only = false;
wkdir = [];

if use_python
    numglist = [1 3000 5000];
    memmlist = [16 32 64 128];
    neededmem = memmlist(sum(sce.NumGenes > numglist));
    [yesgohead, prepare_input_only] = gui.i_memorychecked(neededmem);
    if ~yesgohead, return; end

    extprogname = 'py_scTenifoldXct';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end

    if ~prepare_input_only
        if ~gui.i_setpyenv([], [], FigureHandle), return; end
    end
end

% ── 3. Sample (batch) selection ───────────────────────────────────────────
[~, cL_batch] = findgroups(string(sce.c_batch_id));
[j1, j2, ~, ~] = aaa(cL_batch, string(sce.c_batch_id), FigureHandle);
if isempty(j1) || isempty(j2)
    gui.myWarndlg(FigureHandle, ['All cells have the same BATCH_ID. ' ...
        'Two samples are required.']);
    return;
end
sce1 = copy(sce); sce1.selectcells(j1);
sce2 = copy(sce); sce2.selectcells(j2);

if sce1.NumCells < 50 || sce2.NumCells < 50
    if ~strcmp(gui.myQuestdlg(FigureHandle, ...
            'One sample contains too few cells (n < 50). Continue?', ...
            '', [], [], 'error'), 'Yes')
        return;
    end
end

% ── 4. Cell-type selection ────────────────────────────────────────────────
[~, cL_ct] = findgroups(string(sce.c_cell_type_tx));
[~, ~, celltype1, celltype2] = aaa(cL_ct, string(sce.c_cell_type_tx), FigureHandle);
if isempty(celltype1) || isempty(celltype2)
    gui.myWarndlg(FigureHandle, ['All cells are the same type. ' ...
        'Two different cell types are required.']);
    return;
end
celltype1 = string(celltype1);
celltype2 = string(celltype2);

% ── 5. Direction selection ────────────────────────────────────────────────
a1 = sprintf('%s -> %s', celltype1, celltype2);
a2 = sprintf('%s -> %s', celltype2, celltype1);

twosided = false;
answer = gui.myQuestdlg(FigureHandle, ...
    'Select direction: Source (ligand) -> Target (receptor)', ...
    '', {'Both', a1, a2}, 'Both');
if isempty(answer), return; end
switch answer
    case 'Both'
        ct1 = celltype1; ct2 = celltype2; twosided = true;
    case a1
        ct1 = celltype1; ct2 = celltype2;
    case a2
        ct1 = celltype2; ct2 = celltype1;
    otherwise
        return;
end

% ── 6. Run selected implementation ───────────────────────────────────────
Tres = [];
switch midx

    case MTHD_NN   % ── Path A: Neural Network ─────────────────────────────
        X_s1 = sce1.X(:, sce1.c_cell_type_tx == ct1);
        X_t1 = sce1.X(:, sce1.c_cell_type_tx == ct2);
        X_s2 = sce2.X(:, sce2.c_cell_type_tx == ct1);
        X_t2 = sce2.X(:, sce2.c_cell_type_tx == ct2);
        g    = sce1.g;
        try
            Tres = ten.xct.xctmain2_nn(X_s1, X_t1, X_s2, X_t2, g, ...
                'twosided', twosided);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_SPEC   % ── Path B: Spectral + PCNet ───────────────────────────
        try
            Tres = ten.sctenifoldxct2(sce1, sce2, ct1, ct2, twosided);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_LITE   % ── Path B Lite: Spectral + Pearson ────────────────────
        X_s1 = sce1.X(:, sce1.c_cell_type_tx == ct1);
        X_t1 = sce1.X(:, sce1.c_cell_type_tx == ct2);
        X_s2 = sce2.X(:, sce2.c_cell_type_tx == ct1);
        X_t2 = sce2.X(:, sce2.c_cell_type_tx == ct2);
        g    = sce1.g;
        try
            Tres = ten.xct.xctmain2(X_s1, X_t1, X_s2, X_t2, g, ...
                'twosided', twosided);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_PYTHON   % ── Python (original) ──────────────────────────────
        % Propagate cell-type labels expected by py_scTenifoldXct2
        sce1.c_cell_type_tx = sce1.c_cell_type_tx;
        sce2.c_cell_type_tx = sce2.c_cell_type_tx;
        try
            [Tres, ~] = run.py_scTenifoldXct2(sce1, sce2, ct1, ct2, twosided, ...
                wkdir, true, prepare_input_only, FigureHandle);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

end

% ── 7. Post-process: add direction column ─────────────────────────────────
T = [];
if twosided && iscell(Tres)
    T1 = Tres{1};
    T2 = Tres{2};
    if istable(T1) && height(T1) > 0
        a  = sprintf("%s -> %s", ct1, ct2);
        T1 = addvars(T1, repelem(a, height(T1), 1), 'Before', 1);
        T1.Properties.VariableNames{'Var1'} = 'direction';
    end
    if istable(T2) && height(T2) > 0
        a  = sprintf("%s -> %s", ct2, ct1);
        T2 = addvars(T2, repelem(a, height(T2), 1), 'Before', 1);
        T2.Properties.VariableNames{'Var1'} = 'direction';
    end
    T = [T1; T2];
elseif istable(Tres) && height(Tres) > 0
    T = Tres;
    a = sprintf("%s -> %s", ct1, ct2);
    T = addvars(T, repelem(a, height(T), 1), 'Before', 1);
    T.Properties.VariableNames{'Var1'} = 'direction';
end

% ── 8. Annotate with extended L-R database and export ────────────────────
if ~isempty(T)

    mfolder = fileparts(mfilename('fullpath'));
    load(fullfile(mfolder, '..', 'assets', 'Ligand_Receptor', ...
         'Ligand_Receptor_more.mat'), 'ligand', 'receptor');
    A = [upper(string(T.ligand)), upper(string(T.receptor))];
    B = [upper(string(ligand)),   upper(string(receptor))];
    knownpair = ismember(A, B, 'rows');
    assert(length(knownpair) == height(T));
    T = [T, table(knownpair)];

    if use_python && ~isempty(wkdir) && isfolder(wkdir)
        % Python path: save to working directory
        b = matlab.lang.makeValidName(string(datetime));
        b = fullfile(wkdir, b + ".txt");
        writetable(T, b);
        answer2 = gui.myQuestdlg(FigureHandle, ...
            sprintf('Result saved in %s', b), ...
            '', {'Export result...', 'Locate result file...'}, ...
            'Export result...');
        switch answer2
            case 'Locate result file...'
                winopen(wkdir); pause(2)
                if strcmp(gui.myQuestdlg(FigureHandle, 'Export result to other format?'), 'Yes')
                    gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table', ...
                        [], [], FigureHandle);
                end
            case 'Export result...'
                gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table', ...
                    [], [], FigureHandle);
            otherwise
                winopen(wkdir);
        end
        return;
    end

    % MATLAB paths: save to temp + offer export
    [a, b] = pkg.i_tempdirfile("sctenifoldxct2");
    writetable(T, b);
    answer2 = gui.myQuestdlg(FigureHandle, sprintf('Result saved in %s', b), '', ...
        {'Export result...', 'Locate result file...'}, 'Export result...');
    switch answer2
        case 'Locate result file...'
            winopen(a); pause(2)
            if strcmp(gui.myQuestdlg(FigureHandle, 'Export result to other format?'), 'Yes')
                gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table', ...
                    [], [], FigureHandle);
            end
        case 'Export result...'
            gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table', ...
                [], [], FigureHandle);
        otherwise
            winopen(a);
    end

else
    if use_python && prepare_input_only
        if strcmp(gui.myQuestdlg(FigureHandle, ...
                'Input files prepared. Open working folder?', ''), 'Yes')
            winopen(wkdir);
        end
    else
        gui.myHelpdlg(FigureHandle, 'No differential ligand-receptor pairs identified.');
    end
end

end


% ── Helper: select 2 groups from a list ───────────────────────────────────
function [i1, i2, cL1, cL2] = aaa(listitems, ci, FigureHandle)
    i1 = []; i2 = []; cL1 = []; cL2 = [];
    n = length(listitems);
    if n < 2, return; end

    if gui.i_isuifig(FigureHandle)
        [indx, tf] = gui.myListdlg(FigureHandle, listitems, ...
            'Select two groups:', listitems([n-1, n]));
    else
        [indx, tf] = listdlg('PromptString', {'Select two groups:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'InitialValue', [n-1, n], 'ListSize', [220, 300]);
    end

    if tf == 1
        if numel(indx) ~= 2
            gui.myErrordlg(FigureHandle, 'Please select exactly 2 groups.');
            return;
        end
        cL1 = listitems(indx(1));
        cL2 = listitems(indx(2));
        i1  = ci == cL1;
        i2  = ci == cL2;
    end
end
