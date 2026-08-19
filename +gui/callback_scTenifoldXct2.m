function callback_scTenifoldXct2(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

% ── 1. Method selection ───────────────────────────────────────────────────
% All three MATLAB paths build PCNet networks through net.pcrnet, so all three
% need Statistics and Machine Learning Toolbox. Only Path A needs Deep Learning
% on top of that.
has_dlt = license('test', 'Neural_Network_Toolbox') && ...
          (exist('dlarray', 'builtin') || exist('dlarray', 'file'));
has_stats = license('test', 'Statistics_Toolbox') && ...
          (exist('zscore', 'file') > 0);

methodlist = {
'Neural Network       (Path A — reproduces the published training loop)';
'Spectral             (Path B — toolbox settings: PCNet, 0.75 edge filter)';
'Spectral, reference  (Path B — PCNet nComp=5, unfiltered, as core.py)';
'Python               (scTenifoldXct2 package, original)'
};
MTHD_NN     = 1;
MTHD_SPEC   = 2;
MTHD_REF    = 3;
MTHD_PYTHON = 4;

def_midx = MTHD_NN;
if ~has_dlt, def_midx = MTHD_SPEC; end
if ~has_stats, def_midx = MTHD_PYTHON; end

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

% Block every MATLAB path if Statistics and Machine Learning Toolbox is absent.
% This used to be advertised as a toolbox-free option, which stopped being true
% when the within-type networks moved from Pearson co-expression to PCNet.
if ismember(midx, [MTHD_NN, MTHD_SPEC, MTHD_REF]) && ~has_stats
    gui.myErrordlg(FigureHandle, ...
        ['Statistics and Machine Learning Toolbox is required for the MATLAB ' ...
         'implementations, which build their within-type networks with PCNet ' ...
         '(net.pcrnet), but it was not found on this system. Select the ' ...
         'Python implementation instead.'], ...
        'scTenifoldXct2:noStats');
    return;
end

use_python = (midx == MTHD_PYTHON);

% ── 1b. Correspondence matrix (W12) mode ─────────────────────────────────
% The two modes do not answer the same question, so this is a real choice
% rather than a tuning knob, and it is worth stating plainly before asking.
% Only the MATLAB paths expose it: the Python package always rebuilds its own
% correspondence from expression (query_DB defaults to None) and never consults
% the L-R database during alignment, so there is nothing to select there.
w12mode = "lr";
if ~use_python
    OPT_LR = 'Ligand-receptor (default)';
    OPT_OUTER = 'Outer product';
    answerw = gui.myQuestdlg(FigureHandle, ...
        ['scTenifoldXct2 joins each sample''s two gene networks through a ' ...
         'correspondence matrix W12. How that matrix is built decides which ' ...
         'question the analysis answers, so the two modes are not ' ...
         'interchangeable.' newline newline ...
         'Ligand-receptor: one entry per known L-R pair from the built-in ' ...
         'database, so the alignment is driven by curated biology. This is ' ...
         'the toolbox default.' newline newline ...
         'Outer product: a dense expression-derived score over every gene ' ...
         'pair, reproducing the published Python implementation. The L-R ' ...
         'database takes no part in the alignment and is used only to rank ' ...
         'the results afterwards.' newline newline ...
         'The same mode is applied to both samples, which is what makes their ' ...
         'distances comparable in the differential test.'], ...
        'Correspondence matrix (W12)', {OPT_LR, OPT_OUTER}, OPT_LR);
    if isempty(answerw), return; end
    if strcmp(answerw, OPT_OUTER), w12mode = "outer"; end
end

% ── 1c. Parallel computing ───────────────────────────────────────────────
% net.pcrnet regresses each gene on the principal components of the others.
% That loop dominates the runtime and runs as a parfor when asked. Follows the
% idiom in callback_scTenifoldNet1lite: a GPU wins if present, because
% net.pcrnet cannot combine parfor with gpuArray.
useparallel = false;
if ~use_python
    if pkg.i_usegpu(sce.X)
        disp('GPU detected — using CUDA GPU acceleration.');
    else
        answerp = gui.myQuestdlg(FigureHandle, ...
            ['Build the gene regulatory networks with parallel computing?' ...
             newline newline ...
             'Usually NOT worth it. net.pcrnet runs one SVD per gene, and ' ...
             'the serial loop is already parallel: multithreaded BLAS spreads ' ...
             'each SVD across every core, worth 3.0x at 600 genes and 5.7x ' ...
             'at 1856. A parfor fragments that work and fights it.' ...
             newline newline ...
             'Measured on 20 cores, 784 cells, against serial with no pool: ' ...
             'roughly break-even at 600 genes, 0.5x at 1200 and 0.3x at ' ...
             '1856 - about three times slower. A Threads pool is no better. ' ...
             'Worth trying only for small gene sets.'], ...
            'Parallel Computing', ...
            {'Not use parallel', 'Use parallel'}, 'Not use parallel');
        if isempty(answerp), return; end
        useparallel = strcmp(answerp, 'Use parallel');
    end
end

% ── 2. Python-only setup ──────────────────────────────────────────────────
prepare_input_only = false;
wkdir = [];

if use_python
    numglist = [1 3000 5000];
    memmlist = [16 32 64 128];
    neededmem = memmlist(sum(sce.NumGenes > numglist));
    [yesgohead, prepare_input_only] = gui.i_memorychecked(neededmem, FigureHandle);
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
                'twosided', twosided, 'w12mode', w12mode, 'useparallel', useparallel);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_SPEC   % ── Path B: spectral, toolbox GRN settings ─────────────
        try
            Tres = ten.sctenifoldxct2(sce1, sce2, ct1, ct2, twosided, ...
                'w12mode', w12mode, 'useparallel', useparallel);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_REF   % ── Path B: spectral, reference GRN settings ────────────
        X_s1 = sce1.X(:, sce1.c_cell_type_tx == ct1);
        X_t1 = sce1.X(:, sce1.c_cell_type_tx == ct2);
        X_s2 = sce2.X(:, sce2.c_cell_type_tx == ct1);
        X_t2 = sce2.X(:, sce2.c_cell_type_tx == ct2);
        g    = sce1.g;
        try
            Tres = ten.xct.xctmain2(X_s1, X_t1, X_s2, X_t2, g, ...
                'twosided', twosided, 'w12mode', w12mode, 'useparallel', useparallel);
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
