function callback_scTenifoldXct(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

if ~gui.gui_showrefinfo('scTenifoldXct [PMID:36787742]', FigureHandle), return; end

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
'Python               (scTenifoldXct package, original)'
};
MTHD_NN     = 1;
MTHD_SPEC   = 2;
MTHD_REF    = 3;
MTHD_PYTHON = 4;

def_midx = MTHD_NN;
if ~has_dlt, def_midx = MTHD_SPEC; end
if ~has_stats, def_midx = MTHD_PYTHON; end

[midx, tf] = gui.myListdlg(FigureHandle, methodlist, ...
'Select scTenifoldXct implementation:', methodlist(def_midx));
if ~tf || isempty(midx), return; end

% Block Path A if Deep Learning Toolbox is absent
if midx == MTHD_NN && ~has_dlt
    gui.myErrordlg(FigureHandle, ...
        ['Deep Learning Toolbox is required for the Neural Network ' ...
         '(Path A) implementation but was not found on this system. ' ...
         'Please select a different method.'], ...
        'scTenifoldXct:noDLT');
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
        'scTenifoldXct:noStats');
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
        ['scTenifoldXct joins the two cell types'' gene networks through a ' ...
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
         'Choose Outer product for parity with the reference, or keep ' ...
         'Ligand-receptor for a database-driven analysis.'], ...
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

% ── 2. Python-only setup (wkdir + pyenv) ─────────────────────────────────
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

% ── 3. Cell-type grouping selection ──────────────────────────────────────
[thisc, clabel] = gui.i_select1class(sce, false, ...
'Select grouping variable (cell type):', 'Cell Type', FigureHandle);
if isempty(thisc), return; end

if ~strcmp(clabel, 'Cell Type')
    if ~strcmp(gui.myQuestdlg(FigureHandle, ...
            ['You selected grouping variable other than ''Cell Type''.' ...
            ' Continue?'], '', [], [], 'warning'), 'Yes')
        return;
    end
end

[c, cL] = findgroups(string(thisc));
[idx] = gui.i_selmultidialog(cL, [], FigureHandle);
if isempty(idx), return; end
if numel(idx) < 2
    gui.myWarndlg(FigureHandle, ['Need at least 2 cell groups to ' ...
        'perform cell-cell interaction analysis.']);
    return;
end
if numel(idx) ~= 2
    gui.myWarndlg(FigureHandle, ...
        sprintf(['Need exactly 2 cell groups for cell-cell interaction ' ...
        'analysis. You selected %d.'], numel(idx)));
    return;
end

i1 = idx(1);
i2 = idx(2);

a1 = sprintf('%s -> %s', cL{i1}, cL{i2});
a2 = sprintf('%s -> %s', cL{i2}, cL{i1});

% ── 4. Direction selection ────────────────────────────────────────────────
twosided = false;
answer = gui.myQuestdlg(FigureHandle, ['Select direction: ' ...
'Source (ligand) -> Target (receptor)'], '', ...
{'Both', a1, a2}, 'Both');
if isempty(answer), return; end
switch answer
    case 'Both'
        x1 = i1; x2 = i2; twosided = true;
    case a1
        x1 = i1; x2 = i2;
    case a2
        x1 = i2; x2 = i1;
    otherwise
        return;
end

% ── 5. Prepare filtered SCE ──────────────────────────────────────────────
sce.c_batch_id = thisc;
sce.c_batch_id(c == x1) = "Source";
sce.c_batch_id(c == x2) = "Target";
sce.c_cell_type_tx = string(cL(c));

idx = c == x1 | c == x2;
sce = sce.selectcells(idx);

% ── 6. Run selected implementation ───────────────────────────────────────
Tres = [];
switch midx

    case MTHD_NN   % ── Path A: Neural Network ─────────────────────────────
        X_s = sce.X(:, sce.c_batch_id == "Source");
        X_t = sce.X(:, sce.c_batch_id == "Target");
        g   = sce.g;
        try
            Tres = ten.xct.xctmain_nn(X_s, X_t, g, 'twosided', twosided, ...
                'w12mode', w12mode, 'useparallel', useparallel);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_SPEC   % ── Path B: spectral, toolbox GRN settings ─────────────
        try
            Tres = ten.sctenifoldxct(sce, string(cL{x1}), string(cL{x2}), ...
                twosided, 'w12mode', w12mode, 'useparallel', useparallel);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_REF   % ── Path B: spectral, reference GRN settings ────────────
        X_s = sce.X(:, sce.c_batch_id == "Source");
        X_t = sce.X(:, sce.c_batch_id == "Target");
        g   = sce.g;
        try
            Tres = ten.xct.xctmain(X_s, X_t, g, 'twosided', twosided, ...
                'w12mode', w12mode, 'useparallel', useparallel);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_PYTHON   % ── Python (original) ──────────────────────────────
        try
            Tres = run.py_scTenifoldXct(sce, cL{x1}, cL{x2}, twosided, ...
                wkdir, true, prepare_input_only, FigureHandle);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

end

% ── 7. Post-process: unify two-sided results and add direction column ─────
T = [];
if twosided && iscell(Tres)
    T1 = Tres{1};
    T2 = Tres{2};
    if istable(T1)
        a = sprintf("%s -> %s", cL{x1}, cL{x2});
        T1 = addvars(T1, repelem(a, height(T1), 1), 'Before', 1);
        T1.Properties.VariableNames{'Var1'} = 'direction';
    end
    if istable(T2)
        a = sprintf("%s -> %s", cL{x2}, cL{x1});
        T2 = addvars(T2, repelem(a, height(T2), 1), 'Before', 1);
        T2.Properties.VariableNames{'Var1'} = 'direction';
    end
    T = [T1; T2];
elseif istable(Tres)
    T = Tres;
    a = sprintf("%s -> %s", cL{x1}, cL{x2});
    T = addvars(T, repelem(a, height(T), 1), 'Before', 1);
    T.Properties.VariableNames{'Var1'} = 'direction';
end

% ── 8. Annotate with extended L-R database and export ────────────────────
if ~isempty(T)
    mfolder = fileparts(mfilename('fullpath'));
    load(fullfile(mfolder, '..', 'assets', 'Ligand_Receptor', ...
         'Ligand_Receptor_more.mat'), 'ligand', 'receptor');
    A = [string(T.ligand) string(T.receptor)];
    B = [ligand receptor];
    [knownpair] = ismember(A, B, 'rows');
    assert(length(knownpair) == height(T));
    T = [T, table(knownpair)];

    if use_python && ~isempty(wkdir) && isfolder(wkdir)
        % Python path: save to working directory (preserves original behaviour)
        outfile = fullfile(wkdir, "outfile.csv");
        if isfile(outfile)
            answerx = gui.myQuestdlg(FigureHandle, ...
                sprintf('Overwrite %s? Select No to save in a temporary file.', outfile), ...
                '', [], [], 'warning');
        else
            answerx = 'Yes';
        end
        if ~isempty(wkdir) && isfolder(wkdir) && strcmp(answerx, 'Yes')
            writetable(T, outfile);
            if strcmp(gui.myQuestdlg(FigureHandle, ...
                    sprintf('Result saved in %s. Open working folder?', outfile)), 'Yes')
                winopen(wkdir);
            end
            return;
        end
    end

    % MATLAB paths (and Python fallback): save to temp + offer export
    [a, b] = pkg.i_tempdirfile("sctendifoldxct");
    writetable(T, b);
    answer = gui.myQuestdlg(FigureHandle, sprintf('Result saved in %s', b), '', ...
        {'Export result...', 'Locate result file...'}, 'Export result...');
    switch answer
        case 'Locate result file...'
            winopen(a);
            pause(2)
            if strcmp(gui.myQuestdlg(FigureHandle, 'Export result to other format?'), 'Yes')
                gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable', [], [], FigureHandle);
            end
        case 'Export result...'
            gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable');
        otherwise
            winopen(a);
    end

else
    if use_python && prepare_input_only
        if strcmp(gui.myQuestdlg(FigureHandle, ...
                'Input files prepared. Open working folder?'), 'Yes')
            winopen(wkdir);
        end
    else
        gui.myHelpdlg(FigureHandle, 'No ligand-receptor pairs identified.');
    end
end

end
