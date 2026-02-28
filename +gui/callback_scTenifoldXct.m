function callback_scTenifoldXct(src, ~)

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

if ~gui.gui_showrefinfo('scTenifoldXct [PMID:36787742]', FigureHandle), return; end

% ── 1. Method selection ───────────────────────────────────────────────────
% Default to Path A (Neural Network) when Deep Learning Toolbox is available,
% otherwise fall back to Path B (Spectral + PCNet).
has_dlt = license('test', 'Neural_Network_Toolbox') && ...
          (exist('dlarray', 'builtin') || exist('dlarray', 'file'));

methodlist = {
    'Neural Network       (MATLAB Path A — Deep Learning Toolbox)',
    'Spectral + PCNet     (MATLAB Path B)',
    'Spectral + Pearson   (MATLAB Path B Lite — no extra toolbox)',
    'Python               (scTenifoldXct package, original)'
};
MTHD_NN     = 1;
MTHD_SPEC   = 2;
MTHD_LITE   = 3;
MTHD_PYTHON = 4;

def_midx = MTHD_NN;
if ~has_dlt, def_midx = MTHD_SPEC; end

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

use_python = (midx == MTHD_PYTHON);

% ── 2. Python-only setup (wkdir + pyenv) ─────────────────────────────────
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
            Tres = ten.xct.xctmain_nn(X_s, X_t, g, 'twosided', twosided);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_SPEC   % ── Path B: Spectral + PCNet ───────────────────────────
        try
            Tres = ten.sctenifoldxct(sce, string(cL{x1}), string(cL{x2}), twosided);
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end

    case MTHD_LITE   % ── Path B Lite: Spectral + Pearson ────────────────────
        X_s = sce.X(:, sce.c_batch_id == "Source");
        X_t = sce.X(:, sce.c_batch_id == "Target");
        g   = sce.g;
        try
            Tres = ten.xct.xctmain(X_s, X_t, g, 'twosided', twosided);
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
                gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable');
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
