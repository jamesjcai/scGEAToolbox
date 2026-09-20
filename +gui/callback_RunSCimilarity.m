function [needupdatesce] = callback_RunSCimilarity(src, ~)

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

[modeldir] = gui.i_setscimilaritymodelpath(src, [], FigureHandle);
if isempty(modeldir), return; end

label_ints_file = fullfile(modeldir, 'label_ints.csv');
if exist(label_ints_file, "file")
    answer = gui.myQuestdlg(FigureHandle, 'Unconstrained or constrained annotation','', ...
        {'Unconstrained','Constrained'},'Unconstrained');
    switch answer
        case 'Unconstrained'
            target_celltypes = '';
        case 'Constrained'
            T = readtable(label_ints_file, 'ReadVariableNames',true, ...
                'VariableNamingRule', 'modify');
            allcelltypes = natsort(string(T.x0));
            % scimilmodelpath
            % scimiltargetcel
            if ispref('scgeatoolbox', 'scimiltargetcel')
                preselected_celltypes = getpref('scgeatoolbox', 'scimiltargetcel');
            else
                preselected_celltypes = '';
            end
            [idx] = gui.i_selmultidialog(allcelltypes, preselected_celltypes, FigureHandle);
            if isempty(idx), return; end
            if idx == 0, return; end
            target_celltypes = allcelltypes(idx);
            setpref('scgeatoolbox', 'scimiltargetcel', target_celltypes);
        otherwise
            return;
    end
else
    gui.myWarndlg(FigureHandle, "Missing label_ints.csv. Scimilarity model path is invalid.");
    return;
end

% SCimilarity embeds and classifies each cell on its own, so a cluster
% that is plainly one type comes back speckled with a few percent of its
% neighbours'. Two ways out of that, and they are not the same offer:
% MODE_CLUSTER still labels every cell and only summarises the result,
% while MODE_POOLED labels summed profiles and never sends most cells at
% all.
% Asked here rather than after the run, so a dataset with no clustering
% costs nothing and nobody waits on a multi-minute Python call for one
% more prompt.
MODE_PERCELL = 1;
MODE_CLUSTER = 2;
MODE_POOLED = 3;
mode = MODE_PERCELL;

if ~prepare_input_only
    % SingleCellExperiment seeds c_cluster_id to all-ones, so "has been
    % clustered" is more than one distinct id, not a non-empty vector.
    % With one cluster that vote could only paint every cell the same
    % type, so the option is withheld rather than offered and ignored.
    % Pooling makes its own groups and is always available.
    hasclusters = numel(unique(string(sce.c_cluster_id))) > 1;
    choices = {'Per cell', 'Pooled pseudobulk'};
    if hasclusters
        choices = {'Per cell', 'Majority per cluster', ...
            'Pooled pseudobulk'};
    end
    % Three is the ceiling here: GUI.MYQUESTDLG appends its own Cancel on
    % the uifigure branch, and UICONFIRM takes at most four options.
    answer = gui.myQuestdlg(FigureHandle, ...
        ['SCimilarity labels each cell independently. Label every cell, ', ...
        'summarise those labels over the existing clusters, or pool ', ...
        'cells into groups of about 50 and label one summed profile per ', ...
        'pool instead of its cells?'], ...
        'SCimilarity', choices, 'Per cell');
    switch answer
        case 'Per cell'
            mode = MODE_PERCELL;
        case 'Majority per cluster'
            mode = MODE_CLUSTER;
        case 'Pooled pseudobulk'
            mode = MODE_POOLED;
        otherwise
            return;   % dismissed
    end
end

extprogname = 'py_scimilarity';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

% SCimilarity's align_dataset matches against ca.gene_order, which is
% upper-case HGNC symbols, so the gene list does have to be upper-cased
% for the run. But SCE is a handle object shared with the app, so doing it
% in place renamed every gene in the user's live dataset for the rest of
% the session -- including on the prepare-input-only path, which reports
% that it changed nothing. Restore it when this callback returns, by any
% route; the only change it is meant to leave behind is c_cell_type_tx.
originalGeneList = sce.g;
restoreGeneList = onCleanup(@() i_restoregenelist(sce, originalGeneList));
sce.g = upper(sce.g);


if prepare_input_only
    try
        fw = gui.myWaitbar(FigureHandle);
        run.py_scimilarity(sce, modeldir, wkdir, target_celltypes, true, prepare_input_only);
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
        switch mode
            case MODE_POOLED
                [c, pooled] = i_pooledrun(sce, modeldir, wkdir, ...
                    target_celltypes);
            otherwise
                [c, cellstats] = run.py_scimilarity(sce, modeldir, wkdir, ...
                    target_celltypes, true);
        end
        assert(sce.NumCells==numel(c));
        stashname = pkg.i_stashcelltypehistory(sce);
        if mode ~= MODE_POOLED
            % Confidence in each cell's own label, which SCimilarity
            % computes while predicting and this callback used to discard.
            % Stored for every mode, so a label can always be read next to
            % how well the reference atlas agreed about it.
            confidence = i_cellconfidence(cellstats, sce.NumCells);
            if ~isempty(confidence)
                sce.setCellAttribute('scimilarity_confidence', confidence);
            end
        end
        if mode == MODE_CLUSTER
            % Keep the per-cell calls the vote summarises. They are the
            % only record of how decided each cluster was, and a cluster
            % the vote carried 51/49 reads exactly like a unanimous one
            % once c_cell_type_tx has been overwritten.
            sce.setCellAttribute('scimilarity_per_cell', string(c));
            % Weighted by each cell's own confidence, so that eight cells
            % SCimilarity barely decided cannot outvote seven it was sure
            % of. CONFIDENCE is [] when the run produced no statistics,
            % and the vote is then the plain one it always was.
            [c, nchanged] = pkg.i_majorityvote(c, sce.c_cluster_id, ...
                confidence);
        elseif mode == MODE_POOLED
            % Which pool each cell sat in, how confident its pooled label
            % was, how the pool was settled, and any label SCimilarity
            % gave the cell itself. Without them a label inferred from a
            % pool's summed profile is indistinguishable from one
            % SCimilarity gave that cell directly.
            sce.setCellAttribute('scimilarity_per_cell', pooled.percell);
            sce.setCellAttribute('scimilarity_pool', pooled.poolid);
            sce.setCellAttribute('scimilarity_status', pooled.status);
            % Two confidences, two meanings, so two attributes rather than
            % one column that silently changes what it measures:
            % 'scimilarity_confidence' is the pool profile's, carried by
            % every cell in the pool, and 'scimilarity_cell_confidence' is
            % the cell's own, NaN for cells never classified individually.
            sce.setCellAttribute('scimilarity_confidence', pooled.confidence);
            sce.setCellAttribute('scimilarity_cell_confidence', ...
                pooled.percellweight);
        end
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
    msg = sprintf('%d cell type(s) assigned to %d cells by SCimilarity.', ...
        numel(unique(string(c))), sce.NumCells);
    if mode == MODE_CLUSTER
        if isempty(confidence)
            basis = 'the label most of its cells were given';
        else
            basis = 'the label its cells backed with the most confidence';
        end
        msg = sprintf(['%s\n\nEach cluster was collapsed to %s, which ', ...
            'relabelled %d of %d cells. The per-cell labels are kept as ', ...
            'the ''scimilarity_per_cell'' cell attribute.'], ...
            msg, basis, nchanged, sce.NumCells);
    elseif mode == MODE_POOLED
        msg = sprintf(['%s\n\n%d cells were pooled into %d groups and ', ...
            'labelled from %d pseudobulk profiles, equivalent to %.0f%% ', ...
            'of a per-cell run. Median pool confidence %.2f.\n\n%d ', ...
            'pool(s) were unclear and had their %d cells labelled ', ...
            'individually; %d of those stayed mixed and kept their ', ...
            'cells'' own labels.\n\nThe pool id, its confidence, how it ', ...
            'was settled, and any labels SCimilarity returned for ', ...
            'individual cells are kept as the ''scimilarity_pool'', ', ...
            '''scimilarity_confidence'', ''scimilarity_status'' and ', ...
            '''scimilarity_per_cell'' cell attributes.'], ...
            msg, sce.NumCells, pooled.npools, pooled.nprofiles, ...
            100*(pooled.nprofiles + pooled.nlabelled)/sce.NumCells, ...
            pooled.medianconfidence, pooled.nescalated, ...
            pooled.nlabelled, pooled.nmixed);
    end
    gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));
end
end

function i_restoregenelist(sce, g)
% Put back the gene list the callback upper-cased for SCimilarity.
sce.g = g;
end


function confidence = i_cellconfidence(cellstats, ncells)
% Per-cell label confidence as vsAll_weighted, or [] when the run produced
% none - a working folder left by an older release has a script_mat.py that
% writes no stats file. Returning [] rather than ones() keeps a missing
% measurement out of the attribute table entirely, instead of storing a
% column of fabricated certainty.

confidence = [];
if isempty(cellstats) || ~ismember('vsAll_weighted', ...
        cellstats.Properties.VariableNames)
    return;
end
if height(cellstats) ~= ncells
    return;
end
confidence = cellstats.vsAll_weighted;
end


function [labels, pooled] = i_pooledrun(sce, modeldir, wkdir, target_celltypes)
%I_POOLEDRUN Label ~50-cell pools from one pseudobulk profile each.
%
% Two SCimilarity runs at most. The first sends 3 pseudobulk profiles per
% pool - the pool and each of its halves - which is 3 profiles per 50 cells
% rather than 50 cells. The second, only if some pool was unclear, sends
% that pool's actual cells. Each run pays the Python start-up and the model
% load once, which is why neither works pool by pool.
%
% Returns the per-cell labels and a struct of what the caller stores and
% reports: the labels SCimilarity returned for individual cells (percell),
% each cell's pool (poolid), how that pool was settled (status), and the
% counts the summary quotes.

coords = i_poolcoordinates(sce);
[poolid, npools] = pkg.i_partitioncells(coords, TargetSize=50);

[P, poolname] = pkg.i_poolpseudobulk(sce.X, poolid, coords);

% The pseudobulk profiles go through the ordinary path as though they were
% cells: RUN.PY_SCIMILARITY only reads X and g, and script_mat.py aligns to
% the model gene space and then normalises to 1e4 and log1p - which applied
% to a column of summed counts is exactly SCimilarity's utils.get_centroid.
pseudosce = SingleCellExperiment(P, sce.g);
[poollabels, poolstats] = run.py_scimilarity(pseudosce, modeldir, wkdir, ...
    target_celltypes, true);
if numel(poollabels) ~= 3*npools
    error('gui:callback_RunSCimilarity:pseudobulkCount', ...
        'Sent %d pseudobulk profiles and got %d labels back.', ...
        3*npools, numel(poollabels));
end

stage1 = table(poolname, poollabels(1:npools), ...
    i_poolconfidence(poolstats, npools), ...
    [poollabels(npools+1:2*npools), poollabels(2*npools+1:3*npools)], ...
    'VariableNames', {'Pool', 'Label', 'Confidence', 'HalfLabel'});

annotate = @(idx) i_annotatecells(sce, modeldir, wkdir, target_celltypes, idx);
[labels, report, percell, percellweight] = ...
    pkg.i_pooledannotate(poolid, stage1, annotate);

% Broadcast the per-pool verdict back onto cells, so the stored attributes
% can be read next to any other cell attribute without a join.
[~, where] = ismember(string(poolid), report.Pool);
pooled.percell = percell;
pooled.percellweight = percellweight;
pooled.poolid = poolid;
pooled.status = report.Status(where);
pooled.confidence = report.Confidence(where);
pooled.npools = npools;
pooled.nprofiles = 3*npools;
pooled.nlabelled = sum(strlength(percell) > 0);
pooled.nescalated = sum(report.Status ~= "consensus");
pooled.nmixed = sum(report.Status == "mixed");
pooled.medianconfidence = median(report.Confidence, 'omitnan');
end


function [labels, weights] = i_annotatecells(sce, modeldir, wkdir, ...
    target_celltypes, idx)
% The escalation classifier PKG.I_POOLEDANNOTATE calls, which must return a
% weight per label as well as the label. Splitting SCimilarity's stats table
% down to the one column that is a weight belongs here rather than in the
% helper, which has no business knowing what a SCimilarity stats table is.

[labels, stats] = run.py_scimilarity(sce, modeldir, wkdir, ...
    target_celltypes, true, false, idx);
weights = i_cellconfidence(stats, numel(labels));
end


function confidence = i_poolconfidence(poolstats, npools)
% The pool pseudobulks' label confidence, as vsAll_weighted: the winning
% label's share of the 50 reference neighbours under the same 1/distance
% weighting that chose the label.
%
% A working folder left by an older release has a script_mat.py that writes
% no stats file, so RUN.PY_SCIMILARITY hands back an empty table. Treat that
% as "no confidence information" rather than as low confidence: the run then
% rests on the halves-agree test alone, which is the behaviour that existed
% before there was a confidence at all.

if isempty(poolstats) || ~ismember('vsAll_weighted', ...
        poolstats.Properties.VariableNames)
    warning('gui:callback_RunSCimilarity:noStats', ...
        ['SCimilarity returned no prediction statistics, so pools are ', ...
        'settled on the halves-agree test alone. Delete the working ', ...
        'folder to get a current script_mat.py.']);
    confidence = ones(npools, 1);
    return;
end
confidence = poolstats.vsAll_weighted(1:npools);
end


function coords = i_poolcoordinates(sce)
%I_POOLCOORDINATES PCA coordinates to cut the pools in.
%
% A pool is only worth summing if its cells share a type, which is a
% property of the space it was cut in. PCA of log-normalised HVG counts is
% that space, following GUI.CALLBACK_SUBSAMPLECELLS. An existing 2-D
% embedding is deliberately not reused: UMAP and t-SNE are optimised for
% display, and their distortions put unrelated cells side by side, which
% is exactly the error this mode then carries across 50 cells at once.
%
% HVGs first, and SINGLE, because the FULL() below is the one place this
% path can run a large dataset out of memory.

ndim = 30;
[~, Xhvg] = sc_splinefit(sce.X, sce.g, true, false);
Xhvg = Xhvg(1:min(2000, size(Xhvg, 1)), :);
Xnorm = single(full(log1p(pkg.norm_libsize(Xhvg, 1e4))));
ndim = min(ndim, max(1, min(size(Xnorm))-1));
[~, coords] = pca(Xnorm.', NumComponents=ndim);
end
