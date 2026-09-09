function [done] = callback_Harmony(src, ~)
%CALLBACK_HARMONY Batch-correct the current embedding with native Harmony.
%
%   Runs RUN.ML_HARMONY on SCE.S, grouping cells by SCE.C_BATCH_ID. No R or
%   Python installation is involved; this is the MATLAB implementation.
%
%   Harmony expects a principal-component embedding, not a two-dimensional
%   plotting embedding. Correcting UMAP or t-SNE coordinates directly is
%   possible and is what the R and Python menu entries have always done,
%   but the result is only as good as the two dimensions it is given, so
%   the user is offered the option of correcting the PCs and re-embedding.
%
%   See also RUN.ML_HARMONY, GUI.CALLBACK_HARMONYR, GUI.CALLBACK_HARMONYPY.

done = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

if numel(unique(sce.c_batch_id)) < 2
    gui.myWarndlg(FigureHandle, ...
        'No batch effect (all cells have the same SCE.C_BATCH_ID)');
    return;
end

batchid = sce.c_batch_id;
if ~isnumeric(batchid)
    batchid = findgroups(sce.c_batch_id);
end
batchid = batchid(:);

answer = gui.myQuestdlg(FigureHandle, ...
    ['Correct the principal components and re-embed, or correct the ' ...
    'current embedding in place? Correcting the PCs is what Harmony is ' ...
    'designed for; correcting the embedding is faster but works with ' ...
    'only the dimensions on screen.'], 'Harmony', ...
    {'Correct PCs and re-embed', 'Correct current embedding'}, ...
    'Correct PCs and re-embed');
if isempty(answer), return; end
usepcs = strcmp(answer, 'Correct PCs and re-embed');

fw = gui.myWaitbar(FigureHandle);
try
    if usepcs
        ndim = size(sce.s, 2);
        Xn = log1p(sc_norm(sce.X)).';
        if issparse(Xn)
            Xn = full(Xn);
        end
        [~, pcs] = pca(Xn, 'NumComponents', ...
            min(50, min(size(Xn)) - 1));

        % The uncorrected PCs are kept rather than overwritten, because
        % the no-op check below has to compare what Harmony was GIVEN with
        % what it RETURNED. It used to be `pcs = run.ml_Harmony(pcs, ...)`,
        % discarding the input, and the check then compared the freshly
        % embedded result against the OLD sce.s -- two unrelated
        % embeddings, which are never equal. Measured on a single-batch
        % dataset, where there is nothing to correct: Harmony returned the
        % PCs bit-for-bit unchanged, max|pcs - corrected| = 0, and
        % isequal(sce.s, s) was still false, the two differing by 4.647.
        % So a complete no-op was accepted and the user was handed a brand
        % new UMAP presented as batch-corrected.
        pcsCorrected = run.ml_Harmony(pcs, batchid);
        unchanged = isequal(pcs, pcsCorrected);

        % The PCs are already normalised and corrected, so SC_UMAP must not
        % normalise or log-transform them again.
        s = sc_umap(pcsCorrected.', ndim, false, false);
    else
        s = run.ml_Harmony(sce.s, batchid);
        % Harmony worked on sce.s directly here, so comparing the two IS
        % the right test -- it is only the branch above that could not use
        % it.
        unchanged = isequal(sce.s, s);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(s) || unchanged
    gui.myErrordlg(FigureHandle, ...
        ['Harmony returned the embedding unchanged. That happens when ' ...
        'the batches share no cell type, because a cluster drawn from ' ...
        'one batch alone carries no information about how the batches ' ...
        'differ.'], 'Harmony');
    return;
end

sce.s = s;
gui.myGuidata(FigureHandle, sce, src);
done = true;

end
