function [requirerefresh] = callback_GenerateSyntheticCells(src, ptsSelected)

requirerefresh = false;
if nargin < 1, ptsSelected = []; end
if isempty(ptsSelected), return; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

nBrushed = nnz(ptsSelected);
nTotal = sce.NumCells;  % gene-by-cell matrix assumption

msg = sprintf(['Generate %d synthetic cells based on %d brushed cells using generative Fourier transformer [PMID: 39843603]? ' ...
               'Dataset size will become %d cells.'], ...
               nBrushed, nBrushed, nTotal + nBrushed);

answer = gui.myQuestdlg(FigureHandle, msg);
if ~strcmp(answer,'Yes'), return; end


oldcn = sce.NumCells;
oldgn = sce.NumGenes;

% The copy "Cancel Changes" restores is taken before any edit. The
% is_synthetic attribute used to be added first, so it survived a cancel.
if sce.NumCells*sce.NumGenes < 4e8
    sceori = copy(sce);
else
    answer = gui.myQuestdlg(FigureHandle, ...
        'You are about to change the SCE data. This cannot be undone.');
    if ~strcmp(answer, 'Yes'), return; end
    sceori = [];
end

tag = 'is_synthetic';
issynthetic = sce.getCellAttribute(tag);
if isempty(issynthetic)
    issynthetic = false(sce.NumCells, 1);
    sce.setCellAttribute(tag, issynthetic);
end

% Generation and merge: no progress bar or try before, so a failure was a
% raw "Error while evaluating Menu Callback".
fw = gui.myWaitbar(FigureHandle, [], false, 'Generating synthetic cells...');
try
    sub_sce = copy(sce).selectcells(ptsSelected);
    x = sce.X(:, ptsSelected);
    x = sc_ifft(x);
    sub_sce.X = x;
    s = sub_sce.s;
    r = trnd(100, size(s)) * 0.75;
    sub_sce.s = s + r;  % * (max(range(s)));

    sub_sce.setCellAttribute(tag, true(sub_sce.NumCells, 1));

    % Each synthetic cell keeps the batch of the cell it came from.
    % SC_MERGESCES relabels a single-batch dataset 1/2 (original vs
    % synthetic), which overwrote the original batch IDs; IS_SYNTHETIC
    % already says which cells are new.
    batchids = [sce.c_batch_id(:); sub_sce.c_batch_id(:)];
    sce = sc_mergesces({sce, sub_sce}, "intersect", true, true);
    if numel(batchids) == sce.NumCells
        sce.c_batch_id = batchids;
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    if ~isempty(sceori)
        gui.myGuidata(FigureHandle, copy(sceori), src);   % undo the attribute
    end
    return;
end
gui.myWaitbar(FigureHandle, fw);

% Only the synthetic cells are nudged off the cells they copy. This loop
% used to add noise to every row of every stored embedding, so the real
% cells moved too.
f = fieldnames(sce.struct_cell_embeddings);
newrows = (oldcn + 1):sce.NumCells;
for i = 1:numel(f)
    s = sce.struct_cell_embeddings.(f{i});
    if ~isempty(s) && size(s, 1) == sce.NumCells && ~isempty(newrows)
        s(newrows, :) = s(newrows, :) + trnd(100, numel(newrows), size(s, 2)) * 0.75;
        sce.struct_cell_embeddings.(f{i}) = s;
    end
end

newcn = sce.NumCells;
newgn = sce.NumGenes;

if oldcn-newcn==0 && oldgn-newgn==0
        gui.myHelpdlg(FigureHandle, "No cells are synthesized.");
        return;
    end
if ~isempty(sceori)
        answer = gui.myQuestdlg(FigureHandle, ...
            sprintf('%d cells will be added.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
                newcn - oldcn, oldgn, oldcn, newgn, newcn),'', ...
                {'Accept Changes', 'Cancel Changes'}, 'Accept Changes');
        if ~strcmp(answer, 'Accept Changes')
            sce = copy(sceori);
        else
            requirerefresh = true;
        end
    else
        requirerefresh = true;
    end


gui.myGuidata(FigureHandle, sce, src);

end
