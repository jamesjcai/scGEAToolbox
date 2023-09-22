function xcallback_MultiEmbeddingViewer(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
if isempty(sce.struct_cell_embeddings)
    sce.struct_cell_embeddings = struct('tsne', [], 'umap', [], 'phate', []);
end
s1 = sce.struct_cell_embeddings.tsne;
s2 = sce.struct_cell_embeddings.umap;
s3 = sce.struct_cell_embeddings.phate;
if ~isempty(s1) && ~isempty(s2)
    gui.sc_multiembeddings(s1, s2, 'tSNE', 'UMAP');
elseif ~isempty(s1) && ~isempty(s3)
    gui.sc_multiembeddings(s1, s3, 'tSNE', 'PHATE');
elseif ~isempty(s2) && ~isempty(s3)
    gui.sc_multiembeddings(s2, s3, 'UMAP', 'PHATE');
else
    warndlg('This function requires two embeddings (tSNE, UMAP or PHATE) precomputed.');
    end
end
