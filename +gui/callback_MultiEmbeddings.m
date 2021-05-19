function callback_MultiEmbeddings(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    if isempty(sce.struct_cell_embeddings)
        sce.struct_cell_embeddings=struct('tsne',[],'umap',[],'phate',[]);
    end    
    s1=sce.struct_cell_embeddings.tsne;
    % s2=sce.struct_cell_embeddings.umap;
    s3=sce.struct_cell_embeddings.phate;
    if ~isempty(s1) && ~isempty(s3)
        gui.sc_multiembeddings(s1,s3,'tSNE','PHATE');
    end
end