function [reduction] = run_umap_lite_super(inData, k)
    umap = UMAP;
    umap.n_components = k;
    umap.setMethod('MEX');
    reduction = umap.fit_transform(inData);
end
