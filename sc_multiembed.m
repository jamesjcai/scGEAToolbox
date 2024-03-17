function [sce] = sc_multiembed(sce, dometaviz)
if nargin<2, dometaviz=false; end
if ~isa(sce, 'SingleCellExperiment')
    error('sce should be a SingleCellExperiment object.'); 
end
methodtag = {'tsne2d', 'tsne3d', 'umap2d', 'umap3d',...
    'phate2d', 'phate3d', 'metaviz2d', 'metaviz3d'};


validmethodtag = fieldnames(sce.struct_cell_embeddings);
assert(all(ismember(methodtag,validmethodtag)));

if ~dometaviz, methodtag = methodtag(1:6); end
forced = true;
usehvgs = true;
K = 2000;
for k=1:length(methodtag)
    ndim = 2 + contains(methodtag{k},'3d');
    try
        sce = sce.embedcells(methodtag{k}, ...
            forced, usehvgs, ndim, K, [], false);
    catch ME
        disp(ME.message);
    end
end
end
