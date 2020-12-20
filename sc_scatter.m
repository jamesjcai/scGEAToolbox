function sc_scatter(X,genelist,s,c)

if isa(X,'SingleCellExperiment')
    sc_scatter_sce(X);
    return;
end
if nargin<4 || isempty(c), c=ones(size(X,2),1); end
if nargin<3 || isempty(s), s=randn(size(X,2),3); end
if nargin<2 || isempty(genelist)
    genelist=string((1:size(X,1))');
end

sce=SingleCellExperiment(X,genelist,s,c);
sc_scatter_sce(sce);

end
