function e_writeh5(X,genelist,filename)

% if isa(X,'SingleCellExperiment')
%     genelist=X.g;
%     X=X.X;
% end
% if nargin<2, genelist=string([1:size(X,1)].'); end
if issparse(X)
    X=full(X);
end
h5create(filename, '/X', size(X));
h5write(filename, '/X', X);

if ~isempty(genelist)
    h5create(filename, '/g', size(genelist),'Datatype','string');
    h5write(filename, '/g', genelist);
end
