pyrun("from scimilarity import CellAnnotation")
pyrun("ca = CellAnnotation(model_path=a)", a='y:\jcai\models\model_v1.1');
% pyrun("ca = CellAnnotation(model_path = 'y:\\jcai\\models\\model_v1.1')")
targetg = string(pyrun("a = ca.gene_order","a"))';
load test\testdata.mat
[y, idx] = ismember(sce.g, targetg);
X = zeros(numel(targetg), size(sce.X, 2));
idx = idx(y);
x = full(sce.X(y,:));
X(idx,:) = x;
X = single(log1p(pkg.norm_libsize(X, 10000)));

pyrun("from scipy.sparse import csr_matrix")
pyrun("X = csr_matrix(a)", a=X');
embeddings = pyrun("embeddings = ca.get_embeddings(X)", "embeddings");
    
% pyrun("predictions = ca.get_predictions_knn(embeddings)")
