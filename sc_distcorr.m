function sc_distcorr(X, g, gset)

Xm = sc_impute(X, "type", "MAGIC");
idx = ismember(g, gset);
x = Xm(idx, :);
D = pkg.e_distcorrmtx(x);     % nxn matrix, n genes
R = corr(x');
mdl = fitlm(R(:), D(:));
res = mdl.Residuals.Raw;
res = reshape(res, size(R));


