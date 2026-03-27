function A = distcorrnet(X)
% Construct co-expression network using pairwise distance correlation
%
% A = net.distcorrnet(X)
%
% X - genes x cells expression matrix
% A - genes x genes distance correlation matrix
%
% Distance correlation detects nonlinear dependence (symmetric).
% ref: Szekely et al., Ann. Statist. 2007.

arguments
    X double
end

X = double(full(X));
A = pkg.e_distcorrmtx(X);
end
