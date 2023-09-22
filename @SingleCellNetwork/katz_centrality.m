function s = katz_centrality(obj, a, b)
%Katz centrality (Katz Status Index)
% The Katz centrality for node i is:
% x(i)=alpha * sum(A(ij)*x(j), j) + beta
% where A is the adjacency matrix of the graph G with eigenvalues lambda. The parameter beta controls the initial centrality and alpha < 1/lambda(max).
%ref: https://rdrr.io/cran/centiserve/man/katzcent.html
%input: A - The alpha parameter, which must be between 0.0-0.2. The default is 0.1.
% attenuation factor A
%     alpha=0.1
%     beta=1.0,
%     max_iter=1000,
%     tol=1.0e-6,
%     nstart=None,
%     normalized=True,
%     weight=None,
% ref: https://networkx.org/documentation/stable/_modules/networkx/algorithms/centrality/katz.html

if nargin < 2, a = 0.1; end
if ~(a >= 0 && a <= 0.2), error('alpha parameter'); end
if nargin < 3, b = 1.0; end
A = obj.A;
leadingeigenvalue = max(eig(A));
n = size(A, 1);
a = 0.9 * (1 / max(eig(A)));
s = b * (inv(eye(n)-a*A')) * ones(n, 1);
end
