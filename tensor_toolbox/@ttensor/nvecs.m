function u = nvecs(X, n, r, opts)
%NVECS Compute the leading mode-n vectors for a ttensor.
%
%   U = NVECS(X,n,r) computes the r leading eigenvalues of Xn*Xn'
%   (where Xn is the mode-n matricization of X), which provides
%   information about the mode-n fibers. In two-dimensions, the r
%   leading mode-1 vectors are the same as the r left singular vectors
%   and the r leading mode-2 vectors are the same as the r right
%   singular vectors.
%
%   U = NVECS(X,n,r,OPTS) specifies options:
%   OPTS.eigsopts: options passed to the EIGS routine [struct('disp',0)]
%   OPTS.flipsign: make each column's largest element positive [true]
%
%   Examples
%   X = ttensor(tensor(ones(2,3,4)), 2*ones(1,2), 3*ones(2,3), 4*ones(4,4));
%   nvecs(X, 2, 2)
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','nvecs_doc.html')))">Documentation page for n-vecs</a>
%
%   See also TTENSOR, TENMAT, EIGS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if ~exist('opts', 'var')
    opts = struct;
end

if isfield(opts, 'eigsopts')
    eigsopts = opts.eigsopts;
else
    eigsopts.disp = 0;
end

% Compute inner product of all n-1 factors
V = cell(ndims(X), 1);
for i = 1:ndims(X)
    if i == n
        V{i} = X.u{i};
    else
        V{i} = X.u{i}' * X.u{i};
    end
end

% Form H
H = ttm(X.core, V);

if isa(H, 'sptensor')
    HnT = double(sptenmat(H, n, 't'));
else
    H = full(H);
    HnT = double(tenmat(H, n, 't'));
end
G = X.core;
if isa(G, 'sptensor')
    GnT = double(sptenmat(G, n, 't'));
else
    G = full(G);
    GnT = double(tenmat(G, n, 't'));
end

% Compute Xn * Xn'
Y = HnT' * GnT * X.u{n}';

[u, ~] = eigs(Y, r, 'LM', eigsopts);

if isfield(opts, 'flipsign')
    flipsign = opts.flipsign;
else
    flipsign = true;
end

if flipsign
    % Make the largest magnitude element be positive
    [~, loc] = max(abs(u));
    for i = 1:r
        if u(loc(i), i) < 0
            u(:, i) = u(:, i) * -1;
        end
    end
end
