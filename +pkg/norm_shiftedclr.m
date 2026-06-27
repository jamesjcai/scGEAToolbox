function Xclr = norm_shiftedclr(X)
% PFlog1pPF
% ref: https://www.biorxiv.org/content/10.1101/2022.05.06.490859v3
% X is genes x cells


% PF normalization
d = full(sum(X,1));                     % cell depths
sf = mean(d) ./ d;

Xpf = X * spdiags(sf(:),0,length(sf),length(sf));

% log1p
Xpf(Xpf~=0) = log1p(Xpf(Xpf~=0));

% CLR centering
cell_mean = mean(Xpf,1);

Xclr = full(Xpf) - cell_mean;



% Xclr = PFlog1pPF(X, 1);





end

function Z = PFlog1pPF(X, c)
% PFlog1pPF  Shifted CLR normalization (Booeshaghi et al. 2022)
% Input X: genes x cells (sparse or dense)
% c: pseudocount (default 1); equivalent to scale factor K = 1/c

if nargin < 2, c = 1; end
if ~issparse(X), X = sparse(double(X)); end

% Step 1: PF — scale each cell to mean depth
%         (mean(d) cancels in CLR centering but improves numerical range)
d    = full(sum(X, 1));              % 1×N
sf   = mean(d) ./ d;                % 1×N scale factors
U    = X * spdiags(sf(:), 0, numel(sf), numel(sf));   % sparse G×N

% Step 2: log(u + c) — only on nonzero entries (sparse-safe)
%         For c=1: log1p on proportions; zeros map to log(c) after
%         centering, which is correct per Proposition 1.3
if c == 1
    [ii, jj, vv] = find(U);
    vv = log1p(vv);                 % log(u_i + 1) for nonzeros
    V  = sparse(ii, jj, vv, size(U,1), size(U,2));
    % NB: zero entries contribute log1p(0)=0, consistent with sparse storage
else
    V = full(U);
    V = log(V + c);                 % must densify for c != 1
end

% Step 3: CLR centering — subtract per-cell mean (second PF in log space)
%         mean() over all G genes including structural zeros
cell_mean = full(mean(V, 1));       % 1×N
Z = full(V) - cell_mean;            % G×N dense; sum(Z,1) == 0 per cell

end