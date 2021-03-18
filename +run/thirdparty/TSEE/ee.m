% [X,E,A,T] = ee(Wp,Wn,d,l[,opts])
%
% Fast training of nonlinear embeddings using the spectral direction for the
% Elastic Embedding (EE) algorithm.
%
% In:
%   Wp, Wn: NxN positive and negative weight matrices.
%   d: dimension of the low-dim space (integer >= 1).
%   l: Hx1 list of homotopy parameters in increasing order.
%   opts: structure of optional parameters, with the following fields:
%     X0: N x d matrix, initial low-dim coordinates. Default: random.
%     pattern: NxN binary sparse matrix (desired sparsity pattern in the
%       spectral direction). Default: none.
%     tol: minimum relative distance between consecutive X. Default: 1e-3.
%     maxit: maximum number of iterations (for each value of l). Default: 100.
%     runtime: maximum total allowed runtime in seconds. Default: Inf.
% Out:
%   X, E, A, T: (H+1)x1 (cell) array where element i+1 corresponds to the ith
%      l value of the low-dim coordinates, error curve, step sizes and runtime
%      per iteration, respectively.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2012 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [X,E,A,T] = ee(Wp,Wn,d,l,opts)

N = size(Wp,1);
% ---------- Argument defaults ----------
if ~exist('opts','var') || isempty(opts) opts = []; end;

if ~isfield(opts,'X0') opts.X0 = 1e-5*randn(N,d); end;
if ~isfield(opts,'pattern') opts.pattern = []; end;
if ~isfield(opts,'tol') || isempty(opts.tol) opts.tol = 1e-3; end
if ~isfield(opts,'maxit') || isempty(opts.maxit) opts.maxit = 200; end
if ~isfield(opts,'runtime') || isempty(opts.runtime) opts.runtime = Inf; end
% ---------- End of "argument defaults" ----------

tic
% make sure that the data are symmetric and have zeros along the diagonal
Wp = (Wp+Wp')/2; Wp(1:N+1:N^2) = 0;
Wn = (Wn+Wn')/2; Wn(1:N+1:N^2) = 0;

% compute graph Laplacian
Dp = diag(sparse(sum(Wp,2)));
Lp4 = 4*(Dp-Wp);

% smallest nonzero element of diag(Lp4)
mDiagLp = Lp4(find(Dp>0,1,'first'));

% compute Cholesky factor of the (sparsified) graph Laplacian
if ~isempty(opts.pattern)
  [R,~,S] = chol(sparse(opts.pattern.*Lp4) + 1e-10*mDiagLp*speye(N,N),'upper');
else
  R = chol(Lp4 + 1e-10*mDiagLp*speye(N,N),'upper');
  S = speye(N,N);
end
St = S'; Rt = R';

Xold = opts.X0;
for i = 1:length(l)
  [e(1), ker] = ee_error(Xold,Wp,Wn,l(i));
  j = 1; a(1) = 1; t(1) = toc;
  convcrit = (opts.maxit>=1) & (t(1)<opts.runtime);
  
  while convcrit
    WWn = l(i)*Wn.*ker; DDn = diag(sparse(sum(WWn,2)));
    G = (Lp4-4*(DDn-WWn))*Xold;   % gradient
    P = -S*(R\(Rt\(St*G))); 	    % spectral direction
    [X,e(j+1),ker,a(j+1)] = eels(Xold,Wp,Wn,l(i),P,e(j),G,a(j)); % line search
    
    convcrit = (j<opts.maxit) & (t(end)<opts.runtime) & ...
      (norm(X-Xold,'fro')>opts.tol*norm(Xold,'fro'));

    Xold = X;
    j = j+1;
    t(j) = toc;
  end
      
  % save parameters in an array (non homotopy) or in a cell array (for homotopy)
  if(length(l)==1)
    E = e; A = a; T = t;
  else
    E{i} = e; A{i} = a; T{i} = t;
  end
end