% [X,e,ker,alpha] = eels(X,Wp,Wn,l,P,ff,G[,alpha0,rho,c])
% Backtracking line search for EE
%
% Reference: procedure 3.1, p. 41ff in:
%   Nocedal and Wright: "Numerical Optimization", Springer, 1999.
%
% In:
%   X: NxL matrix, the current iterate.
%   Wp,Wn: as in ee.m.
%   l: current homotopy parameter value.
%   P: Nxd matrix, the search direction.
%   ff: value at X of the error function.
%   G: Nxd matrix, gradient at X of the error function.
%   alpha0: initial step size. Default: 1.
%   rho: rate of decrease of the step size. Default: 0.8.
%   c: Wolfe condition. Default: 1e-1.
% Out:
%   X,e,ker: new iterate, error function value and exp(-sqd).
%   alpha: step size.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.
% Copyright (c) 2012 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [X,e,ker,alpha] = eels(X,Wp,Wn,l,P,ff,G,alpha0,rho,c)
% ---------- Argument defaults ----------
if ~exist('alpha0','var') || isempty(alpha0) alpha0 = 1; end
if ~exist('rho','var') || isempty(rho) rho = 0.8; end
if ~exist('c','var') || isempty(c) c = 1e-1; end
% ---------- End of "argument defaults" ----------

alpha = alpha0;
tmp = c*G(:)'*P(:);
[e,ker] = ee_error(X+alpha*P,Wp,Wn,l);
while e > ff + alpha*tmp
  alpha = rho*alpha;
  [e,ker] = ee_error(X+alpha*P,Wp,Wn,l);
end
X = X + alpha*P;