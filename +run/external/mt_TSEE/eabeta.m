% [b,W] = eabeta(d2,b0,logK,B) Gaussian EAs: beta and affinities
%
% Computes the values of beta and the corresponding Gaussian affinities for
% one point. It does root finding with Newton's method embedded in a bisection
% loop to ensure global convergence.
%
% In:
%   d2: 1 x k vector of square distances to the k nearest neighbors.
%   b0: initial value of beta.
%   logK: log of the perplexity K. 
%   B: 1x2 vector of lower and upper bounds on log(beta).
% Out:
%   b: log(beta) value.
%   W: 1 x k vector of EAs computed with this beta. 

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [b,W] = eabeta(d2,b0,logK,B)

% No need to change these ever, probably:
maxit = 20;		% Max. no. iterations before a bisection
tol = 1e-10;		% Convergence tolerance to stop iterating

if (b0<B(1) || b0>B(2)) b=(B(1)+B(2))/2; else b=b0; end

i = 1;			% maxit counter
% Inside the loop, pbm is a flag to detect any numerical problems (zero
% gradient, infinite function value, etc.).
while 1
  bE = exp(b); pbm = 0;
  
  % Compute the function value: m0, m1v, m1 are O(N)
  ed2 = exp(-d2*bE); m0 = sum(ed2);
  if m0<realmin		% Numerical error
    e = -logK; pbm = 1;
  else
    m1v = ed2.*d2/m0; m1 = sum(m1v); e = bE*m1 + log(m0) - logK;
  end
  
  if abs(e) < tol break; end
  
  % Very narrow bounds, no need to iterate. This can happen if K is very small.
  if B(2)-B(1) < 10*eps break; end
  
  % Update the bounds
  if (e<0 && b<=B(2))
    B(2) = b;
  elseif (e>0 && b>=B(1))
    B(1) = b;
  end
  % on the right term, the first = 0 means the denominator of pn is 0;
  % the second = 0 means e is an infinite value, which means at least one
  % probability is 0
  % the third = 0 means the entropy is smaller than 0, which cannot hold.
  % the forth = 0 means the entropy of a distribution on k poins is larger
  % than the value of log(k), which cannot hold neither. ---AN
  pbm = pbm || isinf(e) || e<-logK || e>log(length(d2))-logK;
  
  if ~pbm
    if i==maxit		% Exceeded maxit, bisection step
      b = (B(1)+B(2))/2; i=1; continue;
    end
    % Compute the gradient: m2 is O(N)
    eg2 = bE^2; m2 = m1v*d2'; m12 = m1^2-m2; 
    g = eg2*m12;
    if g==0 pbm=1; end
  end
  
  % If there was a problem with the function, do bisection with old bounds.
  % If there was a problem with the gradient, do bisection with new bounds.
  if pbm 
    % If there is a numerical problem on both ends of the bounds, return 
    % current value. Practically this happens only when K is very small.
    esqd1 = exp(-d2*exp(B(1))); esqd2 = exp(-d2*exp(B(2)));
    if sum(esqd1+esqd2) < 2*sqrt(realmin) break; end
    b = (B(1)+B(2))/2; i=1; continue; 
  end
  
  % Newton step ok, update bounds
  p = -e/g; b = b + p;
  
  if (b<B(1) || b>B(2))	% Out of bounds, bisection step
    b = (B(1)+B(2))/2; i = 0;
  end
  i=i+1;
end

W = ed2/m0;		% Affinities
