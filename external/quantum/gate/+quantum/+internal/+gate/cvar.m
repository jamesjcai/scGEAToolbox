function out = cvar(p, v, alpha)
% Internal use only. Returns the Conditional Value-at-Risk for probability
% vector p with values v using the confidence level 0 < alpha <= 1.
%
%   Copyright 2025 The MathWorks, Inc.

% References:
% [1] Barkoutsos, Panagiotis Kl, et al. "Improving variational quantum 
% optimization using CVaR." Quantum 4 (2020).

p = p(:);
v = v(:);

if alpha==1
    out = p.'*v;
    return
end

[sv, idx] = sort(v);
sp = p(idx);

cdf = cumsum(sp);
cut = find(cdf >= alpha, 1);

% Probability at the cut index is adjusted such that the total probability
% equals alpha, i.e., sum(w)=alpha.
w = sp;
w(cut) = alpha - cdf(cut) + sp(cut);

% Return weighted average over sorted values using the normalized
% probabilities
out = (w(1:cut) / alpha).' * sv(1:cut);

end