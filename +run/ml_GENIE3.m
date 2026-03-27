function [A, s, G] = ml_GENIE3(X, genelist, donorm, plotit)
% GENIE3: GEne Network Inference with Ensemble of trees
%
% Thin wrapper — implementation moved to net.genie3.
% See also: net.genie3, sc_grn

if nargin < 2, genelist = []; end
if nargin < 3, donorm = false; end
if nargin < 4, plotit = false; end

[A, s, G] = net.genie3(X, genelist, donorm, plotit);
end
