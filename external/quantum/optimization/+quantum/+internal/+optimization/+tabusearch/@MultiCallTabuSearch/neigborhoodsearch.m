function [idxBestVar, doLocalSearch, iter] = neigborhoodsearch(L, linTerm, numVar, fbar, bestf)
%NEIGBORHOODSEARCH Find search direction
%
%   [IDXBESTVAR, DOLOCALSEARCH, ITER] = NEIGBORHOODSEARCH(METHOD, L,
%   LINTERM, NUMVAR, FBAR, BESTF) finds the next direction for
%   TabuSearch to move in. The inputs are:
%
%         L: TabuTenure
%   LINTERM: Current linear term for the model
%    NUMVAR: Number of variables in the problem
%      FBAR: Possible best function value before a call to local search
%     BESTF: Best function value found after a call to local search
%      ITER: Number of iterations
%
%   The outputs are:
%
%      IDXBESTVAR: Index of the variable to flip
%   DOLOCALSEARCH: Indicate if local search is performed
%            ITER: Number of neighborhood iterations

%   Copyright 2022-2023 The MathWorks, Inc.

doLocalSearch = false;
idxNonTabu = L == 0;
idxBestVar = idxNonTabu & (fbar + linTerm > bestf);
idxBestVar = find(idxBestVar, 1);
if isempty(idxBestVar)
    % thisLinearTerm = linTerm(idxNonTabu);
    % [~, idxMaxInNonTabu] = max(thisLinearTerm);
    % idxNonTabu = find(idxNonTabu);
    % idxBestVar = idxNonTabu(idxMaxInNonTabu);
    %
    % The following lines of code are a more efficient implementation of
    % the desired algorithm lines above
    linTerm(~idxNonTabu) = -inf;
    [~, idxBestVar] = max(linTerm);
    iter = numVar - sum(~idxNonTabu);
else
    doLocalSearch = true;
    iter = sum(idxNonTabu(1:idxBestVar));
end


