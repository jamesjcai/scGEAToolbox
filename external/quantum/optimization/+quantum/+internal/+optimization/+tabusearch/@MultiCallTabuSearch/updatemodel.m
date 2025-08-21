function [Cp, dp] = updatemodel(Cp, dp, idx)
%UPDATEMODEL Update the Palubeckis form of the QUBO model
%
%   [CP, DP] = UPDATEMODEL(CP, DP, IDX) updates the Paulbeckis' form of the
%   QUBO model. The model is updated by mapping the IDX-th variable to
%   zero.

%   Copyright 2022 The MathWorks, Inc.

% Negate linear term for this column
dp(idx) = -dp(idx);

% Get column
tmp = getColumnValues(Cp, idx);

% Get row indices for this column
thisRowIdx = getRowIndexForColumn(Cp, idx);

% Update linear term
dp(thisRowIdx) = dp(thisRowIdx) + 2*tmp;

% Negate row and coulmn in quadratic term
Cp = negateRowAndColumn(Cp, idx);