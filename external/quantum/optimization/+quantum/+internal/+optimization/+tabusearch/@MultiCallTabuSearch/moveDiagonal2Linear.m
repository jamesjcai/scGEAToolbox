function [C, d] = moveDiagonal2Linear(C)
%MOVEDIAGONAL2LINEAR Move diagonal terms to linear coefficient vector
%
%   [C, D] = MOVEDIAGONAL2LINEAR(C) moves the diagonal terms of the
%   QUBO, C, to the linear term, D,

%   Copyright 2022 The MathWorks, Inc.

% Palubeckis assumes the diagonal elements of C are zero. Move
% any non-zero diagonal elements into d and zero out the
% diagonal of C
numVar = size(C, 1);
d = diag(C);
C(1:numVar+1:end) = 0;

end
