function x = initialPoint(numVar)
%INITIALPOINT Generates a random initial point for Tabu Search

%X = INITIALPOINT(NUMVAR) generates a random vector of length numVar

%   Copyright 2022 The MathWorks, Inc.

        x = randi([0 1], [numVar 1]);

end

