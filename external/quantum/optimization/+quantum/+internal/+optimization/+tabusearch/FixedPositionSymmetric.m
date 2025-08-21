classdef (Abstract) FixedPositionSymmetric
%FIXEDPOSITIONSYMMETRIC Abstract fixed position symmetric matrix
%
%   FIXEDPOSITIONSYMMETRIC is an abstract base class for fixed position
%   NumVariables-by-NumVariables symmetric matrices used in the TabuSearch
%   algorithm.

%   Copyright 2022-2023 The MathWorks, Inc.

    properties (Abstract)
        NumVariables
    end

    methods (Abstract)
        
        % Method to calculate C = C.*(1 - 2*(x - x').^2), required by the
        % TabuSearch algorithm
        calculateQuadraticTransform

        % Get the values in a specific column
        getColumnValues

        % Get the indices of the non-zero elements in a specific column
        getRowIndexForColumn

        % Get the values for the specified row indices in a specific column
        getValuesInColumn

        % Get the values for the specified column indices in a specific row
        getValuesInRow

        % Negate the i-th row and column of M, that is,
        % C = -C(:, i)
        % C = -C(i, :)
        negateRowAndColumn

        % Sum across elements of selected columns
        sumSelectedColumns

        % Create structure of fixed information that is independent of the
        % values in the matrix
        createFixedStruct

    end

end