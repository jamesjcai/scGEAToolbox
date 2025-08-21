classdef DenseFixedPositionSymmetric < quantum.internal.optimization.tabusearch.FixedPositionSymmetric
%DENSEFIXEDPOSITIONSYMMETRIC Dense fixed position symmetric matrix
%
%   DENSEFIXEDPOSITIONSYMMETRIC objects hold a dense symmetric matrix for
%   use with TabuSearch. The positions of the non-zero elements are assumed
%   not to change.
%   
%   quantum.internal.optimization.tabusearch.DenseFixedPositionSymmetric(C)
%   constructs a DENSEFIXEDPOSITIONSYMMETRIC object from the symmetric
%   matrix C. 
%
%   NOTE: It is assumed that C is symmetric and this is not checked.

%   Copyright 2022-2023 The MathWorks, Inc.

    properties
        NumVariables
        Values
    end

    % Constructor
    methods

        function obj = DenseFixedPositionSymmetric(C)

            % Store dense symmetric matrix
            obj.Values = C;

            % Store the number of variables
            obj.NumVariables = size(C, 1);

        end

    end

    % Main task. Negate i-th row and column of a symmetric sparse matrix.
    % Do this alot.
    methods

        function obj = negateRowAndColumn(obj, idx)
            % Main task. Negate i-th row and column of a fixed position
            % symmetric sparse matrix, that is
            % C = -C(:, i);
            % C = -C(i, :)

            % Values in i-th row and column
            tmp = obj.Values(:, idx);

            % Negate column
            obj.Values(:, idx) = -tmp;

            % Negate row
            obj.Values(idx, :) = -tmp;

        end

        function colVals = getColumnValues(obj, colIdx)
            % Get the values in a specific column
            
            colVals = obj.Values(:, colIdx);
        end

        function obj = calculateQuadraticTransform(obj, x)
            % Calculate C = C.*(1 - 2*(x - x').^2), required by the
            % TabuSearch algorithm
            
            obj.Values = obj.Values.*(1 - 2*(x - x').^2);
        end

        function out = sumSelectedColumns(obj, colIdx)
            % Sum across elements of selected columns
            
            out = sum(obj.Values(:, colIdx), 2);
        end
        
        function nonZeroRowIdx = getRowIndexForColumn(obj, ~)
            % Get the indices of the non-zero elements in a specific column
            
            nonZeroRowIdx = 1:obj.NumVariables;
        end

        function val = getValuesInRow(obj, idxRow, idxCol)
            % Get the values for the specified column indices in a specific
            % row
            
            val = obj.Values(idxRow, idxCol);
        end

        function val = getValuesInColumn(obj, idxRow, idxCol)
            % Get the values for the specified row indices in a specific
            % column
            
            val = obj.Values(idxRow, idxCol);
        end

        function fixedStruct = createFixedStruct(obj)
            % Create a structure containing information that is fixed
            % and independent of the values in the matrix

            % TODO: Convert to strings
            fixedStruct.Type = 'dense';
            fixedStruct.NumVariables = obj.NumVariables;

        end

    end

end