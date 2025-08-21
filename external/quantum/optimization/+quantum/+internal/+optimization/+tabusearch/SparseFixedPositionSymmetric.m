classdef SparseFixedPositionSymmetric < quantum.internal.optimization.tabusearch.FixedPositionSymmetric
%SPARSEFIXEDPOSITIONSYMMETRIC Sparse fixed position symmetric matrix
%
%   SPARSEFIXEDPOSITIONSYMMETRIC objects hold a sparse symmetric matrix for
%   use with TabuSearch. The positions of the non-zero elements are assumed
%   not to change.
%   
%   quantum.internal.optimization.tabusearch.SparseFixedPositionSymmetric(C)
%   constructs a SPARSEFIXEDPOSITIONSYMMETRIC object from the symmetric
%   matrix C. 
%
%   NOTE: It is assumed that C is symmetric and this is not checked.
%
%   Simple example:
%   
%   Given:
%
%   C =
% 
%      0     1     0     2     0
%      1     0     3     0     4
%      0     3     0     5     0
%      2     0     5     0     1
%      0     4     0     1     0
%
%   obj = quantum.internal.optimization.tabusearch.SparseFixedPositionSymmetric(C);
%
%   We have non-zero elements (transposed, just to save space in the help)
%
%   obj.RowIndex'= [2 4 1 3 5 2 4 1 3 5 2 4]
%   obj.Values'  = [1 2 1 3 4 3 5 2 5 1 4 1]
%
%   As the data is stored in compressed column form, obj.Jc(i) indicates
%   the start index in obj.Values/obj.RowIndex of the i-th column,
%   obj.Jc(i+1)-1 denotes the end index of the i-th column.
% 
%   obj.Jc' = [1 3 6 8 11 13]
%
%   Finally, obj.RowMap provides indices into values in row order, rather
%   than column order
%
%   obj.RowMap = [3 8 1 6 11 4 9 2 7 12 5 10]
%
%   As the matrix is symmetric, obj.Jc can be used to determine the start
%   and end indices of the i-th row. E.g., 
%   The third row has a start index of 6 and end index of 7. 
%   So the indices of the third row values = obj.RowMap(6:7) = [4 9]
%   The third row values = obj.Values([4 9]) = [3 5]

%   Copyright 2022-2023 The MathWorks, Inc.

    properties
        % The sparse fixed position symmetric matrix has NNZ non-zero
        % elements

        % NNZ-by-1 vector of values in the matrix
        Values

        % NNZ-by-1 vector of row indices for each non-zero value
        RowIndex

        % Start and end indices of the columns
        Jc

        % Indices into values in row order
        RowMap

        % Number of variables
        NumVariables
    end

    % Constructor
    methods

        function obj = SparseFixedPositionSymmetric(C)

            % Store the number of variables
            obj.NumVariables = size(C, 1);

            % Store sparse quadratic matrix as vectors and indices
            [obj.RowIndex, thisColumnIndex, obj.Values] = find(C);

            % Build the Jc, as you'd have it in sparse on the inside.
            % Elements in col k live in Jc(k):(Jc(k+1)-1) positions of
            % Values/RowIndex.
            obj.Jc = [1; cumsum(accumarray(thisColumnIndex,1,[size(C,2),1]))+1];

            % ColumnIndex already comes out sorted
            [~, obj.RowMap] = sort(obj.RowIndex);

        end

    end

    methods

        function obj = negateRowAndColumn(obj, idx)
            % Main task. Negate i-th row and column of a fixed position
            % symmetric sparse matrix, that is
            % C = -C(:, i);
            % C = -C(i, :)

            % Don't build the whole column index, just use the endpoints,
            % since MATLAB will very happily NOT build the whole colon object.
            c1 = (obj.Jc(idx));
            c2 = (obj.Jc(idx+1)-1);

            % Values in idx-th column
            tmp = obj.Values(c1:c2);

            % Negate column
            obj.Values(c1:c2) = -tmp;

            % Indices of thisRow in obj.Values
            rowIndex = obj.RowMap(c1:c2);

            % Negate row
            obj.Values(rowIndex) = -tmp;

        end

        function nonZeroRowIdx = getRowIndexForColumn(obj, colIdx)
            % Get the indices of the non-zero elements in a specific column

            c1 = obj.Jc(colIdx);
            c2 = obj.Jc(colIdx+1)-1;
            nonZeroRowIdx = obj.RowIndex(c1:c2);
        end

        function colVals = getColumnValues(obj, colIdx)
            % Get the values in a specific column

            c1 = obj.Jc(colIdx);
            c2 = obj.Jc(colIdx+1)-1;
            colVals = obj.Values(c1:c2);
        end

        function obj = calculateQuadraticTransform(obj, x)
            % Calculate C = C.*(1 - 2*(x - x').^2), required by the
            % TabuSearch algorithm

            rowIdx = obj.RowIndex;
            N = numel(obj.Jc)-1;
            colIdx = repelem(1:N, diff(obj.Jc));
            obj.Values = obj.Values.*(1 - 2*(x(rowIdx) - x(colIdx)).^2);
        end

        function out = sumSelectedColumns(obj, colIdx)
            % Sum across elements of selected columns

            out = zeros(obj.NumVariables, 1);
            numSpecifiedCols = numel(colIdx);
            for i = 1:numSpecifiedCols
                thisIdx = obj.Jc(colIdx(i)):obj.Jc(colIdx(i)+1)-1;
                thisRowMap = obj.RowMap(thisIdx);
                out(obj.RowIndex(thisIdx)) = out(obj.RowIndex(thisIdx)) + obj.Values(thisRowMap);
            end
        end

        function val = getValuesInRow(obj, idxRow, idxCol)
            % Get the values for the specified column indices in a specific
            % row

            % By symmetry, the row and column indices can be interchanged
            rowStart = obj.Jc(idxRow);
            rowEnd = obj.Jc(idxRow+1)-1;              
            allColumnIndex = obj.RowIndex;
            
            % We can get the values as a column
            isColVals = false;
            val = getSubsetValuesInColumn(obj, allColumnIndex, rowStart, ...
                rowEnd, idxCol, isColVals);

            % Transpose values back to a row
            val = val';

        end

        function val = getValuesInColumn(obj, idxRow, idxCol)
            % Get the values for the specified row indices in a specific
            % column

            colStart = obj.Jc(idxCol);
            colEnd = obj.Jc(idxCol+1)-1;  
            isColVals = true;
            val = getSubsetValuesInColumn(obj, obj.RowIndex, colStart, ...
                colEnd, idxRow, isColVals);

        end

        function val = getSubsetValuesInColumn(obj, nonZeroRowIndex, ...
                colStart, colEnd, idxRow, isColVals)

            % Initialize values
            val = zeros(numel(idxRow), 1);

            % Actual row indices of the non-zero row entries in the idxCol-th column
            thisRowIdx = nonZeroRowIndex(colStart:colEnd);

            % Find non-zero row indices in idxRow
            [~, idxRowVals1, idxRowVals2] = intersect(thisRowIdx, idxRow);
            nonZeroValsColIdx = colStart:colEnd;
            idxRemainingVals = nonZeroValsColIdx(idxRowVals1);

            % Extract non-zero values
            if isColVals
                nonZeroVals = obj.Values(idxRemainingVals);
            else
                nonZeroVals = obj.Values(obj.RowMap(idxRemainingVals));
            end

            % Insert it into the column
            val(idxRowVals2) = nonZeroVals;

        end

        function fixedStruct = createFixedStruct(obj)
            % Create a structure containing information that is fixed
            % and independent of the values in the matrix. Note that the
            % indices are shifted to be zero-based.

            % TODO: Convert to string
            fixedStruct.Type = 'sparse';
            fixedStruct.NumVariables = obj.NumVariables;
            fixedStruct.Ir = obj.RowIndex - 1;
            fixedStruct.Jc = obj.Jc - 1;
            fixedStruct.RowMap = obj.RowMap - 1;
            fixedStruct.Nnz = fixedStruct.Jc(end);

        end
        
        

    end
end
