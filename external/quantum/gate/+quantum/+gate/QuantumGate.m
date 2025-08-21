classdef QuantumGate < matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
    %QUANTUMGATE Abstract superclass for SimpleGate and CompositeGate
    %
    % This class is for internal use only and will change in a future release.
    % Do not use this class.
    %
    %   See also quantum.gate.SimpleGate, quantum.gate.CompositeGate

    %   Copyright 2021-2022 The MathWorks, Inc.

    % Currently some methods of quantumCircuit are just accessing internal
    % properties of the gates. That should be replaced with additional
    % methods, which will be added here, too.

    methods
        function M = getMatrix(obj)
            %GETMATRIX  Unitary matrix representation of a gate.
            %
            %   M = getMatrix(G) returns the unitary matrix representing
            %   gate G. This matrix has size 2^numQubits -by- 2^numQubits,
            %   where numQubits is the largest qubit index in G.
            %
            %   Example:
            %       % 2-by-2 matrix representation of a gate
            %       g = hGate(1);
            %       M = getMatrix(g)
            %
            %       % 4-by-4 matrix representation of a gate
            %       g = hGate(2);
            %       M = getMatrix(g)
            %
            %       % Matrix representation of a quantum circuit
            %       g = hGate(1);
            %       M = getMatrix(quantumCircuit(g, 2))
            %
            %       % Matrix representation of a composite gate
            %       circ = quantumCircuit([hGate(1) cxGate(1,2)]);
            %       cg = compositeGate(circ, 2:3);
            %       M = getMatrix(cg)
            %
            %   See also quantumCircuit/getMatrix, quantum.gate.SimpleGate,
            %   quantum.gate.CompositeGate,
            %
            arguments
                obj (1,1) quantum.gate.QuantumGate
            end
            M = getMatrix(quantumCircuit(obj));
        end
    end

    methods (Hidden)
        function plot(varargin)
            error(message('quantum:gates:PlotNotSupported'));
        end
    end

    methods (Static, Access = protected)
        function g = getDefaultScalarElement()
            g = quantum.gate.SimpleGate.MissingGate();
        end
    end

    methods (Abstract, Access = protected)
        str = getOneLineDisplayString(gate)
    end

    methods (Sealed, Access = protected)
        function displayNonScalarObject(gates)
            if isvector(gates)
                fprintf( '  %s\n\n', getVectorHeader(gates));
                displayOneLines(gates)
            else
                fprintf( '  %s\n', getArrayHeader(gates))
            end
            f = format;
            if f.LineSpacing == "loose"
                fprintf(newline)
            end
        end

        function displayEmptyObject(gates)
            fprintf( '  %s\n', getArrayHeader(gates))
            f = format;
            if f.LineSpacing == "loose"
                fprintf(newline)
            end
        end
    end

    methods (Sealed, Access = private)
        function header = getVectorHeader(gates)
            % getVectorHeader   Return the header to be displayed for a
            % vector of gates
            sizeString = join(string(size(gates)), ...
                matlab.internal.display.getDimensionSpecifier);
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(gates);
            header = getString(message('quantum:gates:VectorHeader', sizeString, className));
        end

        function header = getArrayHeader(gates)
            % getArrayHeader   Return the header to be displayed for an
            % array of gates of size >= 2
            sizeString = join(string(size(gates)), ...
                matlab.internal.display.getDimensionSpecifier);
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(gates);
            header = getString(message('quantum:gates:ArrayHeader', sizeString, className));
        end

        function displayOneLines(gates)
            % displayOneLines   Display one line for each gate in the
            % vector

            c = cell(length(gates), 6);
            for ii=1:length(gates)
                c(ii, :) = [{string(ii)} getOneLineDisplayString(gates(ii))];
            end
            c(:, [3 4 6]) = iMakeString(c(:, [3 4 6]));

            c = reshape([c{:}], size(c));

            c = [["Id", "Gate", "Control", "Target", "Angle", "NumGates"]; c];

            % Remove "Angle" and/or "NumGates" column if it has no entries
            for jj=[6 5]
                if all(c(2:end, jj) == "")
                    c(:, jj) = [];
                end
            end

            for ii=1:size(c, 2)
                c(:, ii) = pad(c(:, ii));
            end
            
            c(1,:) = "<strong>"+c(1,:)+"</strong>";
            fprintf( '       %s\n', strjoin(c(1, :), "   "));

            for idx=1:numel(gates)
                fprintf( '       %s\n', strjoin(c(idx+1, :), "   "));
            end
        end
    end
end


function vec = iMakeString(vec)
for ii=1:numel(vec)
    vec{ii} = string(vec{ii});
    if isempty(vec{ii})
        vec{ii} = "";
    elseif ~isscalar(vec{ii})
        vec{ii} = "[" + strjoin(vec{ii}, ",") + "]";
    end
end
end
