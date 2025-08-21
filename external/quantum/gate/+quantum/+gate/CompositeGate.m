classdef (Sealed, ...
        InferiorClasses = {?matlab.ui.Figure, ?matlab.ui.container.Tab, ?matlab.ui.container.Panel, ...
        ?matlab.graphics.axis.Axes, ?matlab.ui.control.UIAxes, ...
        ?matlab.ui.container.GridLayout, ?matlab.graphics.layout.TiledChartLayout}) ...
        CompositeGate  < quantum.internal.gate.gates.Gate & quantum.gate.QuantumGate
    %COMPOSITEGATE Composite gate for quantum computing
    %
    %   Use compositeGate creation function to create a
    %   quantum.gate.CompositeGate object. This object fulfills the purpose
    %   of a subfunction in classical programming: It contains a set of
    %   inner gates acting on a small set of qubits, and a mapping from
    %   this small set of qubits to the qubits of the circuit that contains
    %   this composite gate.
    %
    %   Example:
    %       % Construct a CompositeGate from a quantumCircuit:
    %       innerGates = [hGate(1); cxGate(1, 2)];
    %       innerCirc = quantumCircuit(innerGates, Name="bell");
    %       outerGates = [hGate(1:4)
    %                     compositeGate(innerCirc, [1 3])
    %                     compositeGate(innerCirc, [2 4])];
    %       outerCirc = quantumCircuit(outerGates);
    %       plot(outerCirc)
    %       % Click on a CompositeGate block in the plot - a new figure
    %       % showing the internal gates of the block will appear.
    %
    %   CompositeGate properties:
    %      Name             - Name of the composite gate.
    %      TargetQubits     - Qubits of outer circuit that the inner gates are mapped to.
    %      ControlQubits    - Control qubits of the composite gate, which are always empty.
    %      Gates            - Array of inner gates.
    %
    %   CompositeGate methods:
    %      plot             - Plot a CompositeGate.
    %      inv              - Return inverse of a CompositeGate.
    %      getMatrix        - Unitary matrix representation of a CompositeGate.
    %
    %   CompositeGate creation functions:
    %      compositeGate       Construct a Composite gate
    %      qftGate             Quantum Fourier transform gate
    %      mcxGate             Multi-controlled X gate
    %
    %   See also quantumCircuit, quantum.gate.SimpleGate

    %   Copyright 2021-2024 The MathWorks, Inc.

    properties
        %NAME - Name of the composite gate
        %   Name is a string giving the name of the composite gate.
        %   Default is "".
        %   This is also used as a label in the circuit plot of the circuit
        %   block containing the composite gate.
        %
        %   See also quantumCircuit, quantum.gate.CompositeGate
        Name(1, 1) string = "";
    end

    properties (SetAccess = private, GetAccess = public)
        %CONTROLQUBITS - Control qubits of the composite gate
        %   ControlQubits of the CompositeGate, which are always an empty
        %   row vector.
        %
        %   See also quantumCircuit, quantum.gate.CompositeGate
        ControlQubits = zeros(1, 0)
    end

    properties (Dependent, SetAccess = private, GetAccess = public)
        %TARGETQUBITS - Qubits that the composite gate acts on
        %   TargetQubits is a row vector of qubit indices, indicating
        %   which qubits the gate acts on. This is a mapping from the inner
        %   qubit indices used by the Gates in the CompositeGate, to the
        %   outer qubit indices of the circuit containing the
        %   CompositeGate.
        %
        %   See also quantumCircuit, quantum.gate.CompositeGate
        TargetQubits
    end

    properties (Dependent, SetAccess = private, GetAccess = public)
        %GATES - Array of inner gates
        %   GATES is a vector containing all the inner gates of the
        %   composite gate. The elements of this vector are of type
        %   SimpleGate or CompositeGate.
        %
        %   The qubits used in the Gates array must not be larger than
        %   length(TargetQubits).
        %
        %   See also quantumCircuit, quantum.gate.CompositeGate
        Gates
    end

    properties(Dependent, Hidden, SetAccess = private, GetAccess = public)
        nameForPlot;
    end

    properties(Access = private)
        TargetQubits_
        Gates_ = quantum.gate.QuantumGate.empty(0, 1);
    end

    methods
        function tqb = get.TargetQubits(obj)
            tqb = obj.TargetQubits_;
        end

        function gates = get.Gates(obj)
            gates = obj.Gates_;
        end

        function obj = set.Gates(obj, gates)
            quantum.internal.gate.checkGates(gates, length(obj.TargetQubits_));
            obj.Gates_ = reshape(gates, [], 1);
        end

        function obj = set.Name(obj, new_name)
            new_name = string(new_name);
            quantum.internal.gate.checkName(new_name);
            obj.Name = new_name;
        end
    end

    methods(Access=protected)
        function str = getOneLineDisplayString(obj)
            if strlength(obj.Name) == 0
                str = {"composite", obj.ControlQubits, obj.TargetQubits, ...
                    "", numel(obj.Gates)};
            else
                str = {obj.Name, obj.ControlQubits, obj.TargetQubits, ...
                    "", numel(obj.Gates)};
            end
        end
    end

    methods
        function obj = CompositeGate(c, targetQubits, NameValueArgs)
            arguments
                c
                targetQubits {mustBeTargetQubits}
                NameValueArgs.Name {mustBeTextScalar}
            end

            if isa(c, 'quantum.gate.QuantumGate')
                c = quantumCircuit(c, length(targetQubits));
            elseif ~isa(c, 'quantumCircuit')
                error(message('quantum:CompositeGate:invalidFirstInput'));
            end

            if length(targetQubits) ~= c.NumQubits
                error(message('quantum:CompositeGate:InvalidNumQubits',length(targetQubits), c.NumQubits));
            end

            if isfield(NameValueArgs, 'Name')
                c.Name = NameValueArgs.Name;
            end

            obj.TargetQubits_ = targetQubits(:)';
            obj.Gates = c.Gates;
            obj.Name = c.Name;
        end

        function n = get.nameForPlot(obj)
            if strlength(obj.Name) == 0
                n = "CG";
            else
                n = obj.Name;
            end
        end

        function obj = inv(obj)
            %INV  Return inverse of CompositeGate.
            %   invcg = INV(cg) returns a new CompositeGate invcg with
            %   inverse behavior of cg. This reverses the order of the
            %   gates and replaces each original gate with its inverse.
            %
            %   Example:
            %       % Construct a composite gate and find its inverse
            %       gates = [hGate(1); sGate(2); rxGate(1, pi/3)];
            %       circ = quantumCircuit(gates, Name="myCircuit");
            %       cg = compositeGate(circ, [4 2]);
            %       invCG = inv(cg);
            %       cg.Gates
            %       invCG.Gates
            %       invCG.Name
            %
            %   See also quantumCircuit/inv, quantum.gate.CompositeGate

            arguments
                obj (1,1) quantum.gate.CompositeGate
            end

            lG = obj.Gates;
            lG = flip(lG);
            for ii=1:length(lG)
                lG(ii) = inv(lG(ii));
            end
            obj.Gates = lG;

            % Adjust the name of the CompositeGate by adding / removing "_inv".
            if strlength(obj.Name) ~= 0
                if endsWith(obj.Name, "_inv")
                    n = char(obj.Name);
                    n(end-3:end) = [];
                    obj.Name = string(n);
                else
                    obj.Name = obj.Name + "_inv";
                end
            end
        end

        function cp = plot(obj, varargin)
            %PLOT  Plots a CompositeGate
            %
            %   PLOT(cg) plots the CompositeGate cg.
            %
            %   PLOT(cg, QubitBlocks=blocks) distinguishes blocks of qubits
            %   by separating them with a red dashed line in the plot. The
            %   input blocks should be a vector, where each element is the
            %   size of a block. The block sizes must sum up to the length
            %   of cg.TargetQubits.
            %
            %   PLOT(___, NumRows=nr) specifies the number of rows to use
            %   when wrapping a CompositeGate over multiple rows. By
            %   default, NumRows is determined based on the input circuit.
            %
            %   PLOT(___, QubitLabelLocation=loc) specifies the location of labels for
            %   the qubit lines as one of the following values
            %        'left'    - Qubit lines are labeled on the left.
            %        'right'   - Qubit lines are labeled on the right.
            %        'none'    - Qubit lines are not labeled.
            %        'both'    - Qubit lines are labeled on the left and right.
            %   By default, this value is chosen based on the input
            %   composite gate.
            %
            %   PLOT(parent,cg,___) plots CompositeGate cg in the figure,
            %   panel, or tab specified by parent.
            %
            %   C = PLOT(...) returns a QuantumCircuitChart object.
            %   Use the methods and properties of this object to inspect and
            %   adjust the plotted composite gate.
            %
            %   Example:
            %       % Construct and plot CompositeGate.
            %       gates = [hGate(1); cxGate(1, 2)];
            %       circ = quantumCircuit(gates);
            %       cg = compositeGate(circ, [4 2]);
            %       plot(cg)
            %       % Plotting a CompositeGate shows the inner gates. 
            %       % Mapping to the outer circuit is only shown when
            %       % plotting a quantumCircuit containing the CompositeGate.
            %
            %   See also quantum.gate.CompositeGate, quantumCircuit/plot

            args = varargin;
            nameOffset = 1;

            % Check if the first input argument is a graphics object to use as parent.
            if isa(obj,'matlab.graphics.Graphics')
                % plot(parent,cg,___)
                args = [args(:),{'Parent',obj}];
                obj = args{1};
                args(1) = [];
                nameOffset = 2;
            end

            if ~isa(obj, 'quantum.gate.CompositeGate') || ~isscalar(obj)
                error(message('quantum:quantumCircuit:plotMustBeCompositeGate', nameOffset))
            end

            numQubits = length(obj.TargetQubits)+length(obj.ControlQubits);

            if obj.Name == ""
                title = "CompositeGate";
            else
                % Escape underscore since title is in 'tex' mode:
                title = "CompositeGate: " + replace(obj.Name, '_', '\_');
            end

            quantum.internal.gate.checkQubitBlocksForPlot(numQubits, args)

            p = quantum.gate.QuantumCircuitChart('NumQubits', numQubits, 'Gates', obj.Gates, ...
                'Title', title, args{:});

            if nargout > 0
                cp = p;
            end
        end
    end

    methods (Hidden)
        function qb = getQubits(obj)
            qb = obj.TargetQubits;
        end

        function obj = setQubits(obj, qb)
            mustBeTargetQubits(qb)
            assert(length(qb) == length(obj.TargetQubits_)); % Assert enough for this hidden method
            obj.TargetQubits_ = qb;
        end

        function inputState = applyToState(obj,inputState,numQubits)
            % Iteratively implements applytoState method of each gate

            if numel(inputState) < 2^length(obj.TargetQubits)
                error(message("quantum:CompositeGate:tooFewQubits", length(obj.TargetQubits)));
            end

            for i=1:numel(obj.Gates)
                g = obj.Gates(i);
                qb = getQubits(g);
                g = setQubits(g, obj.TargetQubits(qb));
                inputState = applyToState(g, inputState,numQubits);
            end
        end

        % Generate QASM code for this gate
        function [instr, def, map] = getQASM(obj, Unpack, map, varargin)

            qasm_name = lower(obj.Name);
            body = "";
            child_gate_defs = "";
            
            for i = 1:length(obj.Gates)
                g = obj.Gates(i);
                if Unpack
                    % Unpacked: Gates are inlined and no definition blocks are used.
                    qb = getQubits(g);
                    [child_instr, ~, ~] = getQASM(setQubits(g, obj.TargetQubits(qb)), Unpack, map, varargin{:});
                else
                    % Defined (default): Gate definition blocks are used with their instructions
                    [child_instr, child_def, map] = getQASM(g, Unpack, map, varargin{:});
                    child_gate_defs = child_gate_defs + child_def;
                end
                body = body + child_instr;
            end

            if Unpack
                instr = body;
                def = "";
                return
            else

                % NOTE: Depending on what qubits the children gates act on,
                % qargs could contain input qubits that aren't used within the
                % body...not a problem for now, but might be a limitiation
                % when/if vendors support custom gate definitions. Although, I
                % haven't seen anything in the OpenQASM spec preventing this,
                % so its not really a concern.

                % Format as definition
                body = strtrim(regexprep(body, ["[","]"], ""));
                qargs = join("q"+string((1:length(obj.TargetQubits))-1)+"", ",");

                if isKey(map, body)
                    % Gates are defined once
                    qasm_name = map(body);
                    def = "";
                else
                    defined_gates = string(values(map));
                    if strcmp(qasm_name, "")
                        qasm_name = quantum.internal.gate.makeUniqueName(defined_gates);
                    else
                        if any(strcmp(defined_gates, qasm_name))
                            qasm_name = quantum.internal.gate.makeUniqueName(defined_gates, qasm_name);
                        end
                    end
                    map(body) = qasm_name;
                    % NOTE: Does not support symbolic parameters in gate definitions
                    def = child_gate_defs + sprintf("gate %s %s {\n%s\n}\n", qasm_name, qargs, body);
                end
                instr = quantum.internal.gate.formatQASMInstruction(obj, qasm_name);
            end 
        end

        function [types, ctrls, trgts, angles] = getProperties(obj)
            [types, ctrls, trgts, angles] = quantum.internal.gate.getProperties(obj.Gates);
            ctrls = obj.TargetQubits(ctrls);
            trgts = obj.TargetQubits(trgts);
        end
    end

    properties(Constant, Access='private')
        % Version of the serialization and deserialization
        % format. This is used for managing forward compatibility. Value is
        % saved in 'versionSavedFrom' when an instance is serialized.
        %
        %   1.0 : original shipping version
        version = 1.0;
    end
    methods (Hidden)
        function s = saveobj(g)
            s = struct('Name', g.Name, ...
                'TargetQubits', g.TargetQubits, ...
                'Gates', g.Gates, ...
                'versionSavedFrom', quantum.gate.CompositeGate.version, ...
                'minCompatibleVersion', 1);
        end
    end
    methods(Hidden, Static)
        function g = loadobj(s)
            if quantum.gate.CompositeGate.version < s.minCompatibleVersion
                id = 'quantum:CompositeGate:IncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.CompositeGate', getString(message(id))));
                warning(id, loadWarningString);

                g = quantum.gate.CompositeGate(quantumCircuit(0), zeros(1, 0));
                return
            end

            try
                g = compositeGate(s.Gates, s.TargetQubits, Name=s.Name);
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.CompositeGate', err.message));
                warning(err.identifier, loadWarningString);

                g = quantum.gate.CompositeGate(quantumCircuit(0), zeros(1, 0));
                return
            end
        end
    end
end

function mustBeTargetQubits(qubits)
if ~isnumeric(qubits) || ~isreal(qubits) || ~allfinite(qubits) || ...
        ~(isvector(qubits) || (size(qubits, 1) == 0 && size(qubits, 2) == 0))
    error(message("quantum:CompositeGate:invalidTargetQubits"));
end
if ~(all(floor(qubits) == qubits) && all(qubits >= 1))
    error(message("quantum:CompositeGate:invalidTargetQubits"));
end

if any(diff(sort(qubits)) == 0)
    error(message("quantum:gates:matchingQubits"));
end
end
