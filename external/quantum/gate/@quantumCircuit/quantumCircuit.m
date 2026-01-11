classdef (Sealed, ...
        InferiorClasses = {?matlab.ui.Figure, ?matlab.ui.container.Tab, ?matlab.ui.container.Panel, ...
        ?matlab.graphics.axis.Axes, ?matlab.ui.control.UIAxes, ...
        ?matlab.ui.container.GridLayout, ?matlab.graphics.layout.TiledChartLayout}) ...
        quantumCircuit
    %QUANTUMCIRCUIT Quantum computing circuit
    %   C = quantumCircuit(numQubits) builds a quantum computing circuit,
    %   which is a quantumCircuit object, with numQubits qubits and no
    %   gates.
    %
    %   C = quantumCircuit(gates) takes an array of quantum gates and
    %   builds a quantumCircuit with those gates and as many qubits as used
    %   in those gates.
    %
    %   C = quantumCircuit(gates, numQubits) takes an array of quantum
    %   gates and builds a quantumCircuit with those gates and numQubits
    %   qubits. This syntax errors if any gate uses a qubit with index
    %   larger than numQubits.
    %
    %   C = quantumCircuit(___, Name=name) additionally specifies a name
    %   for the circuit. The default name is "".
    %
    %   Example:
    %       % Construct a quantumCircuit, plot the circuit
    %       % and simulate running the circuit:
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       plot(circ)
    %       s = simulate(circ)
    %       histogram(s)
    %
    %   quantumCircuit properties:
    %      NumQubits        - Number of qubits in the circuit.
    %      Gates            - Array of gates in the circuit.
    %      Name             - Name of the circuit.
    %
    %   quantumCircuit methods:
    %      plot             - Plot a quantumCircuit.
    %      simulate         - Simulate a quantumCircuit.
    %      run              - Run a quantumCircuit on a quantum device.
    %
    %      getMatrix        - Unitary matrix representation of a quantumCircuit.
    %      inv              - Inverse of a quantumCircuit.
    %
    %      generateQASM     - Generate OpenQASM code for a quantumCircuit.
    %      unpack           - Unpack CompositeGate objects in a quantumCircuit.
    %
    %   See also quantum.gate.SimpleGate, quantum.gate.CompositeGate, quantum.gate.QuantumState

    %   Copyright 2021-2024 The MathWorks, Inc.

    properties(Dependent)
        %NUMQUBITS - Number of qubits in the circuit
        %   NUMQUBITS is the number of qubits in the circuit.
        %
        %   See also QUANTUMCIRCUIT
        NumQubits

        %GATES - Array of gates in the circuit
        %   GATES is a vector containing all the gates in the circuit. The
        %   elements of this vector are of type SimpleGate or CompositeGate.
        %
        %   See also quantumCircuit, quantum.gate.SimpleGate, quantum.gate.CompositeGate
        Gates
    end

    properties
        %NAME - Name of the circuit
        %   NAME is a string giving the name of the circuit. Default is "".
        Name {mustBeTextScalar} = ""
    end

    properties(Access = private)
        % properties for calculating circuit
        NumQubits_ = 0; % integer number of qubits in circuit
        Gates_ = quantum.gate.QuantumGate.empty(0, 1); % cell array of gate objects for user-defined circuit
    end

    methods
        function nq = get.NumQubits(obj)
            nq = obj.NumQubits_;
        end

        function obj = set.NumQubits(obj, numQubits)
            checkNumQubits(numQubits);
            quantum.internal.gate.checkGates(obj.Gates_, numQubits);
            obj.NumQubits_ = numQubits;
        end

        function gates = get.Gates(obj)
            gates = obj.Gates_;
        end

        function obj = set.Gates(obj, gateList)
            quantum.internal.gate.checkGates(gateList, obj.NumQubits_);
            obj.Gates_ = reshape(gateList, [], 1);
        end

        function obj = set.Name(obj, new_name)
            new_name = string(new_name);
            quantum.internal.gate.checkName(new_name);
            obj.Name = new_name;
        end

    end

    methods
        % constructor
        function obj = quantumCircuit(gates, numQubits, varargin)

            firstInputNumQubits = isnumeric(gates);

            if firstInputNumQubits
                checkNumQubits(gates);
                obj.NumQubits_ = gates;
                % Gates_ default value is empty, no need to set
            elseif isa(gates, 'quantum.gate.QuantumGate')
                if ~isvector(gates)
                    error(message('quantum:quantumCircuit:gatesMustBeVector'));
                end
                obj.Gates_ = reshape(gates, [], 1);
                % NumQubits_ will be set below
            else
                error(message('quantum:quantumCircuit:invalidFirstInput'));
            end

            if ~firstInputNumQubits && nargin > 1 && isnumeric(numQubits)
                % quantumCircuit(gates, numQubits, ['Name', n])
                checkNumQubits(numQubits);
                gatesMaxQubit = quantum.internal.gate.getMaxQubit(gates);
                if gatesMaxQubit > numQubits
                    error(message('quantum:quantumCircuit:invalidQubitInGate', gatesMaxQubit, numQubits));
                end
                obj.NumQubits_ = numQubits;
            else
                % quantumCircuit(gates, ['Name', n]), quantumCircuit(numQubits, ['Name', n])
                if ~firstInputNumQubits
                    obj.NumQubits_ = quantum.internal.gate.getMaxQubit(gates);
                end
                if nargin > 1
                    % Second input is option, not numQubits
                    varargin = [{numQubits} varargin];
                end
            end

            obj.Name = getNVP(varargin{:});
        end
    end

    properties(Constant, Access=private)
        % Version of the serialization and deserialization
        % format. This is used for managing forward compatibility. Value is
        % saved in 'versionSavedFrom' when an instance is serialized.
        %
        %   1.0 : original shipping version
        %   2.0 : support array
        version = 2.0;
    end
    methods (Hidden)
        function [types, ctrls, trgts, angles] = getProperties(obj)
            [types, ctrls, trgts, angles] = quantum.internal.gate.getProperties(obj.Gates);
        end

        function s = saveobj(c)
            % This is valid for the array case but only sees a scalar instance.
            s = struct('NumQubits', c.NumQubits, ...
                'Gates', c.Gates, ...
                'Name', c.Name, ...
                'versionSavedFrom', quantumCircuit.version, ...
                'minCompatibleVersion', 1);
        end
    end
    methods(Hidden, Static)
        function c = loadobj(s)
            % This is valid for the array case but only sees a scalar instance.
            if quantumCircuit.version < s.minCompatibleVersion
                id = 'quantum:quantumCircuit:IncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantumCircuit', getString(message(id))));
                warning(id, loadWarningString);

                c = quantumCircuit(0);
                return
            end

            try
                c = quantumCircuit(s.Gates, s.NumQubits, Name=s.Name);
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantumCircuit', err.message));
                warning(err.identifier, loadWarningString);

                c = quantumCircuit(0);
                return
            end
        end
    end
end

function name = getNVP(NameValueArgs)
% Helper to use argument block for NVP
arguments
    NameValueArgs.Name {mustBeTextScalar} = ""
end
name = NameValueArgs.Name;
end

function checkNumQubits(numQubits)
if ~(isnumeric(numQubits) && isscalar(numQubits) && isreal(numQubits) && ...
        floor(numQubits) == numQubits && isfinite(numQubits) && numQubits >= 0)
    error(message("quantum:quantumCircuit:invalidNumQubits"))
end
end
