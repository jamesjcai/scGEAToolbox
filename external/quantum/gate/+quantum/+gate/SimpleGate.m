classdef (Sealed) SimpleGate < quantum.internal.gate.gates.Gate & quantum.gate.QuantumGate
    %SIMPLEGATE Simple gate for quantum computing
    %
    %   Use gate creation functions (e.g., hGate, cxGate, swapGate) to
    %   construct a SimpleGate object. These objects can be found in the
    %   Gates property of a quantumCircuit object.
    %
    %   Example:
    %       % Construct a SimpleGate of type "h" on qubit 2:
    %       hGate(2)
    %
    %       % Construct a SimpleGate of type "cx":
    %       cxGate(2, 3)
    %
    %       % Construct a SimpleGate of type "rx":
    %       rxGate(2, pi/3)
    %
    %       % Construct a SimpleGate of type "swap":
    %       swapGate(2, 3)
    %
    %   SimpleGate properties:
    %      Type             - Type of the gate.
    %      ControlQubits    - Qubits the gate is controlled by.
    %      TargetQubits     - Qubits the gate acts on.
    %      Angles           - Angles configuring the gate.
    %
    %   SimpleGate methods:
    %      inv              - Return inverse of a SimpleGate object.
    %      getMatrix        - Unitary matrix representation of gate.
    %
    %   SimpleGate creation functions:
    %      hGate       Hadamard gate
    %      idGate      Identity gate
    %      xGate       Pauli X gate
    %      yGate       Pauli Y gate
    %      zGate       Pauli Z gate
    %      sGate       S gate
    %      siGate      Inverse S gate
    %      tGate       T gate
    %      tiGate      Inverse T gate
    %      chGate      Controlled Hadamard gate
    %      cxGate      Controlled Pauli X gate
    %      cyGate      Controlled Pauli Y gate
    %      czGate      Controlled Pauli Z gate
    %      cnotGate    CNOT or Controlled Pauli X gate
    %      swapGate    Swap gate
    %      rxGate      X-axis rotation gate
    %      ryGate      Y-axis rotation gate
    %      rzGate      Z-axis rotation gate
    %      r1Gate      Z-axis rotation gate with global phase
    %      crxGate     Controlled X-axis rotation gate
    %      cryGate     Controlled Y-axis rotation gate
    %      crzGate     Controlled Z-axis rotation gate
    %      cr1Gate     Controlled Z-axis rotation gate with global phase
    %      ccxGate     Controlled-controlled Pauli X gate
    %      rxxGate     Ising XX coupling gate
    %      ryyGate     Ising YY coupling gate
    %      rzzGate     Ising ZZ coupling gate
    %
    %   See also quantumCircuit, quantum.gate.CompositeGate

    %   Copyright 2021-2023 The MathWorks, Inc.

    properties(Dependent, SetAccess=private)
        %TYPE - Type of the gate
        %   Type is a string that represents the type of the gate. For
        %   example, Type is "cx" for a controlled X gate.
        %
        %   See also quantum.gate.SimpleGate, hGate, cxGate
        Type

        %CONTROLQUBITS - Control qubits of the gate
        %   ControlQubits is a row vector of qubit indices, indicating
        %   which qubits the gate is controlled by. The length of
        %   ControlQubits depends on the Type of the gate - for some gates,
        %   ControlQubits is empty.
        %
        %   See also quantum.gate.SimpleGate, cryGate, ccxGate
        ControlQubits

        %TARGETQUBITS - Qubits the gate acts on
        %   TargetQubits is a row vector of qubit indices, indicating
        %   which qubits the gate acts on. The length of
        %   TargetQubits depends on the Type of the gate.
        %
        %   See also quantum.gate.SimpleGate, hGate, swapGate
        TargetQubits
    end


    properties(Dependent)
        %ANGLES - Angles configuring the gate
        %   Angles is a row vector of angles by which a gate rotates or
        %   applies phase shifts to its target qubit(s). The length of
        %   Angles depends on the Type of the gate - for many gates, Angles
        %   is empty.
        %
        %   See also quantum.gate.SimpleGate, rxGate, rxxGate
        Angles
    end

    properties (Access = private)
        InnerGate; % Contains the numerics / code generation parts
    end

    properties(Dependent, Hidden, SetAccess = private, GetAccess = public)
        nameForPlot;
    end

    methods(Access=protected)
        function str = getOneLineDisplayString(obj)
            if ~isempty(obj.Angles)
                str = {obj.Type, obj.ControlQubits, obj.TargetQubits, ...
                    quantum.internal.gate.makeAngleString(obj.Angles, "pi"), []};
            else
                str = {obj.Type, obj.ControlQubits, obj.TargetQubits, "", []};
            end
        end
    end

    methods
        function n = get.nameForPlot(obj)
            n = obj.InnerGate.nameForPlot;
        end

        function q = get.TargetQubits(obj)
            q = getTargetQubits(obj.InnerGate);
        end

        function q = get.ControlQubits(obj)
            q = getControlQubits(obj.InnerGate);
        end

        function a = get.Angles(obj)
            a = angles(obj.InnerGate);
        end

        function t = get.Type(obj)
             t = obj.InnerGate.type;
        end

        function obj = set.Angles(obj, ang)
            currentAngles = angles(obj.InnerGate);
            if ~isnumeric(ang) || ~isreal(ang) || ~isvector(ang) || ...
                    length(ang) ~= length(currentAngles)
                error(message('quantum:gates:SimpleGateInvalidAngles'))
            end
            ang = reshape(double(full(ang)), 1, []);
            obj.InnerGate = setAngles(obj.InnerGate, ang);
        end
    end

    methods
        % constructor
        function obj = SimpleGate(type, varargin)
            arguments
                type {mustBeTextScalar}
            end
            arguments(Repeating)
                varargin
            end
            [nrControl, nrTarget, nrAngles, gateHandle] = processType(type);

            narginchk(nrControl+nrTarget+nrAngles+1, nrControl+nrTarget+nrAngles+1)
            for ii=1:nrControl+nrTarget
                checkQubit(varargin{ii});
            end
            for ii=nrControl+nrTarget+1:nrControl+nrTarget+nrAngles
                checkAngles(varargin{ii});
            end

            n = 1;
            for ii=1:length(varargin)
                if ~isscalar(varargin{ii})
                    % Already checked that this is a vector in checkQubit
                    if n == 1
                        n = length(varargin{ii});
                    elseif length(varargin{ii}) ~= n
                        error(message("quantum:gates:dimMismatch"));
                    end
                end
            end

            if n==0
                obj = quantum.gate.SimpleGate.empty(0,1);
            elseif n==1
                obj.InnerGate = gateHandle(varargin{:});
            elseif n>1
                for ii=1:length(varargin)
                    if isscalar(varargin{ii})
                        varargin{ii} = repmat(varargin{ii}, n, 1);
                    end
                end

                c = cellfun(@(x) x(1), varargin, 'UniformOutput', false);
                obj.InnerGate = gateHandle(c{:});
                obj = repmat(obj, n, 1);
                for jj=2:n

                    c = cell(length(varargin),1);
                    for ii = 1:length(varargin)
                        c{ii} = varargin{ii}(jj);
                    end

                    g = gateHandle(c{:});
                    obj(jj).InnerGate = g;
                end
            end
        end

        function obj = inv(obj)
            %INV  Return inverse of SimpleGate
            %
            %   invg = INV(g) returns a new SimpleGate invg with inverse behavior
            %   of g. This is also a SimpleGate, but may have a different type or
            %   different angles.
            %
            %   Example:
            %       % Compute inverse of some SimpleGate objects
            %
            %       % The inverse of rxGate rotates in the opposite direction
            %       inv(rxGate(3, pi/3))
            %
            %       % The sGate has an associated siGate which is its inverse
            %       inv(sGate(1))
            %
            %       % The hGate is its own inverse
            %       inv(hGate(2))
            %
            %   See also quantum.gate.SimpleGate, quantumCircuit/inv

            arguments
                obj (1,1) quantum.gate.SimpleGate
            end
            obj.InnerGate = inv(obj.InnerGate);
        end
    end

    methods (Hidden) %(Access = ?quantumCircuit)
        function inputState = applyToState(obj,inputState,numQubits)
            inputState = applyToState(obj.InnerGate, inputState, numQubits);
        end

        % Generate QASM code for this gate
        function [instr, def, map] = getQASM(obj, varargin)
            [instr, def, map] = getQASM(obj.InnerGate, varargin{:});
        end

        function obj = setQubits(obj, qubits)
            obj.InnerGate = setQubits(obj.InnerGate, qubits);
        end

        function qb = getQubits(obj)
            qb = getQubits(obj.InnerGate);
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
            s = struct('Type', g.Type, ...
                'ControlQubits', g.ControlQubits, ...
                'TargetQubits', g.TargetQubits, ...
                'Angles', g.Angles, ...
                'versionSavedFrom', quantum.gate.SimpleGate.version, ...
                'minCompatibleVersion', 1);
        end
    end
    methods(Hidden, Static)
        function g = loadobj(s)
            if quantum.gate.SimpleGate.version < s.minCompatibleVersion
                id = 'quantum:gates:SimpleGateIncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.SimpleGate', getString(message(id))));
                warning(id, loadWarningString);

                g = quantum.gate.SimpleGate.MissingGate();
                return
            end

            try
                args = [num2cell(s.ControlQubits), num2cell(s.TargetQubits), num2cell(s.Angles)];

                g = quantum.gate.SimpleGate(s.Type, args{:});
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.SimpleGate', err.message));
                warning(err.identifier, loadWarningString);

                g = quantum.gate.SimpleGate.MissingGate();
                return
            end
        end

        function g = MissingGate()
            g = quantum.gate.SimpleGate('H', 1); % arbitrary SimpleGate to modify
            g.InnerGate = quantum.internal.gate.gates.Missing;
        end
    end
end

function checkQubit(qubits)
if ~isnumeric(qubits) || ~isreal(qubits) || ~allfinite(qubits) || ...
    ~(isvector(qubits) || (size(qubits, 1) == 0 && size(qubits, 2) == 0))
    error(message("quantum:gates:invalidQubitVector"));
end
if ~(all(floor(qubits) == qubits) && all(qubits >= 1))
    error(message("quantum:gates:invalidQubitVector"));
end
end

function checkAngles(qubits)
if ~isnumeric(qubits) || ~isreal(qubits) || ~allfinite(qubits) || ...
    ~(isvector(qubits) || (size(qubits, 1) == 0 && size(qubits, 2) == 0))
    error(message("quantum:gates:invalidAngle"));
end
end

function [Nctrl, Ntrgt, Nangle, constr] = processType(type)
type = lower(char(type));
switch type
    case 'ccx'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CCX.numInputs;
        constr = @quantum.internal.gate.gates.CCX;
    case 'ch'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CH.numInputs;
        constr = @quantum.internal.gate.gates.CH;
    case 'crx'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CRX.numInputs;
        constr = @quantum.internal.gate.gates.CRX;
    case 'cry'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CRY.numInputs;
        constr = @quantum.internal.gate.gates.CRY;
    case 'crz'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CRZ.numInputs;
        constr = @quantum.internal.gate.gates.CRZ;
    case 'cr1'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CR1.numInputs;
        constr = @quantum.internal.gate.gates.CR1;
    case 'cx'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CX.numInputs;
        constr = @quantum.internal.gate.gates.CX;
    case 'cy'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CY.numInputs;
        constr = @quantum.internal.gate.gates.CY;
    case 'cz'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.CZ.numInputs;
        constr = @quantum.internal.gate.gates.CZ;
    case 'h'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.H.numInputs;
        constr = @quantum.internal.gate.gates.H;
    case 'id'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.ID.numInputs;
        constr = @quantum.internal.gate.gates.ID;
    case 'x'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.X.numInputs;
        constr = @quantum.internal.gate.gates.X;
    case 'y'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.Y.numInputs;
        constr = @quantum.internal.gate.gates.Y;
    case 'z'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.Z.numInputs;
        constr = @quantum.internal.gate.gates.Z;
    case 's'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.S.numInputs;
        constr = @quantum.internal.gate.gates.S;
    case 'si'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.SI.numInputs;
        constr = @quantum.internal.gate.gates.SI;
    case 't'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.T.numInputs;
        constr = @quantum.internal.gate.gates.T;
    case 'ti'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.TI.numInputs;
        constr = @quantum.internal.gate.gates.TI;
    case 'rx'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RX.numInputs;
        constr = @quantum.internal.gate.gates.RX;
    case 'ry'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RY.numInputs;
        constr = @quantum.internal.gate.gates.RY;
    case 'rz'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RZ.numInputs;
        constr = @quantum.internal.gate.gates.RZ;
    case 'rxx'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RXX.numInputs;
        constr = @quantum.internal.gate.gates.RXX;
    case 'ryy'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RYY.numInputs;
        constr = @quantum.internal.gate.gates.RYY;
    case 'rzz'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.RZZ.numInputs;
        constr = @quantum.internal.gate.gates.RZZ;
    case 'r1'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.R1.numInputs;
        constr = @quantum.internal.gate.gates.R1;
    case 'swap'
        [Nctrl, Ntrgt, Nangle] = quantum.internal.gate.gates.SWAP.numInputs;
        constr = @quantum.internal.gate.gates.SWAP;
    otherwise
        error(message("quantum:gates:GateNotSupported"))
end
end
