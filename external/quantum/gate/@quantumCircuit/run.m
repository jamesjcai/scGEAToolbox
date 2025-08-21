function task = run(circuit, device, options)
%RUN  Run a quantumCircuit on a quantum device
%
%   task = RUN(circ,dev) runs the quantumCircuit circ remotely on a quantum
%   device dev. Input dev must be a QuantumDeviceAWS or QuantumDeviceIBM
%   object. The result is a QuantumTaskAWS or QuantumTaskIBM object, which
%   can be used to monitor the task and retrieve its result.
%
%   task = RUN(circ, dev, Name=Value) specifies additional options for
%   running the quantumCircuit circ on the quantum device dev. The available
%   options depend on the type of device:
%
%   QuantumDeviceAWS:
%
%   task = RUN(..., NumShots=n) runs the circuit n times remotely on
%   the quantum device. The default number is NumShots = 100.
%
%   QuantumDeviceIBM:
%
%   task = RUN(..., NumShots=n) runs the circuit n times remotely on
%   the quantum device. The default number is NumShots = 100.
%
%   task = RUN(..., UseErrorMitigation=tf) toggles error mitigation on the
%   measurement result. The default value of tf is true.
%
%   Example:
%       % Construct and run a quantumCircuit on a device using AWS.
%       % NOTE: This example requires you to have
%       % access to AWS Braket. When you start running the circuit on a
%       % device, you will be charged to your AWS account.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       dev = quantum.backend.QuantumDeviceAWS("Lucy");
%       task = run(circ, dev)
%       task.Status
%       wait(task)
%       measurement = fetchOutput(task)
%
%   Example:
%       % Construct and run a quantumCircuit on a device using IBM.
%       % NOTE: This example requires you to have
%       % access to Qiskit Runtime. When you start running the circuit on a
%       % device, you will be charged to your IBM account.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       dev = quantum.backend.QuantumDeviceIBM("ibmq_belem");
%       task = run(circ, dev)
%       task.Status
%       wait(task)
%       measurement = fetchOutput(task)
%
%   See also quantumCircuit, quantumCircuit/simulate,
%   quantum.backend.QuantumDeviceAWS, quantum.backend.QuantumTaskAWS,
%   quantum.backend.QuantumDeviceIBM, quantum.backend.QuantumTaskIBM

%   Copyright 2022-2024 The MathWorks, Inc.

arguments
    circuit quantumCircuit
    device {mustBeA(device, 'quantum.backend.QuantumDevice')}
end
arguments(Repeating)
    options
end

task = sendCircuit(device, circuit, options{:});
end

