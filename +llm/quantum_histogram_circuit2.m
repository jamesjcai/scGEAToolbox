function circuit = generate_quantum_histogram_circuit_matlab(histogram, num_qubits)
    % GENERATE_QUANTUM_HISTOGRAM_CIRCUIT_MATLAB - Generate a quantum circuit to produce a given histogram
    % using MATLAB's Quantum Computing Support Package
    %
    % Inputs:
    %   histogram - array of frequencies/probabilities for each bitstring
    %   num_qubits - number of qubits in the circuit
    %
    % Output:
    %   circuit - quantum circuit object from MATLAB's package
    %
    % Requires: MATLAB Quantum Computing Support Package
    
    % Check if the quantum computing package is available
    if ~exist('quantumCircuit', 'file')
        error('MATLAB Quantum Computing Support Package is required. Install from Add-On Explorer.');
    end
    
    % Ensure histogram is normalized (probabilities sum to 1)
    histogram = histogram / sum(histogram);
    
    % Check if histogram size matches 2^num_qubits
    if length(histogram) ~= 2^num_qubits
        error('Histogram size must be 2^num_qubits');
    end
    
    % Create quantum circuit
    circuit = quantumCircuit(num_qubits);
    
    % Method 1: Use MATLAB's built-in state preparation if available
    if exist('statePreparation', 'file')
        % Use MATLAB's optimized state preparation
        target_amplitudes = sqrt(histogram);
        circuit = statePreparation(circuit, target_amplitudes);
    else
        % Method 2: Manual implementation using MATLAB's quantum gates
        circuit = implement_state_preparation(circuit, histogram, num_qubits);
    end
    
    % Display circuit information
    fprintf('Generated quantum circuit with %d qubits\n', num_qubits);
    fprintf('Circuit depth: %d\n', circuitDepth(circuit));
    fprintf('Total gates: %d\n', numGates(circuit));
    
    % Verify the circuit by simulation
    if nargout > 0 || true  % Always verify
        verify_circuit(circuit, histogram);
    end
end

function circuit = implement_state_preparation(circuit, histogram, num_qubits)
    % Manual implementation of state preparation using MATLAB's quantum gates
    
    target_amplitudes = sqrt(histogram);
    
    % Recursive state preparation approach
    circuit = recursive_state_preparation(circuit, target_amplitudes, 1:num_qubits, []);
end

function circuit = recursive_state_preparation(circuit, amplitudes, qubits, controls)
    % Recursive function to prepare quantum state
    
    if length(qubits) == 1
        % Base case: single qubit
        if length(amplitudes) == 2
            % Calculate rotation angle
            if abs(amplitudes(1)) < 1e-10
                theta = pi/2;
            elseif abs(amplitudes(2)) < 1e-10
                theta = 0;
            else
                theta = 2 * atan(abs(amplitudes(2)) / abs(amplitudes(1)));
            end
            
            % Apply rotation
            if abs(theta) > 1e-10
                if isempty(controls)
                    circuit = addGate(circuit, ryGate(theta), qubits(1));
                else
                    circuit = addGate(circuit, ryGate(theta), qubits(1), controls);
                end
            end
        end
        return;
    end
    
    % Recursive case: multiple qubits
    n = length(qubits);
    half_size = 2^(n-1);
    
    % Split amplitudes into two halves
    upper_half = amplitudes(1:half_size);
    lower_half = amplitudes(half_size+1:end);
    
    % Calculate the probability of being in upper vs lower half
    p_upper = sum(abs(upper_half).^2);
    p_lower = sum(abs(lower_half).^2);
    
    % Normalize the conditional states
    if p_upper > 1e-10
        upper_half = upper_half / sqrt(p_upper);
    end
    if p_lower > 1e-10
        lower_half = lower_half / sqrt(p_lower);
    end
    
    % Calculate rotation angle for the most significant qubit
    if p_upper < 1e-10
        theta = pi/2;
    elseif p_lower < 1e-10
        theta = 0;
    else
        theta = 2 * atan(sqrt(p_lower) / sqrt(p_upper));
    end
    
    % Apply rotation to the most significant qubit
    if abs(theta) > 1e-10
        if isempty(controls)
            circuit = addGate(circuit, ryGate(1, theta), qubits(1));
        else
            circuit = addGate(circuit, ryGate(1, theta), qubits(1), controls);
        end
    end
    
    % Recursively prepare the conditional states
    remaining_qubits = qubits(2:end);
    
    % Prepare state conditioned on qubit being 0
    if p_upper > 1e-10
        circuit = recursive_state_preparation(circuit, upper_half, remaining_qubits, [controls, qubits(1)]);
    end
    
    % Prepare state conditioned on qubit being 1
    if p_lower > 1e-10
        % Add X gate to flip the control condition
        circuit = addGate(circuit, xGate(), qubits(1), controls);
        circuit = recursive_state_preparation(circuit, lower_half, remaining_qubits, [controls, qubits(1)]);
        % Flip back
        circuit = addGate(circuit, xGate(), qubits(1), controls);
    end
end

function verify_circuit(circuit, target_histogram)
    % Verify the circuit by simulating it
    try
        % Create initial state (all qubits in |0>)
        num_qubits = circuit.NumQubits;
        
        % Simulate the circuit
        if exist('simulate', 'file')
            result = simulate(circuit);
            probabilities = abs(result.StateVector).^2;
        else
            % Alternative simulation method
            sim = quantumSimulator('Method', 'statevector');
            result = run(sim, circuit);
            probabilities = abs(result.StateVector).^2;
        end
        
        % Compare with target
        fprintf('Verification Results:\n');
        fprintf('Target vs Achieved Probabilities:\n');
        for i = 1:length(target_histogram)
            fprintf('State %d: %.4f vs %.4f (error: %.4f)\n', i-1, target_histogram(i), probabilities(i), abs(target_histogram(i) - probabilities(i)));
        end
        
        % Calculate fidelity
        fidelity = sum(sqrt(target_histogram .* probabilities));
        fprintf('Fidelity: %.6f\n', fidelity);
        
        % Plot comparison
        plot_histogram_comparison_matlab(target_histogram, probabilities);
        
    catch ME
        fprintf('Circuit verification failed: %s\n', ME.message);
    end
end

function plot_histogram_comparison_matlab(target, achieved)
    % Plot comparison between target and achieved histograms
    figure;
    
    x = 0:length(target)-1;
    
    % Create subplot for target histogram
    subplot(2,1,1);
    bar(x, target, 'FaceColor', [0.3, 0.6, 0.8]);
    title('Target Histogram');
    xlabel('Bitstring (decimal)');
    ylabel('Probability');
    grid on;
    
    % Create subplot for achieved histogram
    subplot(2,1,2);
    bar(x, achieved, 'FaceColor', [0.8, 0.3, 0.3]);
    title('Achieved Histogram');
    xlabel('Bitstring (decimal)');
    ylabel('Probability');
    grid on;
    
    % Add error subplot
    figure;
    error_values = abs(target - achieved);
    bar(x, error_values, 'FaceColor', [0.8, 0.8, 0.3]);
    title('Absolute Error');
    xlabel('Bitstring (decimal)');
    ylabel('|Target - Achieved|');
    grid on;
end

function display_circuit(circuit)
    % Display the quantum circuit
    if exist('drawCircuit', 'file')
        figure;
        drawCircuit(circuit);
        title('Quantum Circuit for Histogram Generation');
    else
        % Alternative display method
        disp(circuit);
    end
end

%% Example usage with MATLAB Quantum Computing Support Package:

% Example 1: Simple 2-qubit histogram
fprintf('=== Example 1: 2-qubit system ===\n');
histogram_2q = [0.1, 0.3, 0.4, 0.2];
circuit_2q = generate_quantum_histogram_circuit_matlab(histogram_2q, 2);

% Example 2: 3-qubit histogram with more complex distribution
fprintf('\n=== Example 2: 3-qubit system ===\n');
histogram_3q = [0.05, 0.15, 0.1, 0.2, 0.25, 0.1, 0.1, 0.05];
circuit_3q = generate_quantum_histogram_circuit_matlab(histogram_3q, 3);

% Example 3: Uniform distribution (for comparison)
fprintf('\n=== Example 3: Uniform distribution ===\n');
histogram_uniform = ones(1, 8) / 8;
circuit_uniform = generate_quantum_histogram_circuit_matlab(histogram_uniform, 3);

% Display circuit information
fprintf('\n=== Circuit Analysis ===\n');
if exist('circuitDepth', 'file')
    fprintf('Circuit depth: %d\n', circuitDepth(circuit_3q));
end
if exist('numGates', 'file')
    fprintf('Number of gates: %d\n', numGates(circuit_3q));
end

% Optional: Display the circuit diagram
% display_circuit(circuit_3q);