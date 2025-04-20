# Updated imports for current Qiskit structure
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
import numpy as np

def oracle_for_state(circuit, qubits, target_state):
    """
    Parameterized oracle that marks a specific target state.
    
    Args:
        circuit: The quantum circuit
        qubits: List of qubit indices
        target_state: Binary string representing the target state (e.g. '011')
    """
    # Apply X gates to qubits that should be 0 in the target state
    for i, bit in enumerate(target_state):
        if bit == '0':
            circuit.x(qubits[i])
    
    # Apply multi-controlled phase flip (Z) to the last qubit
    if len(qubits) == 3:
        # For 3 qubits, we can use the CCZ operation (or an equivalent)
        # We'll implement it using CCX (Toffoli) with the help of an h gate
        circuit.h(qubits[-1])
        circuit.ccx(qubits[0], qubits[1], qubits[2])
        circuit.h(qubits[-1])
    elif len(qubits) > 3:
        # Implement multi-controlled Z for larger systems
        circuit.h(qubits[-1])
        circuit.mcx(qubits[:-1], qubits[-1])
        circuit.h(qubits[-1])
    
    # Undo the X gates
    for i, bit in enumerate(target_state):
        if bit == '0':
            circuit.x(qubits[i])
    
    return circuit

def diffuser(circuit, qubits):
    """
    Diffusion operator that implements the reflection about the average.
    Also known as the Grover diffusion operator.
    """
    # Transform to Hadamard basis
    circuit.h(qubits)
    
    # Apply phase flip to |0⟩ (mark |00...0> state)
    circuit.x(qubits)  # Convert |0⟩ to |1⟩
    
    # Apply multi-controlled Z to mark |11...1> state
    circuit.h(qubits[-1])
    if len(qubits) == 3:
        circuit.ccx(qubits[0], qubits[1], qubits[2])
    else:
        circuit.mcx(qubits[:-1], qubits[-1])
    circuit.h(qubits[-1])
    
    # Undo the X gates
    circuit.x(qubits)
    
    # Transform back to computational basis
    circuit.h(qubits)
    
    return circuit

def grover_search(n_qubits=3, target_state='011', iterations=1):
    """
    Implements Grover's algorithm for a specified number of qubits and iterations.
    
    Args:
        n_qubits: Number of qubits in the circuit
        target_state: The binary string representing the marked state
        iterations: Number of Grover iterations
    
    Returns:
        QuantumCircuit: Compiled circuit
    """
    # Validate target_state
    if len(target_state) != n_qubits:
        raise ValueError(f"Target state length ({len(target_state)}) must match number of qubits ({n_qubits})")
    
    # Initialize circuit with qubits and classical bits for measurement
    qc = QuantumCircuit(n_qubits, n_qubits)
    
    # Create equal superposition
    qc.h(range(n_qubits))
    
    # Apply Grover iterations
    for _ in range(iterations):
        # Apply oracle with the specified target state
        oracle_for_state(qc, range(n_qubits), target_state)
        
        # Apply diffusion operator
        diffuser(qc, range(n_qubits))
    
    # Measure all qubits
    qc.measure(range(n_qubits), range(n_qubits))
    
    return qc

def run_grover(circuit, shots=4096):
    """
    Executes the Grover circuit and returns results.
    
    Args:
        circuit: The quantum circuit to execute
        shots: Number of repetitions of the circuit
        
    Returns:
        dict: Measurement counts
    """
    try:
        simulator = AerSimulator()
        result = simulator.run(circuit, shots=shots).result()
        counts = result.get_counts()
        return counts
    except Exception as e:
        print(f"Simulation failed with error: {e}")
        try:
            # Try an alternative approach if available
            from qiskit import Aer
            backend = Aer.get_backend('qasm_simulator')
            job = backend.run(circuit, shots=shots)
            result = job.result()
            counts = result.get_counts()
            return counts
        except Exception as e2:
            print(f"Alternative simulation approach failed with error: {e2}")
            return {"Error": "Could not run simulation with any available method"}

# Function to verify correct behavior
def verify_grover(target_state='011', n_qubits=3, iterations=1, shots=4096):
    """Run Grover's algorithm and verify it finds the correct target state."""
    print(f"Running Grover's algorithm to find state: |{target_state}⟩")
    
    # Create and run circuit
    grover_circuit = grover_search(n_qubits=n_qubits, target_state=target_state, iterations=iterations)
    counts = run_grover(grover_circuit, shots=shots)
    
    # Display the results
    print("Measurement counts:", counts)
    print(f"Expected state: {target_state} (decimal {int(target_state, 2)})")
    
    # Calculate success metrics
    if "Error" not in counts:
        total_shots = sum(counts.values())
        target_count = counts.get(target_state, 0)
        target_probability = target_count / total_shots
        
        # Calculate the theoretical probability
        N = 2**n_qubits  # number of states
        theta = np.arcsin(np.sqrt(1/N))
        theoretical_prob = np.sin((2*iterations + 1) * theta)**2
        
        print(f"Measured probability for target state: {target_probability:.4f}")
        print(f"Theoretical probability for target state: {theoretical_prob:.4f}")
        
        # Plot histogram
        plot_histogram(counts)
        plt.title(f"Grover's Algorithm Results (Target: {target_state})")
        plt.show()
        
        # Check if the target state has the highest probability
        most_probable_state = max(counts, key=counts.get)
        if most_probable_state == target_state:
            print("Success! Target state has the highest probability.")
        else:
            print(f"Warning: State |{most_probable_state}⟩ has higher probability than target state |{target_state}⟩")
            print("This suggests an issue with the oracle or diffuser implementation.")
            
            # For diagnostic purposes, print bit-flipped state
            bit_flipped = ''.join('1' if bit == '0' else '0' for bit in target_state)
            if most_probable_state == bit_flipped:
                print(f"Note: The highest probability state |{most_probable_state}⟩ is the bit-flip of target |{target_state}⟩")
                print("This suggests the oracle may be marking the wrong state.")
    else:
        print("Simulation failed.")

# Test with target state '011'
verify_grover(target_state='011', iterations=1)