from qiskit import qiskit
from math import pi

def unitary_operator_c(theta : float):
    """
    Returns a controlled Phase-Gate for the QFT

    Parameters: 
    theta: float
        is the value(theta) the the Phase-Gate applys 
    """
    return (qiskit.circuit.library.standard_gates.PhaseGate(theta,str(theta))).control(1)

def apply_unitary_operator(quantum_circuit : qiskit.circuit, theta : float, control : int, target : int, exponent : int):
    """
    Applies a controlled Phase-Gate applied on a circuit.
    The phasefactor is theta / 2**exponent

    Used for the QFT for each allies Phase-Gate

    Parameters:
    quantum_circuit: qiskit.circuit
        the circuit on which the gate will be applied
    theta: float
        is the value(theta) the the Phase-Gate applys 
    control: int
        position of the controll Qubit in the circuit
    target: int 
        position of the target of the Phase-Gate 
    exponent: int
        exponent to the power of two
    """
    quantum_circuit.append(unitary_operator_c(theta/(2**exponent)),[control, target])

def QFT_Gate(amount_qubits : int,m : int, inverse = False, swaps = False, MSB_first = True):
    """
    Returns the (approximative) Quantum-Fourier-Transformation as a gate. 

    Parameters:
    amount_qubits: int
        how many qubits are used for the QFT
    m: int
        approximation factor
    inverse: Bool
        True if the QFT should be applied inverse
    swaps: Bool
        True if Swap-Gates should be used to restore the order of the qubits
    MSB_first: Bool
        True if the first Qubit(the one with the lowest Value) is first
    
    
    
    """
    quantum_circuit = qiskit.QuantumCircuit(amount_qubits)
    if not MSB_first:
        measurement_qubits = list(reversed(range(amount_qubits)))
    else: 
        measurement_qubits = range(amount_qubits)
    if not inverse:
        index = 0
        for target in measurement_qubits:
            exponent = 2
            quantum_circuit.h(target)
            for control in measurement_qubits[index + 1:]:
                apply_unitary_operator(quantum_circuit, 2 * pi, control, target, exponent)
                exponent+=1
                if (m != -1) and exponent-1>m:
                    break 
            index+=1
        if swaps:
            for x in range(0, len(measurement_qubits)//2):
                quantum_circuit.swap(measurement_qubits[x], measurement_qubits[len(measurement_qubits) - 1 - x])
                x+=1
        quantum_circuit = quantum_circuit.to_gate()
        quantum_circuit.name = "  QFT"
    else:
        if swaps:
            for x in range(0, len(measurement_qubits)//2):
                quantum_circuit.swap(measurement_qubits[x], measurement_qubits[len(measurement_qubits) - 1 - x])
                x+=1
        index = len(measurement_qubits) -1
        for target in reversed(measurement_qubits):
            #quantum_circuit.barrier()
            exponent = len(measurement_qubits[index + 1:]) + 1
            for control in reversed(measurement_qubits[index + 1:]):
                if m == -1 or exponent-2 < m:
                    apply_unitary_operator(quantum_circuit, -2 * pi, control, target, exponent)
                exponent -= 1
            quantum_circuit.h(target)
            index-=1
        quantum_circuit = quantum_circuit.to_gate()
        quantum_circuit.name = "  iQFT"
    return quantum_circuit

def iterative_QFT(target: int, quantum_circuit ,k: range, controlling_qubit: int, m: int):
    """
    Applies one iteration of the iterative Quantum-Fourier-Transformation

    Parameters:
    quantum_circuit: qiskit.circuit
        The Quantum Circuite on which the iterative QFT will be applies on 
    controlling_qubit : int
        Position of the control Qubit, on which the iterative QFT will be applied on
    target: int
        Position of the classical controll Bits
    k: int
        Value of the controlling_qubit Qubit in the quantum_circuit, 
        defines which Phase-Gates will be applied on the controlling_qubit
    m: int
        approximation factor
    controlling_qubit : int
        Position of the classical 
    
    """
    if target != 0:
        for i in list(reversed(range(target))):
            if (not m == -1) and i+1>m:
                continue
            quantum_circuit.p((-pi/(2**(i+1))),controlling_qubit).c_if(k[target - (i+1)],1)
    quantum_circuit.h(controlling_qubit)