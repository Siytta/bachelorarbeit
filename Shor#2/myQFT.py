from qiskit import qiskit
from numpy import pi

def unitary_operator_c(theta):
    return (qiskit.circuit.library.standard_gates.PhaseGate(theta,str(theta))).control(1)

def apply_unitary_operator(quantum_circuit, theta, control, target, exponent):
    #Rechne theta * exponent da schneller als exponent-Fache anwendung von U
    quantum_circuit.append(unitary_operator_c(theta/(2**exponent)),[control, target])

def myQFT(quantum_circuit, measurement_qubits, inverse = False, swaps = False):
    if not inverse:
        index = 0
        for target in measurement_qubits:
            quantum_circuit.barrier()
            exponent = 2
            quantum_circuit.h(target)
            for control in measurement_qubits[index + 1:]:
                apply_unitary_operator(quantum_circuit, 2 * pi, control, target, exponent)
                exponent+=1
            index+=1
        quantum_circuit.barrier()
        if swaps:
            for x in range(0, len(measurement_qubits)//2):
                quantum_circuit.swap(measurement_qubits[x], measurement_qubits[len(measurement_qubits) - 1 - x])
                x+=1
        quantum_circuit.barrier()
    else:
        if swaps:
            quantum_circuit.barrier()
            for x in range(0, len(measurement_qubits)//2):
                quantum_circuit.swap(measurement_qubits[x], measurement_qubits[len(measurement_qubits) - 1 - x])
                x+=1
        index = len(measurement_qubits) -1
        for target in reversed(measurement_qubits):
            quantum_circuit.barrier()
            exponent = len(measurement_qubits[index + 1:]) + 1
            for control in reversed(measurement_qubits[index + 1:]):
                apply_unitary_operator(quantum_circuit, -2 * pi, control, target, exponent)
                exponent -= 1
            quantum_circuit.h(target)
            index-=1
        quantum_circuit.barrier()