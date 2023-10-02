from qiskit import qiskit
from numpy import pi
from myQFT import QFT_Gate, iterative_QFT
from binToDez import binToDez, dezToBin, mod_exp
from gcd import gcd

def P_Gate(theta: float, controls: int = 0) -> qiskit.circuit.library.standard_gates.PhaseGate:
    """Returns a Phase Gate with *controls* many controlled inputs
        A single Qubit Phase Gate looks like this:
        ( 1        0     )
        ( 0  e^(i*theta) )

    Parameters:
    theta: float
        determines the phase shift
    controls: int
        defines the amount of control bits. Default = 0

    Returns: Phase Gate with *controls* many controlled inputs
    """
    if controls > 0:
        return (qiskit.circuit.library.standard_gates.PhaseGate(theta,str(theta))).control(controls)
    else:
        return (qiskit.circuit.library.standard_gates.PhaseGate(theta,str(theta)))

def A_Gate(addend_bin: list[int]) -> qiskit.circuit.gate:
    """Returns the complete Addition Circuit as an Gate for the quantum register + classial register Addition as descriped by Beauregard(Stêphane Beauregard, Circuit for Shor’s algorithm using 2n+3 qubits, https://arxiv.org/pdf/quant-ph/0205095.pdf) 
        Infact this is a combination of len(a_bin) many Phase Gates with added up theta's depending/controlled by the binary representation of the classically summand(a)
        The Gate operates on n+1 Qubits B contraining an n bit number to prevent overflow. So the input(before QFT) for the MSB of B should be (an extra) |0>

    Parameters:
    a_bin: list[int]
        Binary representation of the classically summand a. The first element of the list needs to be the Least significant Bit & the last element the Most significant Bit. Example [0,0,0,1] to represent the decimal number 8.

    Returns: Complete classic + quantum Addition Circuite as an gate
    """
    A_Gate = qiskit.QuantumCircuit(len(addend_bin))
    theta_list = [0.0]*len(addend_bin)
    for target_bit in range(len(addend_bin)):
        exponent = 1
        for control_bit in reversed(range(target_bit+1)):
            if addend_bin[control_bit] == 1:
                theta_list[target_bit]+= 2*pi/(2**(exponent))
            exponent+=1
    for qubit_index in range(len(addend_bin)):
        A_Gate.append(P_Gate(theta_list[qubit_index]),[qubit_index])
    A_Gate = A_Gate.to_gate()
    A_Gate.name = "  Add(" + str (binToDez(addend_bin) )+ ")"
    return A_Gate

def S_Gate(subtrahend_bin: list[int]) -> qiskit.circuit.gate:
    S_Gate = A_Gate(subtrahend_bin).inverse()
    S_Gate.name = "  Sub(" + str (binToDez(subtrahend_bin) )+ ")"
    return S_Gate

def modular_adder_gate(a_bin: list[int],N_bin: list[int], m: int) -> qiskit.circuit.gate:
    """Returns the modular adder gate. The Gate requiers 2 + len(a_bin) + 1 Qubit. The first two Qubits are both control Qubits, the last one is an anciallary Qubit which needs to be in state |0> in the beginning.

    Parameters:
    a_bin: list[int]
        Binary representation of the classically summand a. The first element of the list needs to be the Least significant Bit & the last element the Most significant Bit. Example [0,0,0,1] to represent the decimal number 8.
    N_bin: list[int]
        Binary representation of the Moduland N / the Prime Factor. The first element of the list needs to be the Least significant Bit & the last element the Most significant Bit. Example [0,0,0,1] to represent the decimal number 8.

    Returns:  modular adder gate
    """
    c_qbits = [0,1]
    b_qbits = list(range(2, len(a_bin)+2))
    cond_qbit = len(a_bin)+2
    m_a_g = qiskit.QuantumCircuit(2 + len(a_bin) + 1) # 2 Control Qubits + qubits for the input + ancillary qubit
    m_a_g.append(A_Gate(a_bin).control(2), c_qbits + b_qbits)
    m_a_g.append(S_Gate(N_bin),b_qbits)
    m_a_g.append(QFT_Gate(len(a_bin),inverse = True, MSB_first = False, m = m), b_qbits)
    m_a_g.cnot(b_qbits[-1],cond_qbit)
    m_a_g.append(QFT_Gate(len(a_bin),inverse = False, MSB_first = False, m = m), b_qbits)
    m_a_g.append(A_Gate(N_bin).control(1), [cond_qbit] + b_qbits)
    m_a_g.append(S_Gate(a_bin).control(2), c_qbits + b_qbits)
    m_a_g.append(QFT_Gate(len(a_bin),inverse = True, MSB_first = False, m = m), b_qbits)
    m_a_g.x(b_qbits[-1])
    m_a_g.cnot(b_qbits[-1],cond_qbit)
    m_a_g.x(b_qbits[-1])
    m_a_g.append(QFT_Gate(len(a_bin),inverse = False, MSB_first = False, m = m), b_qbits)
    m_a_g.append(A_Gate(a_bin).control(2), c_qbits + b_qbits)
    m_a_g = m_a_g.to_gate()
    m_a_g.name = "Add " + str(binToDez(a_bin)) + " Mod " + str(binToDez(N_bin))
    return m_a_g

def cmult_gate(x_bits_amount: int,a_bin: list[int],N_bin: list[int], m : int) -> qiskit.circuit.gate:
    """Returns the controlled multiplier gate. 
    The gate requiers 1 + x_bits_amount + len(a_bit) + 1
    The first Qubit is the controll qubit. The second x_bits_amount Qubits are for the Qubits of the (number) x. The a_bit should be the equal lenght of the b input qubits. The last Qubits is the ancillary Gate of the underlaying modular_adder_gates

    Parameters:
    x_bits_amount: int
       amount as an int of qubits for x
    a_bin: list[int]
        a_bin for the calculaton of the gate. a_bin should be as long as N_bin and as many binary bits as b qubits
    N_bin: list[int]
        N_bin for the calculaton of the gate. N_bin should be as long as a_bin
  
    Returns: controlled multiplier gate
    """
    a_dez = binToDez(a_bin)
    N_dez = binToDez(N_bin)
    c_qbit = [0]
    x_qbits =  list(range(1, x_bits_amount + 1))
    b_qbits = list(range(1+x_bits_amount, 1 + x_bits_amount + len(a_bin)))
    cond_qbit = [1 + x_bits_amount + len(a_bin)]
    cmult_gate = qiskit.QuantumCircuit(1 + x_bits_amount + len(a_bin) + 1)
    cmult_gate.append(QFT_Gate(len(b_qbits),inverse = False, MSB_first = False, m = m), b_qbits)
    for i in x_qbits:
        a_i = ((2**(i - 1)) *  a_dez) % N_dez
        a_i_bin = dezToBin(a_i, len(a_bin))
        cmult_gate.append(modular_adder_gate(a_i_bin, N_bin, m = m), c_qbit + [i] + b_qbits + cond_qbit )
    cmult_gate.append(QFT_Gate(len(b_qbits),inverse = True, MSB_first = False, m = m), b_qbits)
    cmult_gate = cmult_gate.to_gate()
    cmult_gate.name = "cmult " + str(a_dez) + " Mod " + str(N_dez)
    return cmult_gate

def inv_cmult_gate(x_bits_amount: int,a_bin: list[int],N_bin: list[int], m : int) -> qiskit.circuit.gate:
    """Returns the controlled multiplier gate. 
    The gate requiers 1 + x_bits_amount + len(a_bit) + 1
    The first Qubit is the controll qubit. The second x_bits_amount Qubits are for the Qubits of the (number) x. The a_bit should be the equal lenght of the b input qubits. The last Qubits is the ancillary Gate of the underlaying modular_adder_gates

    Parameters:
    x_bits_amount: int
       amount as an int of qubits for x
    a_bin: list[int]
        a_bin for the calculaton of the gate. a_bin should be as long as N_bin and as many binary bits as b qubits
    N_bin: list[int]
        N_bin for the calculaton of the gate. N_bin should be as long as a_bin
  
    Returns: controlled multiplier gate
    """
    inv_cmult_gate = cmult_gate(x_bits_amount,a_bin,N_bin, m = m).inverse()
    inv_cmult_gate.name = "inv cmult " + str(binToDez(a_bin)) + " Mod " + str(binToDez(N_bin))
    return inv_cmult_gate

def U_gate(x_bits_amount: int,a_bin: list[int],N_bin: list[int], m : int) -> qiskit.circuit.gate:
    """Returns the U_a_gate. 

    Qubits:
    Number 0: Controll Qubit
    Number 1 to x_bits_amount: |x>
    Number (x_bits_amount+1) to (x_bits_amount+1+len(a_bin)): |0>
    Number Last: ancillary Gate of the underlaying controlled multiplier gates (|0>)

    Parameters:
    x_bits_amount: int
       amount as an int of qubits for x
    a_bin: list[int]
        a_bin for the calculaton of the gate. a_bin should be as long as N_bin and as many binary bits as b qubits
    N_bin: list[int]
        N_bin for the calculaton of the gate. N_bin should be as long as a_bin

    Returns:  controlled U-Gate
    """
    c_qbit = [0]
    x_qbits =  list(range(1, x_bits_amount + 1))
    b_qbits = list(range(1+x_bits_amount, 1 + x_bits_amount + len(a_bin)))
    cond_qbit = [1 + x_bits_amount + len(a_bin)]
    full_range = c_qbit + x_qbits + b_qbits + cond_qbit
    U_a_gate = qiskit.QuantumCircuit(len(full_range))
    U_a_gate.append(cmult_gate(x_bits_amount,a_bin ,N_bin, m = m), full_range)
    for i in (range(x_bits_amount)):
        U_a_gate.cswap(c_qbit,x_qbits[i],b_qbits[i])
    a_inv_bin = dezToBin(gcd(binToDez(a_bin),binToDez(N_bin)),len(a_bin))
    U_a_gate.append(inv_cmult_gate(x_bits_amount,a_inv_bin ,N_bin, m = m), full_range)
    U_a_gate = U_a_gate.to_gate()
    U_a_gate.name = "  U_" + str(binToDez(a_bin))
    return U_a_gate

def Shor(a: int, N: int, quantum_circuit : qiskit.circuit ,number_shots: int = 1 ,backend: str = 'aer_simulator',m = -1, k = -1)  -> tuple:
    """
    Uses the Shor-Algorithm to find the primefactors of N with the non iterative Quantum-Phase-Estimation

    Parameters: 
    a:int
        in a which is not a factor of N
    N:int
        The Number of which the prime factors will be found
    quantum_circuit: qiskit.circuit
        A qiskit circuit. If this Parameter is = 0 the Function will create a circuit for Shor otherwise 
        the Function will only run the circuit in the Parameter
    number_shots: int
        amount of messurings/runs
    backend: str
        Name of the used backend system. Can be a IBM Quantum Ressource or the Simulator
    m : int 
        approximaion factor of the QFT
    k : int 
        precision (Amount of controll Qubits)
    """
    if quantum_circuit == 0:
        n = N.bit_length()
        c_qbits = range(k)
        ev_qbits = list(range(len(c_qbits),2*n+2+k))
        quantum_circuit = qiskit.QuantumCircuit(2*n+2+k,len(c_qbits)) # 4*n+2 size for complete input bits, 2*n (N_COUNT) size of the control qubits which will also be measured (if k = 2n)
        for c_bit in c_qbits:
            quantum_circuit.h(c_bit)
        quantum_circuit.x(ev_qbits[0])
        for c_bit in c_qbits:
            quantum_circuit.append(U_gate(n,dezToBin(mod_exp(a, 2**c_bit, N),n+1),dezToBin(N,n+1), m = m),[c_bit] + ev_qbits)
        quantum_circuit.append(QFT_Gate(len(c_qbits),inverse = True, MSB_first = False,swaps = True, m = m ), c_qbits)
        quantum_circuit.measure(c_qbits,c_qbits)
    simulator = qiskit.Aer.get_backend(backend)
    sim_result = qiskit.execute(quantum_circuit, backend=simulator, shots=number_shots).result()
    return sim_result.get_counts(quantum_circuit), quantum_circuit

def Shor_sequential(a: int, N: int, quantum_circuit : qiskit.circuit,number_shots: int = 1, backend: str = 'aer_simulator',  k: int = -1, m : int = -1) -> tuple:
    """
    Uses the Shor-Algorithm to find the primefactors of N with the iterative Quantum-Phase-Estimation
    Parameters: 
    a:int
        in a which is not a factor of N
    N:int
        The Number of which the prime factors will be found
    quantum_circuit: qiskit.circuit
        A qiskit circuit. If this Parameter is = 0 the Function will create a circuit for Shor otherwise 
        the Function will only run the circuit in the Parameter
    number_shots: int
        amount of messurings/runs
    backend: str
        Name of the used backend system. Can be a IBM Quantum Ressource or the Simulator
    m : int 
        approximaion factor of the QFT
    k : int 
        precision (Amount of episods of the iterative Quantum-Phase-Estimation)
    """
    if quantum_circuit == 0:
        n = N.bit_length()
        c_qbits = range(k)
        controlling_qubit = qiskit.QuantumRegister(1, "c_qubit")
        qreg = qiskit.QuantumRegister(2*n+2, "qreg")
        creg = qiskit.ClassicalRegister(len(c_qbits), "creg")
        quantum_circuit = qiskit.QuantumCircuit(controlling_qubit, qreg,creg) 
        quantum_circuit.x(1)
        for control in c_qbits:
            quantum_circuit.h(0)
            quantum_circuit.append(U_gate(n,dezToBin(mod_exp(a, 2**(len(c_qbits) - control - 1), N),n+1),dezToBin(N,n+1), m = m),[0] + list(range(1,2*n+3)))
            iterative_QFT(control,quantum_circuit,creg,controlling_qubit, m = m)
            quantum_circuit.measure(controlling_qubit,creg[control])
            quantum_circuit.x(0).c_if(creg[control],1)
    simulator = qiskit.Aer.get_backend(backend)
    sim_result = qiskit.execute(quantum_circuit, backend=simulator, shots=number_shots).result()
    return sim_result.get_counts(quantum_circuit), quantum_circuit