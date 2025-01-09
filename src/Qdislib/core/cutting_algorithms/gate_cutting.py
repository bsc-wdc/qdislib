from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import EstimatorV2 as Estimator

def expec_value(qc):
    observable = SparsePauliOp("Z" * 2)
    #params = [0.1] * qc.num_parameters
    
    exact_estimator = Estimator()
    # The circuit needs to be transpiled to the AerSimulator target
    pass_manager = generate_preset_pass_manager(3, AerSimulator())
    isa_circuit = pass_manager.run(qc)
    pub = (isa_circuit, observable)
    job = exact_estimator.run([pub])
    result = job.result()
    pub_result = result[0]
    exact_value = float(pub_result.data.evs)
    print(exact_value)
    return exact_value

