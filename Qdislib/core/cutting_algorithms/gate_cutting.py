#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# -*- coding: utf-8 -*-

"""Gate cutting algorithms."""

import networkx as nx
import math
import numpy as np

from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import EstimatorV2 as Estimator
from qiskit.primitives import BackendSamplerV2

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.parameter import COLLECTION_IN
from pycompss.api.parameter import COLLECTION_OUT
from pycompss.api.constraint import constraint
from pycompss.api.implement import implement
from pycompss.api.api import compss_barrier

import Qdislib
from Qdislib.utils.graph_qibo import _update_qubits
from Qdislib.utils.graph_qibo import _remove_red_edges
from Qdislib.utils.graph_qibo import _update_qubits, circuit_qibo_to_dag, plot_dag, _max_qubits_graph
from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit,circuit_qiskit_to_dag, _dag_to_circuit_qiskit_subcircuits
from Qdislib.utils.graph_qibo import dag_to_circuit_qibo,circuit_qibo_to_dag, _dag_to_circuit_qibo_subcircuits
from Qdislib.core.cutting_algorithms.wire_cutting import _sum_results
from Qdislib.core.find_cut.find_cut import find_nodes_with_qubit

from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2 as Sampler, EstimatorV2 as Estimator


import qiskit
import qibo
from qibo import gates
from qiskit import transpile


import pickle
import os
import time
import typing

def gate_cutting(
    dag: typing.Any,
    gates_cut: typing.List[typing.Any],
    observables: typing.Optional[typing.Any] = None,
    shots: int = 1024,
    method: str = "automatic",
    sync: bool = True,
    gpu: bool = False,
    gpu_min_qubits: typing.Optional[int] = None,
    qpu: bool = False,
    qpu_dict: typing.Optional[dict] = None,
    max_time_ibm: int = 1800
):
    """
    Apply gate cutting to a quantum circuit to enable distributed smaller executions of subcircuits.

    This function partitions a quantum circuit by cutting a list of specified two-qubit gates,
    evaluates the resulting subcircuits individually (optionally on CPU/GPU or QPU resources),
    and reconstructs the final expectation value from the partial results.
    It supports circuits represented in Qiskit or Qibo formats, as well as their internal DAG representations.

    If no gates are cut, the function directly computes the expectation value of the original circuit.

    Parameters
    ----------
    dag : qiskit.QuantumCircuit or qibo.models.Circuit or networkx.Graph
        The quantum circuit to cut and evaluate. Can be a Qiskit circuit, a Qibo circuit,
        or an already constructed DAG (as a NetworkX graph).
    gates_cut : list
        List of two-qubit gates (edges in the DAG) to cut. Each gate is typically identified by its node ID.
    observables : list or dict, optional
        Observables to measure. If specified, the circuit is transformed into the measurement basis.
        Observables may include the identity ("I").
    shots : int, default=1024
        Number of measurement shots to use for circuit evaluations.
    method : str, default='automatic'
        Simulation method to use. Possible values are backend-dependent (e.g., "statevector", "qasm_simulator").
    sync : bool, default=True
        Whether to synchronize and wait for all computations to finish (useful in distributed environments like PyCOMPSs).
    return_subcircuits : bool, default=False
        If True, also return the list of subcircuits generated after cutting.
    gpu : bool, default=False
        If True, attempt to run circuit evaluations on GPU simulators when possible.
    gpu_min_qubits : int, optional
        Minimum number of qubits needed to offload execution to GPU simulators. Default is 0.
    qpu : bool, default=False
        If True, attempt to run circuit evaluations on real Quantum Processing Units (QPUs).
    qpu_dict : dict, optional
        Dictionary specifying QPU backends and their maximum supported qubits.
        Example: `{"IBM_Quantum": 153, "MN_Ona": 5}`.
    max_time_ibm : int, default=1800
        Maximum allowed execution time (in seconds) for IBM Quantum batches.

    Returns
    -------
    final_recons : float
        The reconstructed expectation value after gate cutting and execution of subcrcuits.

    Notes
    -----
    - If the circuit consists of multiple disconnected components, each component is processed separately.
    - Depending on the device availability and size of the subcircuits, they are evaluated on CPU, GPU, or QPU.
    - Results are automatically scaled by the number of gates cut (following standard gate cutting reconstruction rules).
    - GPU and QPU support is optional and requires appropriate setup (e.g., Qiskit Aer with GPU support, IBM Q account credentials).

    Examples
    --------
    Cutting a simple Qiskit circuit:

    .. code-block:: python

        from qiskit import QuantumCircuit
        from mymodule import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)     # Gate "H_1"
        qc.cz(0, 1) # Gate "CZ_2"
        qc.h(1)     # Gate "H_3"

        # Define the gate to cut (e.g., the CZ gate)
        gates_to_cut = ["CZ_2"]  # You must identify the gate in order of definition

        reconstruction = gate_cutting(qc, gates_to_cut, shots=2048)

    """

    if qpu and "IBM_Quantum" in qpu_dict:
        batch, backend = _check_ibm_qc(max_time_ibm=max_time_ibm)
    else:
        batch, backend = None, None

    if observables:
        dag = _change_basis(dag,observables)

    if qpu_dict is None:
        qpu_dict = {}
    if gpu_min_qubits is None:
        gpu_min_qubits = 0

    if type(dag) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        if observables:
            if "I" in observables:
                dag = circuit_qiskit_to_dag(dag, obs_I=observables)
            else:
                dag = circuit_qiskit_to_dag(dag)
        else:
            dag = circuit_qiskit_to_dag(dag)
    elif type(dag) == qibo.models.Circuit:
        if observables:
            if "I" in observables:
                dag = circuit_qibo_to_dag(dag,obs_I=observables)
            else:
                dag = circuit_qibo_to_dag(dag)
        else:
            dag = circuit_qibo_to_dag(dag)
    else:
        dag = dag

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in nx.connected_components(dag.to_undirected())
        ]
        results = []
        for s in S:
            #num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in gates_cut:
                if s.has_node(c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                graphs = _execute_gate_cutting(dag, tmp_cuts, shots=shots, method=method, gpu=gpu, gpu_min_qubits=gpu_min_qubits,qpu=qpu ,qpu_dict=qpu_dict, batch=batch, backend=backend)
                graphs = _sum_results(graphs)
                results.append(graphs)  # 1/(2**len(tmp_cuts))*sum(graphs)
                #subcircuits.append(subcircuit)
            else:
                _max_qubit = _max_qubits_graph(s)
                s_new, highest_qubit = _update_qubits(s)
                subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
                if qpu and 'MN_Ona' in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                    expected_value = _expec_value_qibo_qpu(subcirc,shots=shots,method=method)
                elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                    expected_value = _expec_value_qiskit_gpu(subcirc, shots=shots, method=method)
                elif qpu and "IBM_Quantum" in qpu_dict and qpu_dict["IBM_Quantum"] >= _max_qubit:
                    expected_value = _expec_value_qiskit_qpu(subcirc, shots=shots, batch=batch, backend=backend)
                else:
                    expected_value = _expec_value_qiskit(subcirc, shots=shots, method=method)
                results.append(expected_value)
        if sync:
            results = compss_wait_on(results)
        # Consider gate cutting within wire cutting:
        # if gate_cutting:
        #     final_recons = 1/2*sum(results)
        # else:
        final_recons = 1 / (2 ** len(gates_cut)) * math.prod(results)
        return final_recons
    else:
        #num_qubits = _max_qubit(dag)
        if gates_cut:
            results = _execute_gate_cutting(dag, gates_cut, shots=shots, method=method, gpu=gpu, gpu_min_qubits=gpu_min_qubits, qpu=qpu, qpu_dict=qpu_dict, batch=batch, backend=backend)
            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(gates_cut)) * sum(results)
        else:
            _max_qubit = _max_qubits_graph(dag)
            s_new, highest_qubit = _update_qubits(dag)
            subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
            if qpu and 'MN_Ona' in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                final_recons = _expec_value_qibo_qpu(subcirc,shots=shots,method=method)
            elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                final_recons = _expec_value_qiskit_gpu(subcirc,shots=shots,method=method)
            elif qpu and "IBM_Quantum" in qpu_dict and qpu_dict["IBM_Quantum"] >= _max_qubit:
                final_recons = _expec_value_qiskit_qpu(subcirc, shots=shots, batch=batch, backend=backend)
            else:
                final_recons = _expec_value_qiskit(subcirc, shots,method=method)
        return final_recons

def _generate_cut(dag, gates_cut):
    dag_copy = dag.copy()
    dag_copy = _remove_red_edges(dag_copy)
    for index, gate_name in enumerate(gates_cut, start=1):
        target_qubits = dag_copy.nodes[gate_name]["qubits"]
        # Find predecessor node with qubit 1
        pred_0 = find_nodes_with_qubit(dag_copy, gate_name, qubit=target_qubits[0], direction='predecessor')

        # Find predecessor node with qubit 2
        pred_1 = find_nodes_with_qubit(dag_copy, gate_name, qubit=target_qubits[1], direction='predecessor')

        # Find successor node with qubit 1
        succ_0 = find_nodes_with_qubit(dag_copy, gate_name, qubit=target_qubits[0], direction='successor')

        # Find successor node with qubit 2
        succ_1 = find_nodes_with_qubit(dag_copy, gate_name, qubit=target_qubits[1], direction='successor')

        # Output the results
        #print(f"Predecessor nodes with qubit {target_qubits[0]}: {pred_0}")
        #print(f"Predecessor nodes with qubit {target_qubits[1]}: {pred_1}")
        #print(f"Successor nodes with qubit {target_qubits[0]}: {succ_0}")
        #print(f"Successor nodes with qubit {target_qubits[1]}: {succ_1}")


        dag_copy.remove_node(gate_name)

        dag_copy.add_node(f"SUBS1_{index}", gate='S', qubits=(target_qubits[0],), parameters=())
        dag_copy.add_node(f"SUBS2_{index}", gate='S', qubits=(target_qubits[1],), parameters=())

        if pred_0:
            dag_copy.add_edge(pred_0[0], f"SUBS1_{index}", color="blue")

        if succ_0:
            dag_copy.add_edge(f"SUBS1_{index}", succ_0[0], color="blue")

        if pred_1:
            dag_copy.add_edge(pred_1[0], f"SUBS2_{index}", color="blue")

        if succ_1:
            dag_copy.add_edge( f"SUBS2_{index}", succ_1[0], color="blue")


    return dag_copy


def _decimal_to_base6(num):
    if num == 0:
        return "0"

    base6 = ""
    while num > 0:
        base6 = str(num % 6) + base6
        num //= 6
    return base6


#@constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, graph_components=COLLECTION_OUT)
def _generate_gate_cutting(updated_dag, gates_cut, index, graph_components):

    base6_rep = _decimal_to_base6(index)
    base6_rep = base6_rep.zfill(len(gates_cut))
    list_substitutions = list(map(int, base6_rep))

    #print(list_substitutions)
    for idx2, idx in enumerate(list_substitutions,start=1):
        list_succ1 = list(updated_dag.succ[f'SUBS1_{idx2}'])
        list_pred1 = list(updated_dag.pred[f'SUBS1_{idx2}'])

        list_succ2 = list(updated_dag.succ[f'SUBS2_{idx2}'])
        list_pred2 = list(updated_dag.pred[f'SUBS2_{idx2}'])

        # 1 - Rz(-pi/2) -- Rz(-pi/2)
        if idx == 0:
            updated_dag.nodes[f'SUBS1_{idx2}']['gate'] = 'rz'
            updated_dag.nodes[f'SUBS1_{idx2}']['parameters'] = [-np.pi/2]
            updated_dag.nodes[f'SUBS2_{idx2}']['gate'] = 'rz'
            updated_dag.nodes[f'SUBS2_{idx2}']['parameters'] = [-np.pi/2]


        #2 - Z Rz(-pi/2) -- Z Rz(-pi/2)
        elif idx == 1:
            updated_dag.nodes[f'SUBS1_{idx2}']['gate'] = 'z'
            updated_dag.add_node(f'SUBS11_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS1_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f'SUBS1_{idx2}', succ1)
                updated_dag.add_edge(f'SUBS11_{idx2}', succ1, color="blue")

            updated_dag.add_edge(f'SUBS1_{idx2}', f'SUBS11_{idx2}', color="blue")

            updated_dag.nodes[f'SUBS2_{idx2}']['gate'] = 'z'
            updated_dag.add_node(f'SUBS22_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS2_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f'SUBS2_{idx2}', succ2)
                updated_dag.add_edge(f'SUBS22_{idx2}', succ2, color="blue")

            updated_dag.add_edge(f'SUBS2_{idx2}', f'SUBS22_{idx2}', color="blue")


        #3 - MEASURE Rz(-pi/2) -- Rz(-pi)
        elif idx == 2:
            updated_dag.nodes[f'SUBS1_{idx2}']['gate'] = 'measure'
            updated_dag.add_node(f'SUBS11_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS1_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f'SUBS1_{idx2}', succ1)
                updated_dag.add_edge(f'SUBS11_{idx2}', succ1, color="blue")

            updated_dag.add_edge(f'SUBS1_{idx2}', f'SUBS11_{idx2}', color="blue")

            updated_dag.nodes[f'SUBS2_{idx2}']['gate'] = 'rz'
            updated_dag.nodes[f'SUBS2_{idx2}']['parameters'] = [-np.pi]

        #4 - MEASURE Rz(-pi/2) -- RES
        elif idx == 3:
            updated_dag.nodes[f'SUBS1_{idx2}']['gate'] = 'measure'
            updated_dag.add_node(f'SUBS11_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS1_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f'SUBS1_{idx2}', succ1)
                updated_dag.add_edge(f'SUBS11_{idx2}', succ1, color="blue")

            updated_dag.add_edge(f'SUBS1_{idx2}', f'SUBS11_{idx2}', color="blue")

            updated_dag.remove_node(f'SUBS2_{idx2}')

            if list_succ2 and list_pred2:
                succ2 = list_succ2[0]
                pred2 = list_pred2[0]

                updated_dag.add_edge(pred2, succ2)

        #5 - Rz(-pi) -- MEASURE Rz(-pi/2)
        elif idx == 4:
            updated_dag.nodes[f'SUBS1_{idx2}']['gate'] = 'rz'
            updated_dag.nodes[f'SUBS1_{idx2}']['parameters'] = [-np.pi]

            updated_dag.nodes[f'SUBS2_{idx2}']['gate'] = 'measure'
            updated_dag.add_node(f'SUBS22_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS2_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f'SUBS2_{idx2}', succ2)
                updated_dag.add_edge(f'SUBS22_{idx2}', succ2, color="blue")

            updated_dag.add_edge(f'SUBS2_{idx2}', f'SUBS22_{idx2}', color="blue")

        #6 - RES -- MEASURE Rz(-pi/2)
        elif idx == 5:
            updated_dag.remove_node(f'SUBS1_{idx2}')

            if list_succ1 and list_pred1:
                succ1 = list_succ1[0]
                pred1 = list_pred1[0]
                updated_dag.add_edge(pred1, succ1, color="blue")

            updated_dag.nodes[f'SUBS2_{idx2}']['gate'] = 'measure'
            updated_dag.add_node(f'SUBS22_{idx2}', gate='rz', qubits=updated_dag.nodes[f'SUBS2_{idx2}'].get('qubits'), parameters=([-np.pi/2]))

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f'SUBS2_{idx2}', succ2)
                updated_dag.add_edge(f'SUBS22_{idx2}', succ2, color="blue")

            updated_dag.add_edge(f'SUBS2_{idx2}', f'SUBS22_{idx2}', color="blue")


        else:
            raise TypeError


    #updated_dag = _remove_red_edges(updated_dag)
    for i, c in enumerate(nx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")

    return updated_dag

#@constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign_gate_cutting(expectation_value, index):
    expectation_value = [x for x in expectation_value if x is not None]
    expectation_value = math.prod(expectation_value)
    number = index

    _change_sign = False

    while number != 0:
        digit = number % 6  # Get the last digit
        if digit in {3, 5}:  # Check if the digit is 3, 5, or 7
            _change_sign = not _change_sign  # Flip the sign change flag
        number //= 6 # Move to the next digit
    #print(_change_sign)

    # If _change_sign is True, we flip the sign of the original number
    if _change_sign:
        return -expectation_value
    else:
        return expectation_value


def _execute_gate_cutting(dag, gates_cut, shots=10000, method='automatic', gpu=False, gpu_min_qubits=None, qpu=None, qpu_dict=None, batch=None, backend=None):
    new_dag = _generate_cut(dag,gates_cut)


    S = [new_dag.subgraph(c).copy() for c in nx.connected_components(new_dag.to_undirected())]

    max_qubit_list = []
    for s in S:
        _max_qubit = _max_qubits_graph(s)
        max_qubit_list.append(_max_qubit)
    #plot_dag(new_dag)

    '''list_dags = []
    for i in range(6**len(gates_cut)):
        list_dags.append(new_dag.copy())'''

    #print(list_dags)

    bool_mult_qpu = True

    reconstruction = []
    for index in range(6**len(gates_cut)):

        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(
            copy_graph.to_undirected()
        )

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())
        #print(len(graph_components))

        graph = _generate_gate_cutting(copy_graph, gates_cut,index, graph_components)
        #graph = compss_wait_on(graph)
        #plot_dag(graph)

        #for i in graph_components:
        #    pass
            #plot_dag(i)

        exp_values = []
        for i, s in enumerate(graph_components):
            #print(s)
            #if s.nodes():
            _max_qubit = max_qubit_list[i]
            s_new, highest_qubit = _update_qubits(s)
            if qpu and 'MN_Ona' in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
            else:
                subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
            #subcircuits.append(subcirc)
            #subcirc = compss_wait_on(subcirc)
            #print(subcirc.draw())
            if qpu and 'MN_Ona' in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                exp = _expec_value_qibo_qpu(subcirc,shots=shots,method=method)
            elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                exp = _expec_value_qiskit_gpu(subcirc,shots=shots,method=method)
            elif qpu and "IBM_Quantum" in qpu_dict and qpu_dict["IBM_Quantum"] >= _max_qubit:
                exp = _expec_value_qiskit_qpu(subcirc, shots=shots, batch=batch, backend=backend)
            else:
                exp = _expec_value_qiskit(subcirc, shots=shots,method=method)
            #print(exp)
            exp_values.append(exp)



        #exp_values = compss_wait_on(exp_values)
        #print(exp_values)
        exp = _change_sign_gate_cutting(exp_values, index)
        reconstruction.append(exp)

        if bool_mult_qpu:
            bool_mult_qpu=False
        else:
            bool_mult_qpu= True
        '''if index % 5000 == 0:
            compss_barrier()'''

    #print(reconstruction)
    #if subcircuits is not []:
    #    subcircuits = [item for sublist in subcircuits for item in sublist if item is not None]
    return  reconstruction

@task(returns=1)
def _expec_value_qibo_qpu(subcirc, shots=1024, method='numpy'):

    if subcirc is None:
        return None

    tmp = subcirc[1]
    subcirc = subcirc[0]

    #print(subcirc.draw())
    if tmp:
        obs_I = tmp
    else:
        obs_I = None
    observables = ["Z"] * subcirc.nqubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    #print(obs_I)

    subcirc.add(gates.M(*range(subcirc.nqubits)))

    observables = "".join(observables)
    #print(observables)


    '''counter = 1  # Initialize counter
    while os.path.exists(f"/home/bsc/bsc019635/ona_proves/subcircuits/circuit_{counter}.pkl"):  # Check if file already exists
        counter += 1'''

    unique_code = str(time.time_ns())

    circuit_filename = f"/home/bsc/bsc019635/ona_proves/subcircuits/circuit_{unique_code}.pkl"
    result_filename = f"/home/bsc/bsc019635/ona_proves/subcircuits/result_{unique_code}.pkl"

    # Save the circuit to the unique file
    with open(circuit_filename, "wb") as f:
        pickle.dump(subcirc, f)

    print(f"Circuit saved: {circuit_filename}, waiting for {result_filename}...")

    while not os.path.exists(result_filename):
        time.sleep(1)  # Check every second

    # Load the circuit from the file
    with open(result_filename, "rb") as f:
        result = pickle.load(f)
        print(f"Received result from {result_filename}: {result}")

    freq = result.counts()

    #os.remove(circuit_filename)  # Remove circuit file after processing
    os.remove(result_filename)   # Remove result file after reading

    expectation_value = 0
    for key, value in freq.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value



@constraint(processors=[{"processorType": "CPU", "computingUnits": "1"}])
@task(returns=1)
def _expec_value_qiskit(qc, shots=1024, method='automatic'):

    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    #print(subcirc.draw())
    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    #print(obs_I)

    observables = "".join(observables)
    #print(observables)

    qc, _ = qc
    num_qubits = qc.num_qubits
    observable = SparsePauliOp(observables)
    #params = [0.1] * qc.num_parameters

    qc.measure_all()

    import subprocess
    import os


    try:
        subprocess.check_output('nvidia-smi')
        print('Nvidia GPU detected!')
    except:
        print('No Nvidia GPU in system!')

    '''target_gpus = [os.environ["COMPSS_BINDED_GPUS"]]
    print("TGPU ", target_gpus)
    print("CUDA ", [os.environ["CUDA_VISIBLE_DEVICES"]])
    print("DEVICES ", AerSimulator().available_devices())
    print("METHODS ", AerSimulator().available_methods())'''

    #sampler = BackendSamplerV2(backend = AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    #pass_manager = generate_preset_pass_manager(0, AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    #simulator = AerSimulator(device='GPU', method=method, cuStateVec_enable=True)
    #except Exception: # this command not being found can raise quite a few different errors depending on the configuration
    #print('No Nvidia GPU in system!')
    #sampler = BackendSamplerV2(backend = AerSimulator())
    #pass_manager = generate_preset_pass_manager(0, AerSimulator())
    print(method)
    simulator = AerSimulator(device='CPU', method=method, max_parallel_threads=1, mps_omp_threads=1, mps_parallel_threshold=1)


    #circ = transpile(qc, simulator)


    #exact_estimator = Estimator()
    # The circuit needs to be transpiled to the AerSimulator target

    #isa_circuit = pass_manager.run(qc)
    #pub = (isa_circuit, observable)
    #job = sampler.run([isa_circuit], shots=shots)
    job = simulator.run(qc, shots=shots)
    result = job.result()
    print(f'backend: {result.backend_name}')
    #pub_result = result[0]
    #print("METADATA ", result)
    #counts = pub_result.data.meas.get_counts()
    counts = result.get_counts()
    print(counts)
    expectation_value = 0
    for key, value in counts.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value


@task(returns=1)
def _expec_value_qiskit_qpu(qc, shots=1024, method='automatic', batch=None, backend=None):

    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    #print(subcirc.draw())
    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    #print(obs_I)

    observables = "".join(observables)
    #print(observables)

    qc, _ = qc
    num_qubits = qc.num_qubits
    observable = SparsePauliOp(observables)
    #params = [0.1] * qc.num_parameters

    qc.measure_all()

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)

    sampler = Sampler(mode=batch)
    job = sampler.run([(isa_circuit)])
    print(f"job id: {job.job_id()}")
    job_result = job.result()
    print(job_result)
    pub_result = job_result[0].data.meas.get_counts()
    print(pub_result)

    expectation_value = 0
    for key, value in pub_result.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value



#@implement(source_class="Qdislib.core.cutting_algorithms.gate_cutting", method="_expec_value_qiskit")
@constraint(processors=[{"processorType": "CPU", "computingUnits": "1"},{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1)
def _expec_value_qiskit_gpu_implements(qc, shots=1024, method='automatic'):
    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    #print(subcirc.draw())
    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    #print(obs_I)

    observables = "".join(observables)
    #print(observables)

    qc, _ = qc
    num_qubits = qc.num_qubits
    observable = SparsePauliOp(observables)
    #params = [0.1] * qc.num_parameters

    qc.measure_all()

    import subprocess
    import os


    try:
        subprocess.check_output('nvidia-smi')
        print('Nvidia GPU detected!')
    except Exception: # this command not being found can raise quite a few different errors depending on the configuration
        print('No Nvidia GPU in system!')


    target_gpus = [os.environ["COMPSS_BINDED_GPUS"]]
    print("TGPU ", target_gpus)
    print("CUDA ", [os.environ["CUDA_VISIBLE_DEVICES"]])
    print("DEVICES ", AerSimulator().available_devices())
    print("METHODS ", AerSimulator().available_methods())
    #sampler = BackendSamplerV2(backend = AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    #pass_manager = generate_preset_pass_manager(0, AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    simulator = AerSimulator(device='GPU', method=method, cuStateVec_enable=True)

    #sampler = BackendSamplerV2(backend = AerSimulator())
    #pass_manager = generate_preset_pass_manager(0, AerSimulator())
    #simulator = AerSimulator(device='CPU', method=method)


    circ = transpile(qc, simulator)


    #exact_estimator = Estimator()
    # The circuit needs to be transpiled to the AerSimulator target

    #isa_circuit = pass_manager.run(qc)
    #pub = (isa_circuit, observable)
    #job = sampler.run([isa_circuit], shots=shots)
    job = simulator.run(circ, shots=shots)
    result = job.result()
    #print(f'backend: {result.backend_name}')
    #pub_result = result[0]
    #print("METADATA ", result)
    #counts = pub_result.data.meas.get_counts()
    counts = result.get_counts()
    print(counts)
    expectation_value = 0
    for key, value in counts.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value

@constraint(processors=[{"processorType": "CPU", "computingUnits": "1"},{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1)
def _expec_value_qiskit_gpu(qc, shots=1024, method='automatic'):
    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    #print(subcirc.draw())
    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    #print(obs_I)

    observables = "".join(observables)
    #print(observables)

    qc, _ = qc
    num_qubits = qc.num_qubits
    observable = SparsePauliOp(observables)
    #params = [0.1] * qc.num_parameters

    qc.measure_all()

    import subprocess
    import os


    try:
        subprocess.check_output('nvidia-smi')
        print('Nvidia GPU detected!')
    except Exception: # this command not being found can raise quite a few different errors depending on the configuration
        print('No Nvidia GPU in system!')


    target_gpus = [os.environ["COMPSS_BINDED_GPUS"]]
    print("TGPU ", target_gpus)
    print("CUDA ", [os.environ["CUDA_VISIBLE_DEVICES"]])
    print("DEVICES ", AerSimulator().available_devices())
    print("METHODS ", AerSimulator().available_methods())
    #sampler = BackendSamplerV2(backend = AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    #pass_manager = generate_preset_pass_manager(0, AerSimulator(method='statevector',device='GPU', cuStateVec_enable=True, batched_shots_gpu=True))
    simulator = AerSimulator(device='GPU', method=method, cuStateVec_enable=True)

    #sampler = BackendSamplerV2(backend = AerSimulator())
    #pass_manager = generate_preset_pass_manager(0, AerSimulator())
    #simulator = AerSimulator(device='CPU', method=method)


    circ = transpile(qc, simulator)


    #exact_estimator = Estimator()
    # The circuit needs to be transpiled to the AerSimulator target

    #isa_circuit = pass_manager.run(qc)
    #pub = (isa_circuit, observable)
    #job = sampler.run([isa_circuit], shots=shots)
    job = simulator.run(circ, shots=shots)
    result = job.result()
    #print(f'backend: {result.backend_name}')
    #pub_result = result[0]
    #print("METADATA ", result)
    #counts = pub_result.data.meas.get_counts()
    counts = result.get_counts()
    print(counts)
    expectation_value = 0
    for key, value in counts.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value

def _change_basis(circuit, observables):
    for idx,i in enumerate(observables):
        if i == "X":
            circuit.add(gates.H(idx))
        elif i == "Y":
            circuit.add(gates.SDG(idx))
            circuit.add(gates.H(idx))
        else:
            pass

    return circuit



def _check_ibm_qc(max_time_ibm):
    urls = {
      "http": "http://localhost:44433",
      "https": "http://localhost:44433"
    }

    proxies={"urls": urls}

    service = QiskitRuntimeService(channel="ibm_cloud", token=token, instance=instance, proxies=proxies)

    #backend = service.least_busy(operational=True, simulator=False)
    backend = service.backend("ibm_marrakesh")

    batch = Batch(backend=backend, max_time=max_time_ibm)
    print(batch)
    return batch, backend


def gate_cutting_subcircuits(
    dag: typing.Any,
    gates_cut: typing.List[typing.Any],
    software ="qiskit"
):
    if type(dag) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dag = circuit_qiskit_to_dag(dag)
    elif type(dag) == qibo.models.Circuit:
        dag = circuit_qibo_to_dag(dag)
    else:
        dag = dag

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in nx.connected_components(dag.to_undirected())
        ]
        subcircuits = []
        for s in S:
            tmp_cuts = []
            for c in gates_cut:
                if s.has_node(c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                subcirc = _execute_gate_cutting_subcircuits(dag, tmp_cuts,software=software)
                subcircuits.append(subcirc)
            else:
                s_new, highest_qubit = _update_qubits(s)
                if software == "qibo":
                    subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
                elif software=="qiskit":
                    subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
                else:
                    raise ValueError
                subcircuits.append(subcirc)
        return subcircuits
    else:
        if gates_cut:
            subcirc = _execute_gate_cutting_subcircuits(dag, gates_cut,software=software)
        else:
            s_new, highest_qubit = _update_qubits(dag)
            if software == "qibo":
                subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software=="qiskit":
                subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
        return subcirc

def _execute_gate_cutting_subcircuits(dag, gates_cut,software):
    new_dag = _generate_cut(dag,gates_cut)

    S = [new_dag.subgraph(c).copy() for c in nx.connected_components(new_dag.to_undirected())]

    subcircuits = []
    for index in range(6**len(gates_cut)):

        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(
            copy_graph.to_undirected()
        )

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())

        graph = _generate_gate_cutting(copy_graph, gates_cut,index, graph_components)


        for i, s in enumerate(graph_components):
            s_new, highest_qubit = _update_qubits(s)

            if software == "qibo":
                subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software=="qiskit":
                subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
            subcircuits.append(subcirc)
    return  subcircuits
