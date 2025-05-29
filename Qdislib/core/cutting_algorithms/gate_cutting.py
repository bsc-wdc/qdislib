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


try:
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    from pycompss.api.parameter import COLLECTION_IN
    from pycompss.api.parameter import COLLECTION_OUT
    from pycompss.api.constraint import constraint
    from pycompss.api.implement import implement

    pycompss_available = True
except ImportError:
    print("PyCOMPSs is NOT installed. Proceeding with serial execution (no parallelism).")

    def task(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    def compss_wait_on(obj):
        return obj

    def constraint(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    def implement(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    COLLECTION_IN = COLLECTION_OUT = None
    pycompss_available = False


from Qdislib.utils.graph_qibo import _update_qubits
from Qdislib.utils.graph_qibo import _remove_red_edges
from Qdislib.utils.graph_qibo import (
    _update_qubits,
    circuit_qibo_to_dag,
    plot_dag,
    _max_qubits_graph,
)
from Qdislib.utils.graph_qiskit import (
    dag_to_circuit_qiskit,
    circuit_qiskit_to_dag,
    _dag_to_circuit_qiskit_subcircuits,
)
from Qdislib.utils.graph_qibo import (
    dag_to_circuit_qibo,
    circuit_qibo_to_dag,
    _dag_to_circuit_qibo_subcircuits,
)
from Qdislib.core.cutting_algorithms.wire_cutting import _sum_results
from Qdislib.core.find_cut.find_cut import _find_nodes_with_qubit

from qiskit_ibm_runtime import (
    QiskitRuntimeService,
    Batch,
    SamplerV2 as Sampler,
    EstimatorV2 as Estimator,
)


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
    verbose: bool = False
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
    verbose : bool, default=False
        Receive information about middle steps

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
        batch, backend = _check_ibm_qc()
    else:
        batch, backend = None, None

    if observables:
        dag = _change_basis(dag, observables)

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
                dag = circuit_qibo_to_dag(dag, obs_I=observables)
            else:
                dag = circuit_qibo_to_dag(dag)
        else:
            dag = circuit_qibo_to_dag(dag)
    else:
        dag = dag

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy() for c in nx.connected_components(dag.to_undirected())
        ]
        results = []
        for s in S:
            # num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in gates_cut:
                if s.has_node(c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                graphs = _execute_gate_cutting(
                    dag,
                    tmp_cuts,
                    shots=shots,
                    method=method,
                    gpu=gpu,
                    gpu_min_qubits=gpu_min_qubits,
                    qpu=qpu,
                    qpu_dict=qpu_dict,
                    batch=batch,
                    backend=backend,
                    verbose=verbose
                )
                graphs = _sum_results(graphs)
                results.append(graphs)
            else:
                _max_qubit = _max_qubits_graph(s)
                s_new, highest_qubit = _update_qubits(s)
                subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
                if qpu and "MN_Ona" in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                    expected_value = _expec_value_qibo_qpu(
                        subcirc, shots=shots, method=method
                    )
                elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                    expected_value = _expec_value_qiskit_gpu(
                        subcirc, shots=shots, method=method
                    )
                elif (
                    qpu
                    and "IBM_Quantum" in qpu_dict
                    and qpu_dict["IBM_Quantum"] >= _max_qubit
                ):
                    expected_value = _expec_value_qiskit_qpu(
                        subcirc, shots=shots, batch=batch, backend=backend
                    )
                else:
                    expected_value = _expec_value_qiskit(
                        subcirc, shots=shots, method=method
                    )
                results.append(expected_value)
        if sync:
            results = compss_wait_on(results)          
        final_recons = 1 / (2 ** len(gates_cut)) * math.prod(results)
        if verbose:
            print("Multiply top and bottom expected value of subcircuits")
            for idx, el in enumerate(results,start=1):
                print(f"Subcircuit {idx} expected value: {el}")
            print(f"Reconstructed original circuit expected value: {final_recons}")
        return final_recons
    else:
        if gates_cut:
            results = _execute_gate_cutting(
                dag,
                gates_cut,
                shots=shots,
                method=method,
                gpu=gpu,
                gpu_min_qubits=gpu_min_qubits,
                qpu=qpu,
                qpu_dict=qpu_dict,
                batch=batch,
                backend=backend,
                verbose=verbose
            )
            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(gates_cut)) * sum(results)
            if verbose:
                print("Multiply top and bottom expected value of subcircuits")
                for idx, el in enumerate(results,start=1):
                    print(f"Subcircuit {idx} expected value: {el}")
                print(f"Reconstructed original circuit expected value: {final_recons}")
        else:
            _max_qubit = _max_qubits_graph(dag)
            s_new, highest_qubit = _update_qubits(dag)
            subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
            if qpu and "MN_Ona" in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                final_recons = _expec_value_qibo_qpu(
                    subcirc, shots=shots, method=method
                )
            elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                final_recons = _expec_value_qiskit_gpu(
                    subcirc, shots=shots, method=method
                )
            elif (
                qpu
                and "IBM_Quantum" in qpu_dict
                and qpu_dict["IBM_Quantum"] >= _max_qubit
            ):
                final_recons = _expec_value_qiskit_qpu(
                    subcirc, shots=shots, batch=batch, backend=backend
                )
            else:
                final_recons = _expec_value_qiskit(subcirc, shots, method=method)
            if sync:
                final_recons = compss_wait_on(final_recons)
            if verbose:
                print(f"No CUT performed, expected value circuit: {final_recons}")
        return final_recons


def _generate_cut(dag, gates_cut):
    dag_copy = dag.copy()
    dag_copy = _remove_red_edges(dag_copy)
    for index, gate_name in enumerate(gates_cut, start=1):
        target_qubits = dag_copy.nodes[gate_name]["qubits"]
        pred_0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[0], direction="predecessor"
        )

        pred_1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[1], direction="predecessor"
        )

        succ_0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[0], direction="successor"
        )

        succ_1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[1], direction="successor"
        )

        dag_copy.remove_node(gate_name)

        dag_copy.add_node(
            f"SUBS1_{index}", gate="S", qubits=(target_qubits[0],), parameters=()
        )
        dag_copy.add_node(
            f"SUBS2_{index}", gate="S", qubits=(target_qubits[1],), parameters=()
        )

        if pred_0:
            dag_copy.add_edge(pred_0[0], f"SUBS1_{index}", color="blue")

        if succ_0:
            dag_copy.add_edge(f"SUBS1_{index}", succ_0[0], color="blue")

        if pred_1:
            dag_copy.add_edge(pred_1[0], f"SUBS2_{index}", color="blue")

        if succ_1:
            dag_copy.add_edge(f"SUBS2_{index}", succ_1[0], color="blue")

    return dag_copy


def _decimal_to_base6(num):
    if num == 0:
        return "0"

    base6 = ""
    while num > 0:
        base6 = str(num % 6) + base6
        num //= 6
    return base6


# @constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, graph_components=COLLECTION_OUT)
def _generate_gate_cutting(updated_dag, gates_cut, index, graph_components):

    base6_rep = _decimal_to_base6(index)
    base6_rep = base6_rep.zfill(len(gates_cut))
    list_substitutions = list(map(int, base6_rep))

    # print(list_substitutions)
    for idx2, idx in enumerate(list_substitutions, start=1):
        list_succ1 = list(updated_dag.succ[f"SUBS1_{idx2}"])
        list_pred1 = list(updated_dag.pred[f"SUBS1_{idx2}"])

        list_succ2 = list(updated_dag.succ[f"SUBS2_{idx2}"])
        list_pred2 = list(updated_dag.pred[f"SUBS2_{idx2}"])

        # 1 - Rz(-pi/2) -- Rz(-pi/2)
        if idx == 0:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS1_{idx2}"]["parameters"] = [-np.pi / 2]
            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS2_{idx2}"]["parameters"] = [-np.pi / 2]

        # 2 - Z Rz(-pi/2) -- Z Rz(-pi/2)
        elif idx == 1:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "z"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "z"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        # 3 - MEASURE Rz(-pi/2) -- Rz(-pi)
        elif idx == 2:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS2_{idx2}"]["parameters"] = [-np.pi]

        # 4 - MEASURE Rz(-pi/2) -- RES
        elif idx == 3:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.remove_node(f"SUBS2_{idx2}")

            if list_succ2 and list_pred2:
                succ2 = list_succ2[0]
                pred2 = list_pred2[0]

                updated_dag.add_edge(pred2, succ2)

        # 5 - Rz(-pi) -- MEASURE Rz(-pi/2)
        elif idx == 4:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS1_{idx2}"]["parameters"] = [-np.pi]

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        # 6 - RES -- MEASURE Rz(-pi/2)
        elif idx == 5:
            updated_dag.remove_node(f"SUBS1_{idx2}")

            if list_succ1 and list_pred1:
                succ1 = list_succ1[0]
                pred1 = list_pred1[0]
                updated_dag.add_edge(pred1, succ1, color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        else:
            raise TypeError

    for i, c in enumerate(nx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")

    return updated_dag


# @constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign_gate_cutting(expectation_value, index):
    expectation_value = [x for x in expectation_value if x is not None]
    expectation_value = math.prod(expectation_value)
    number = index

    change_sign = False

    while number != 0:
        digit = number % 6  # Get the last digit
        if digit in {3, 5}:  # Check if the digit is 3, 5
            change_sign = not change_sign  # Flip the sign change flag
        number //= 6  # Move to the next digit

    # If change_sign is True, we flip the sign of the original number
    if change_sign:
        return -expectation_value
    else:
        return expectation_value


def _execute_gate_cutting(
    dag,
    gates_cut,
    shots=10000,
    method="automatic",
    gpu=False,
    gpu_min_qubits=None,
    qpu=None,
    qpu_dict=None,
    batch=None,
    backend=None,
    verbose=False
):
    new_dag = _generate_cut(dag, gates_cut)

    S = [
        new_dag.subgraph(c).copy()
        for c in nx.connected_components(new_dag.to_undirected())
    ]

    max_qubit_list = []
    for s in S:
        _max_qubit = _max_qubits_graph(s)
        max_qubit_list.append(_max_qubit)
    
    bool_mult_qpu = True

    reconstruction = []
    for index in range(6 ** len(gates_cut)):
        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(copy_graph.to_undirected())

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())

        graph = _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        exp_values = []
        for i, s in enumerate(graph_components):
            _max_qubit = max_qubit_list[i]
            s_new, highest_qubit = _update_qubits(s)
            if verbose:
                print(f"SUBCIRCUIT {index*2+i}")
            if qpu and "MN_Ona" in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
            else:
                subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
            if qpu and "MN_Ona" in qpu_dict and qpu_dict["MN_Ona"] >= _max_qubit:
                exp = _expec_value_qibo_qpu(subcirc, shots=shots, method=method)
            elif gpu and gpu_min_qubits <= _max_qubit and _max_qubit <= 30:
                exp = _expec_value_qiskit_gpu(subcirc, shots=shots, method=method)
            elif (
                qpu
                and "IBM_Quantum" in qpu_dict
                and qpu_dict["IBM_Quantum"] >= _max_qubit
            ):
                exp = _expec_value_qiskit_qpu(
                    subcirc, shots=shots, batch=batch, backend=backend
                )
            else:
                exp = _expec_value_qiskit(subcirc, shots=shots, method=method, verbose=verbose)
            exp_values.append(exp)

        exp = _change_sign_gate_cutting(exp_values, index)
        reconstruction.append(exp)

        if bool_mult_qpu:
            bool_mult_qpu = False
        else:
            bool_mult_qpu = True

    return reconstruction


@task(returns=1)
def _expec_value_qibo_qpu(subcirc, shots=1024, method="numpy"):

    if subcirc is None:
        return None

    tmp = subcirc[1]
    subcirc = subcirc[0]

    if tmp:
        obs_I = tmp
    else:
        obs_I = None
    observables = ["Z"] * subcirc.nqubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    subcirc.add(gates.M(*range(subcirc.nqubits)))

    observables = "".join(observables)

    unique_code = str(time.time_ns())

    circuit_filename = (
        f"/home/bsc/bsc019635/ona_proves/subcircuits/circuit_{unique_code}.pkl"
    )
    result_filename = (
        f"/home/bsc/bsc019635/ona_proves/subcircuits/result_{unique_code}.pkl"
    )

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

    os.remove(result_filename)  # Remove result file after reading

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
    return expectation_value


@constraint(processors=[{"processorType": "CPU", "computingUnits": "1"}])
@task(returns=1)
def _expec_value_qiskit(qc, shots=1024, method="automatic",verbose=False):

    if qc is None:
        return None

    if type(qc) != qiskit.circuit.quantumcircuit.QuantumCircuit:
        tmp = qc[1]
        subcirc = qc[0]
    else:
        tmp = False
        subcirc = qc

    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    observables = "".join(observables)
    qc = subcirc
    qc.measure_all()

    import subprocess
    import os

    try:
        subprocess.check_output("nvidia-smi")
        if verbose:
            print("Nvidia GPU detected!")
    except:
        if verbose:
            print("No Nvidia GPU in system!")

    if verbose:
        print("Simulation using CPUs")
        print("Method: ",method)
    simulator = AerSimulator(
        device="CPU",
        method=method,
        max_parallel_threads=1,
        mps_omp_threads=1,
        mps_parallel_threshold=1,
    )

    job = simulator.run(qc, shots=shots)
    result = job.result()
    if verbose:
        print(f"Backend: {result.backend_name}")

    counts = result.get_counts()
    if verbose:
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
    if verbose:
        print("Expected Value: ", expectation_value)
        print("\n")
    return expectation_value


@task(returns=1)
def _expec_value_qiskit_qpu(
    qc, shots=1024, method="automatic", batch=None, backend=None
):

    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    observables = "".join(observables)

    qc, _ = qc
    qc.measure_all()

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)

    sampler = Sampler(mode=batch)
    job = sampler.run([isa_circuit])
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
    # print(expectation_value)
    return expectation_value


# @implement(source_class="Qdislib.core.cutting_algorithms.gate_cutting", method="_expec_value_qiskit")
@constraint(
    processors=[
        {"processorType": "CPU", "computingUnits": "1"},
        {"processorType": "GPU", "computingUnits": "1"},
    ]
)
@task(returns=1)
def _expec_value_qiskit_gpu_implements(qc, shots=1024, method="automatic"):
    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    observables = "".join(observables)

    qc, _ = qc
    qc.measure_all()

    import subprocess
    import os

    try:
        subprocess.check_output("nvidia-smi")
        print("Nvidia GPU detected!")
    except (
        Exception
    ):
        print("No Nvidia GPU in system!")

    simulator = AerSimulator(device="GPU", method=method, cuStateVec_enable=True)

    circ = transpile(qc, simulator)

    job = simulator.run(circ, shots=shots)
    result = job.result()

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
    return expectation_value


@constraint(
    processors=[
        {"processorType": "CPU", "computingUnits": "1"},
        {"processorType": "GPU", "computingUnits": "1"},
    ]
)
@task(returns=1)
def _expec_value_qiskit_gpu(qc, shots=1024, method="automatic"):
    if qc is None:
        return None

    tmp = qc[1]
    subcirc = qc[0]

    if tmp:
        obs_I = tmp
    else:
        obs_I = None

    observables = ["Z"] * subcirc.num_qubits

    if obs_I:
        for element in obs_I:
            observables[element] = "I"

    observables = "".join(observables)
    qc, _ = qc
    qc.measure_all()

    import subprocess
    import os

    try:
        subprocess.check_output("nvidia-smi")
        print("Nvidia GPU detected!")
    except (
        Exception
    ):
        print("No Nvidia GPU in system!")

    simulator = AerSimulator(device="GPU", method=method, cuStateVec_enable=True)

    circ = transpile(qc, simulator)

    job = simulator.run(circ, shots=shots)
    result = job.result()
    
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
    return expectation_value


def _change_basis(circuit, observables):
    for idx, i in enumerate(observables):
        if i == "X":
            circuit.add(gates.H(idx))
        elif i == "Y":
            circuit.add(gates.SDG(idx))
            circuit.add(gates.H(idx))
        else:
            pass

    return circuit


def _check_ibm_qc():
    urls = {"http": "http://localhost:44433", "https": "http://localhost:44433"}

    proxies = {"urls": urls}

    token = os.environ.get("IBM_QUANTUM_TOKEN", "")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE","")
    channel = os.environ.get("IBM_QUANTUM_CHANNEL", "")
    max_time_ibm = os.environ.get("IBM_QUANTUM_MAX_TIME", 60)
    quantum_chip = os.environ.get("IBM_QUANTUM_QPU_NAME", "")

    print(quantum_chip)
    
    if channel == "ibm_quantum":
        service = QiskitRuntimeService(
            channel=channel, token=token, proxies=proxies
        )
    elif channel == "ibm_cloud":
        service = QiskitRuntimeService(
            channel=channel,token=token, instance=instance #, proxies=proxies
        )
    else:
        raise ValueError

    if quantum_chip != "":
        backend = service.backend(quantum_chip)
    else:
        backend = service.least_busy(operational=True, simulator=False)

    print(backend.operation_names)

    batch = Batch(backend=backend, max_time=int(max_time_ibm))
    print(batch)
    return batch, backend


def gate_cutting_subcircuits(
    dag: typing.Any, gates_cut: typing.List[typing.Any], software="qiskit"
): 
    """
    Partition a quantum circuit into subcircuits via gate cutting.

    This function takes a quantum circuit (or its DAG representation) and a list of gates to cut.
    It returns a set of subcircuits resulting from the cutting process. The function supports
    circuits defined using either Qiskit or Qibo.

    Parameters
    ----------
    dag : Union[qiskit.QuantumCircuit, qibo.models.Circuit, nx.Graph]
        The input quantum circuit. It can be a Qiskit `QuantumCircuit`, a Qibo `Circuit`,
        or a preprocessed DAG/networkx graph representation.

    gates_cut : List[Any]
        A list of gate identifiers (typically node references) to be used as cutting points
        in the circuit. These gates are where the circuit will be split into subcircuits.

    software : str, optional
        The target software framework for the subcircuit output. Options are:
        - `"qiskit"` (default): returns subcircuits as Qiskit `QuantumCircuit` objects
        - `"qibo"`: returns subcircuits as Qibo `Circuit` objects

    Returns
    -------
    Union[List[Any], Any]
        A list of subcircuits needed to recostruct the original expected value. The type of each subcircuit depends on the `software` argument.

    Raises
    ------
    ValueError
        If an unsupported `software` value is provided.

    Notes
    -----
    - If the input circuit contains multiple connected components, each component is processed separately.
    - If no cuts are specified, the entire circuit is returned as one subcircuit after qubit reindexing.
    - Internally, Qiskit and Qibo circuits are converted to DAGs using appropriate helper functions.
    - This function relies on internal utilities such as `_execute_gate_cutting_subcircuits`, `_update_qubits`,
      `_dag_to_circuit_qiskit_subcircuits`, and `_dag_to_circuit_qibo_subcircuits`.

    Examples
    --------
    >>> qc = QuantumCircuit(2)
    >>> qc.cx(0, 1)
    >>> subcircuits = gate_cutting_subcircuits(qc, gates_cut=["CX_1"])
    >>> for subcircuit in subcircuits:
    >>>     print(subcircuit)
    """
    if type(dag) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dag = circuit_qiskit_to_dag(dag)
    elif type(dag) == qibo.models.Circuit:
        dag = circuit_qibo_to_dag(dag)
    else:
        dag = dag

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy() for c in nx.connected_components(dag.to_undirected())
        ]
        subcircuits = []
        for s in S:
            tmp_cuts = []
            for c in gates_cut:
                if s.has_node(c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                subcirc = _execute_gate_cutting_subcircuits(
                    dag, tmp_cuts, software=software
                )
                subcircuits.append(subcirc)
            else:
                s_new, highest_qubit = _update_qubits(s)
                if software == "qibo":
                    subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
                elif software == "qiskit":
                    subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
                else:
                    raise ValueError
                subcircuits.append(subcirc)
        return subcircuits
    else:
        if gates_cut:
            subcirc = _execute_gate_cutting_subcircuits(
                dag, gates_cut, software=software
            )
        else:
            s_new, highest_qubit = _update_qubits(dag)
            if software == "qibo":
                subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software == "qiskit":
                subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
        return subcirc


def _execute_gate_cutting_subcircuits(dag, gates_cut, software):
    new_dag = _generate_cut(dag, gates_cut)

    S = [
        new_dag.subgraph(c).copy()
        for c in nx.connected_components(new_dag.to_undirected())
    ]

    subcircuits = []
    for index in range(6 ** len(gates_cut)):

        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(copy_graph.to_undirected())

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())

        graph = _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        for i, s in enumerate(graph_components):
            s_new, highest_qubit = _update_qubits(s)

            if software == "qibo":
                subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software == "qiskit":
                subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
            subcircuits.append(subcirc)
    return subcircuits

def gate_cutting_subcircuit_reconstruction(expected_values, number_cuts, number_components=2):
    """
    Reconstruct the final expectation value from subcircuit results after gate cutting.

    This function performs post-processing of expectation values obtained from
    individually simulated or measured subcircuits. It applies the gate-cutting
    sign correction and computes the final observable estimate according to the
    number of cuts made in the original quantum circuit.

    Parameters
    ----------
    expected_values : List[float]
        A flat list of raw expectation values resulting from the execution
        of subcircuits. It is expected to be of even length, as values are
        grouped in pairs corresponding to measurement outcomes.

    number_cuts : int
        The total number of gate cuts performed in the original circuit.
        This determines the normalization factor used in the final result.
    
    number_components : int
        The number of independent parts the initial subcircuits has.

    Returns
    -------
    float
        The final reconstructed expectation value of the original circuit,
        computed as the sum of sign-adjusted subcircuit values divided by
        2 raised to the power of the number of cuts.

    Notes
    -----
    - The `expected_values` list is divided into pairs, and each pair is processed
      using the `_change_sign_gate_cutting` function, which applies sign corrections
      based on the cut structure.
    - The computation uses `compss_wait_on()` to ensure all results are synchronized
      before aggregation. This is relevant in distributed or parallel execution contexts (using PyCOMPSs).
    - The final result scales by a factor of `1 / (2 ** number_cuts)`, consistent with
      the probabilistic interpretation of cut circuits.

    Examples
    --------
    >>> raw_values = [0.5, -0.5, 0.3, -0.3, 0.5, -0.5, 0.3, -0.3, 0.5, -0.5, 0.3, -0.3]
    >>> number_of_cuts = 1
    >>> reconstructed_value = gate_cutting_subcircuit_reconstruction(raw_values, number_of_cuts, number_components=2)
    >>> print(reconstructed_value)
    0.0
    """
    if number_cuts == 0:
        return math.prod(expected_values)
    else:
        expected_values = [expected_values[i:i+number_components] for i in range(0, len(expected_values), number_components)]
        results = []
        for index, exp_values in enumerate(expected_values):
            result = _change_sign_gate_cutting(exp_values, index)
            results.append(result)
        
        results = compss_wait_on(results)
        final_recons = 1 / (2 ** number_cuts) * sum(results)
        return final_recons

