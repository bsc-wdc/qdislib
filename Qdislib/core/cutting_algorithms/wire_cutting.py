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

"""Wire cutting algorithms."""

import math
import networkx
import qibo
import qiskit
import typing
import time
import pickle
import os
import random

from qibo import models, gates
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import COLLECTION_IN
from pycompss.api.parameter import COLLECTION_OUT
from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
from Qdislib.utils.graph_qibo import dag_to_circuit_qibo, _dag_to_circuit_qibo_subcircuits
from Qdislib.utils.graph_qibo import _max_qubit
from Qdislib.utils.graph_qibo import _update_qubits
from Qdislib.utils.graph_qibo import _remove_red_edges
from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag,_dag_to_circuit_qiskit_subcircuits


def wire_cutting(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    sync: bool = True,
    observables: str = None,
    shots: int = 1024,
    backend: str = "numpy",
    qpu: str = None
):
    """
    Apply wire cutting to a quantum circuit to simplify its structure for distributed or hardware-constrained execution of subcircuits.

    Wire cutting removes specific connections (wires) between qubits, splitting the circuit into independent
    subcircuits. The function then reconstructs the expected value of the original circuit by combining the results
    from the subcircuits.

    Supports circuits created with Qiskit, Qibo, or as a DAG (NetworkX graph).

    Parameters
    ----------
    rand_qc : typing.Any
        Input quantum circuit to be cut. It can be a Qiskit QuantumCircuit, a Qibo Circuit, or a DAG (networkx.Graph).
    cut : typing.List[typing.Any]
        List of qubit wire cuts, specified as pairs of qubits (edges).
    sync : bool, default=True
        Whether to synchronize and wait for all circuit evaluations to complete (important in distributed settings like PyCOMPSs).
    observables : str, optional
        Observables to measure after applying the cuts. If provided, the circuit is first transformed into the measurement basis. By default, it will use Z observables.
    shots : int, default=1024
        Number of measurement shots for each subcircuit evaluation.
    backend : str, default="numpy"
        Backend to simulate or execute the circuits (e.g., "numpy").
    qpu : optional
        If specified (e.g., `"MN_Ona"`), indicates execution on a Quantum Processing Unit.

    Returns
    -------
    final_recons : float
        The reconstructed expectation value of the original quantum circuit after wire cutting.

    Notes
    -----
    - If the input circuit has multiple disconnected components, they are processed independently.
    - Each cut doubles the number of circuit evaluations required.
    - The final expectation value is scaled by \(1/2^{n_{\text{cuts}}}\) where \(n_{\text{cuts}}\) is the number of cuts.
    - GPU/QPU execution is supported if properly configured.
    - If no cuts are provided, the circuit is executed directly.

    Examples
    --------
    Applying wire cutting to a simple quantum circuit:

    .. code-block:: python

        from qiskit import QuantumCircuit
        from mymodule import wire_cutting

        qc = QuantumCircuit(2)
        qc.h(0)     # Gate "H_1"
        qc.cx(0, 1) # Gate "CZ_2"
        qc.measure_all()

        cuts = [("H_1","CZ_2")]  # List of gate tuples
        reconstruction = wire_cutting(qc, cuts, shots=2048, backend="qiskit")
    """
    if observables:
        rand_qc = _change_basis(rand_qc,observables)

    if type(rand_qc) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        if observables:
            if "I" in observables:
                dag = circuit_qiskit_to_dag(rand_qc, obs_I=observables)
            else:
                dag = circuit_qiskit_to_dag(rand_qc)
        else:
            dag = circuit_qiskit_to_dag(rand_qc)
    elif type(rand_qc) == qibo.models.Circuit:
        if observables:
            if "I" in observables:
                dag = circuit_qibo_to_dag(rand_qc,obs_I=observables)
            else:
                dag = circuit_qibo_to_dag(rand_qc)
        else:
            dag = circuit_qibo_to_dag(rand_qc)

    else:
        dag = rand_qc

    if networkx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected())
        ]
        results = []
        for s in S:
            num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in cut:
                if s.has_edge(*c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                graphs = _generate_wire_cutting(s, tmp_cuts, num_qubits=num_qubits, qpu=qpu)
                graphs = _sum_results(graphs)
                results.append(graphs)  # 1/(2**len(tmp_cuts))*sum(graphs)
            else:
                s_new, highest_qubit = _update_qubits(s)
                subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
                if qpu=="MN_Ona":
                    expected_value = _expec_value_qibo_qpu(subcirc, shots, method=backend)
                else:
                    expected_value = _execute_subcircuits(subcirc, shots, backend)
                results.append(expected_value)
        if sync:
            results = compss_wait_on(results)
        # Consider gate cutting within wire cutting:
        # if gate_cutting:
        #     final_recons = 1/2*sum(results)
        # else:
        final_recons = 1 / (2 ** len(cut)) * math.prod(results)
        return final_recons
    else:
        num_qubits = _max_qubit(dag)
        if cut:
            results = _generate_wire_cutting(dag, cut, num_qubits=num_qubits, qpu=qpu)

            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(cut)) * sum(results)
        else:
            s_new, highest_qubit = _update_qubits(dag)
            subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
            if qpu=="MN_Ona":
                final_recons = _expec_value_qibo_qpu(subcirc, shots, method=backend)
            else:
                final_recons = _execute_subcircuits(subcirc,shots,backend)
        return final_recons


def _generate_wire_cutting(
    dag: typing.Any,
    edges_to_replace: typing.List[typing.Tuple[str, str]],
    num_qubits: int,
    shots: int = 10000,
    backend: str = "numpy",
    qpu=None
) -> typing.List[float]:
    """Replace a specific edge in the DAG with a source and end node.

    :param dag: The directed acyclic graph (DAG) to modify.
    :param edge_to_replace: The edge to remove (tuple of nodes).
    :param num_qubits: The current number of qubits in the circuit.
    :return: The updated dag (the modified DAG with new source and end nodes).
    """
    reconstruction = []
    for index, edge_to_replace in enumerate(edges_to_replace, start=1):
        # Extract the nodes of the edge to be replaced
        source, target = edge_to_replace

        # Remove the original edge
        dag.remove_edge(source, target)

        source_gate_info = dag.nodes[source]

        target_gate_info = dag.nodes[target]

        common_qubit = list(
            set(target_gate_info.get("qubits")).intersection(
                set(source_gate_info.get("qubits"))
            )
        )

        successors = []
        # Iterate over all nodes in the graph
        for node in dag.nodes:
            if dag.has_edge(target, node):
                successors.append(node)

        # Include the target node itself
        nodes = [target] + successors

        for successor in nodes:
            qubits = dag.nodes[successor].get("qubits")
            for qubit in qubits:

                if common_qubit[0] is qubit:
                    temp_list = list(dag.nodes[successor].get("qubits"))

                    # Replace the common element with the new value
                    for i, element in enumerate(temp_list):
                        if element == common_qubit[0]:
                            temp_list[i] = num_qubits + index

                    updated_tuple = tuple(temp_list)
                    dag.nodes[successor]["qubits"] = updated_tuple

        dag.add_node(f"O_{index}", gate="S", qubits=common_qubit, parameters=())

        # Add the new end node with the same properties as the target node
        dag.add_node(
            f"PS_{index}", gate="T", qubits=(num_qubits + index,), parameters=()
        )

        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")

        copy_dag = dag.copy()
        red_edges = []
        for ed in dag.edges:
            if dag.get_edge_data(ed[0], ed[1])["color"] == "red":
                red_edges.append(ed)

        copy_dag.remove_edges_from(red_edges)

    graphs = []
    for i in range(8 ** len(edges_to_replace)):
        graphs.append(dag.copy())

    for index, graph in enumerate(graphs, start=0):
        copy_graph = graph.copy()

        copy_graph = _remove_red_edges(copy_graph)

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )

        graph_components = []
        for i in range(num_components):
            graph_components.append(networkx.DiGraph().copy())

        graph = _generate_subcircuits_wire_cutting(
            copy_graph,
            num_qubits + len(edges_to_replace),
            index,
            edges_to_replace,
            graph_components,
        )

        exp_value = []
        for s in graph_components:
            s_new, highest_qubit = _update_qubits(s)
            subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
            if qpu=="MN_Ona":
                expected_value = _expec_value_qibo_qpu(subcirc, shots, method=backend)
            else:
                expected_value = _execute_subcircuits(subcirc,shots, backend)
            exp_value.append(expected_value)
        exp_value = _change_sign(exp_value, index)
        reconstruction.append(exp_value)

    return reconstruction


@task(returns=1, graph_components=COLLECTION_OUT)
def _generate_subcircuits_wire_cutting(
    updated_dag: typing.Any,
    num_qubits: int,
    idx: int,
    edges_to_replace: typing.List[typing.Tuple[str, str]],
    graph_components: typing.List[typing.Any],
) -> typing.Any:
    """Generate the subcircuits from a circuit applying wire cutting.

    :param updated_dag: Given DAG.
    :param num_qubits: Current number of Qbits.
    :param idx: Index.
    :param edges_to_replace: List of edges to replace.
    :param graph_components: List of graph components.
    :raises ValueError: If it is not possible to generate subcircuits.
    :return: The updated DAG.
    """
    base8_rep = oct(idx)[2:]
    base8_rep = base8_rep.zfill(len(edges_to_replace))
    list_substitutions = list(map(int, base8_rep))

    for idx2, index in enumerate(list_substitutions, start=0):
        idx2 = idx2 + 1

        # I 0
        if index == 0:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "Observable I"
            # updated_dag.remove_node(f'O_{idx2}')
            updated_dag.remove_node(f"PS_{idx2}")

        # I 1
        elif index == 1:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "Observable I"
            # updated_dag.remove_node(f'O_{idx2}')
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "X"

        # X +
        elif index == 2:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "H"
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "H"

        # X -
        elif index == 3:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "H"
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "H"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="X",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")

        # Y +i
        elif index == 4:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "SDG"
            updated_dag.add_node(
                f"O2_{idx2}",
                gate="H",
                qubits=updated_dag.nodes[f"O_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"O_{idx2}", f"O2_{idx2}", color="blue")
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "S"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="H",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")

        # Y -i
        elif index == 5:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "SDG"
            updated_dag.add_node(
                f"O2_{idx2}",
                gate="H",
                qubits=updated_dag.nodes[f"O_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"O_{idx2}", f"O2_{idx2}", color="blue")
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "S"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="H",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")
            updated_dag.add_node(
                f"PS3_{idx2}",
                gate="X",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS3_{idx2}", f"PS2_{idx2}", color="blue")

        # Z 0
        elif index == 6:
            updated_dag.remove_node(f"O_{idx2}")
            updated_dag.remove_node(f"PS_{idx2}")

        # Z 1
        elif index == 7:
            updated_dag.remove_node(f"O_{idx2}")
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "X"

        else:
            Exception("Something went wrong preparing the combinations")

    updated_dag = _remove_red_edges(updated_dag)
    for i, c in enumerate(networkx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")

    return updated_dag


@task(returns=1, lst=COLLECTION_IN)
def _sum_results(lst: typing.List[int]) -> int:
    """Calculate the sum of all results.

    :param lst: List of partial results.
    :return: Sum of all partial results.
    """
    return sum(lst)


@task(returns=1)
def _execute_subcircuits(subcirc: typing.Any, shots=1024, backend="numpy") -> float:
    """Execute the given circuit.

    :param subcirc: Circuit to execute.
    :raises ValueError: If there is an unsupported observable.
    :return: The circuit expected value.
    """
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
   
    observables = "".join(observables)
    #print(observables)

    qibo.set_backend(backend)
    subcirc.add(gates.M(*range(subcirc.nqubits)))  #
    result = subcirc(nshots=shots)
    freq = dict(result.frequencies(binary=True))

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

    unique_code = str(time.time_ns())  # Nanosecond precision timestamp

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


@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign(expectation_value, index):
    """Get the product of all expected values and apply change of sign if necessary.

    :param expectation_value: List of expected values.
    :param index: Number to be used to determine if the sign has to be changed.
    :return: The expected value.
    """
    expectation_value = math.prod(expectation_value)
    number = index

    sign_change = False

    while number != 0:
        digit = number % 8  # Get the last digit
        if digit in {3, 5, 7}:  # Check if the digit is 3, 5, or 7
            sign_change = not sign_change  # Flip the sign change flag
        number //= 8  # Move to the next digit

    # If _change_sign is True, we flip the sign of the original number
    if sign_change:
        return -expectation_value
    else:
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


def wire_cutting_subcircuits(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    software = "qiskit"
):
    if type(rand_qc) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dag = circuit_qiskit_to_dag(rand_qc)
    elif type(rand_qc) == qibo.models.Circuit:
        dag = circuit_qibo_to_dag(rand_qc)
    else:
        dag = rand_qc

    if networkx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected())
        ]
        subcircuits = []
        for s in S:
            num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in cut:
                if s.has_edge(*c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                graphs = _generate_wire_cutting_subcircuits(s, tmp_cuts, num_qubits=num_qubits,software=software)
                subcircuits.append(graphs)
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
        num_qubits = _max_qubit(dag)
        if cut:
            subcircuits = _generate_wire_cutting_subcircuits(dag, cut, num_qubits=num_qubits,software=software)

        else:
            s_new, highest_qubit = _update_qubits(dag)
            if software == "qibo":
                subcircuits = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software=="qiskit":
                subcircuits = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
        return subcircuits

def _generate_wire_cutting_subcircuits(
    dag: typing.Any,
    edges_to_replace: typing.List[typing.Tuple[str, str]],
    num_qubits: int,
    software
) -> typing.List[float]:
    """Replace a specific edge in the DAG with a source and end node.

    :param dag: The directed acyclic graph (DAG) to modify.
    :param edge_to_replace: The edge to remove (tuple of nodes).
    :param num_qubits: The current number of qubits in the circuit.
    :return: The updated dag (the modified DAG with new source and end nodes).
    """
    reconstruction = []
    for index, edge_to_replace in enumerate(edges_to_replace, start=1):
        # Extract the nodes of the edge to be replaced
        source, target = edge_to_replace

        # Remove the original edge
        dag.remove_edge(source, target)

        source_gate_info = dag.nodes[source]

        target_gate_info = dag.nodes[target]

        common_qubit = list(
            set(target_gate_info.get("qubits")).intersection(
                set(source_gate_info.get("qubits"))
            )
        )

        successors = []
        # Iterate over all nodes in the graph
        for node in dag.nodes:
            if dag.has_edge(target, node):
                successors.append(node)

        # Include the target node itself
        nodes = [target] + successors

        for successor in nodes:
            qubits = dag.nodes[successor].get("qubits")
            for qubit in qubits:

                if common_qubit[0] is qubit:
                    temp_list = list(dag.nodes[successor].get("qubits"))

                    # Replace the common element with the new value
                    for i, element in enumerate(temp_list):
                        if element == common_qubit[0]:
                            temp_list[i] = num_qubits + index

                    updated_tuple = tuple(temp_list)
                    dag.nodes[successor]["qubits"] = updated_tuple

        dag.add_node(f"O_{index}", gate="S", qubits=common_qubit, parameters=())

        # Add the new end node with the same properties as the target node
        dag.add_node(
            f"PS_{index}", gate="T", qubits=(num_qubits + index,), parameters=()
        )

        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")

        copy_dag = dag.copy()
        red_edges = []
        for ed in dag.edges:
            if dag.get_edge_data(ed[0], ed[1])["color"] == "red":
                red_edges.append(ed)

        copy_dag.remove_edges_from(red_edges)

    graphs = []
    for i in range(8 ** len(edges_to_replace)):
        graphs.append(dag.copy())

    subcircuits = []
    for index, graph in enumerate(graphs, start=0):
        copy_graph = graph.copy()

        copy_graph = _remove_red_edges(copy_graph)

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )

        graph_components = []
        for i in range(num_components):
            graph_components.append(networkx.DiGraph().copy())

        graph = _generate_subcircuits_wire_cutting(
            copy_graph,
            num_qubits + len(edges_to_replace),
            index,
            edges_to_replace,
            graph_components,
        )

        
        for s in graph_components:
            s_new, highest_qubit = _update_qubits(s)
            if software == "qibo":
                subcirc = _dag_to_circuit_qibo_subcircuits(s_new, highest_qubit)
            elif software=="qiskit":
                subcirc = _dag_to_circuit_qiskit_subcircuits(s_new, highest_qubit)
            else:
                raise ValueError
            subcircuits.append(subcirc)
    return subcircuits
