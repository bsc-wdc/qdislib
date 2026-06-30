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

"""
Qdislib qiskit graph utils.

This file contains all auxiliary qiskit graph classes and functions.
"""

from __future__ import annotations  # defer annotation eval so qiskit can be optional

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_qiskit_to_dag",
    "dag_to_circuit_qiskit",
]

import networkx

try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister

    qiskit_available = True
except ImportError:  # qiskit is optional; functions that need it fail when called
    QuantumCircuit = QuantumRegister = ClassicalRegister = None
    qiskit_available = False

from Qdislib.utils.graph_qibo import _update_qubits_serie

from Qdislib.utils.runtime import task


def circuit_qiskit_to_dag(circuit: QuantumCircuit, obs_I=None) -> networkx.DiGraph:
    """Convert a Qiskit quantum circuit into a DAG where each node stores gate information.

    This function takes a :class:`qiskit.QuantumCircuit` object and builds a directed acyclic graph (DAG)
    representation, where each node corresponds to a gate in the circuit and edges capture qubit dependencies
    between gates. Each node stores metadata including gate type, acting qubits, and parameters.

    Optional post-processing allows for the insertion of "Observable I" nodes at sink positions in the graph,
    useful for modeling identity observables in quantum measurement pipelines.

    :param circuit: A :class:`qiskit.QuantumCircuit` object representing the circuit to convert.
    :param obs_I: (Optional) A list of characters (e.g., `['X', 'I', 'Z', ...]`) where each entry corresponds
                  to a qubit, and an `"I"` indicates the insertion of an "Observable I" node at the end of
                  the DAG path for that qubit. Internal use for the workflow when cutting.

    :return: A :class:`networkx.DiGraph` representing the circuit DAG. Each node in the DAG has the following attributes:
             - `gate`: The name of the gate (string).
             - `qubits`: A tuple of qubit indices the gate acts upon.
             - `parameters`: A list of parameters passed to the gate.

    Notes:
        - Barrier operations are ignored and not included in the DAG.
        - Edges are added between nodes that act on shared qubits, in the order of execution.
          - Blue edges represent direct dependency (first match).
          - Red edges represent additional constraints due to qubit reuse without direct lineage.
        - Observable I nodes, if added, are terminal nodes connected to the corresponding final gate on that qubit.

    Example:
        .. code-block:: python

            from qiskit import QuantumCircuit
            import Qdislib as qd

            qc = QuantumCircuit(2)
            qc.h(0)
            qc.cx(0, 1)

            dag = qd.circuit_qiskit_to_dag(qc)
            print(dag.nodes(data=True))
    """
    # Create a directed graph
    dag = networkx.DiGraph()

    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.data, start=1):

        if gate.operation.name != "barrier":

            # Unique identifier for each gate instance
            gate_name = f"{gate.operation.name}_{gate_idx}".upper()

            qubits = ()
            for elem in gate.qubits:
                qubits = qubits + (elem._index,)

            # Add the gate to the DAG, including the gate type, qubits, and parameters
            dag.add_node(
                gate_name,
                gate=gate.operation.name,
                qubits=qubits,
                parameters=gate.operation.params,
            )

            # Connect gates based on qubit dependencies
            for qubit in qubits:
                # Skip the last node since it is the current gate being added
                for pred_gate in reversed(list(dag.nodes)):
                    # Check if the qubit is in the node's qubits
                    if (
                        dag.nodes[pred_gate].get("qubits")
                        and qubit in dag.nodes[pred_gate]["qubits"]
                        and gate_name != pred_gate
                    ):
                        dag.add_edge(pred_gate, gate_name, color="blue")
                        break

            for qubit in qubits:
                for pred_gate in reversed(list(dag.nodes)):
                    if (
                        dag.nodes[pred_gate].get("qubits")
                        and qubit in dag.nodes[pred_gate]["qubits"]
                        and gate_name != pred_gate
                        and not dag.has_edge(pred_gate, gate_name)
                    ):
                        dag.add_edge(pred_gate, gate_name, color="red")
    if obs_I:
        sinks = [node for node in dag.nodes if dag.out_degree(node) == 0]
        for idx, i in enumerate(obs_I):
            if i == "I":
                for elem in sinks:
                    if dag.nodes[elem].get("qubits") == (idx,):
                        dag.add_node(
                            f"OBSI_{idx}",
                            gate="Observable I",
                            qubits=(idx,),
                            parameters=(),
                        )
                        dag.add_edge(elem, f"OBSI_{idx}", color="blue")
    return dag


@task(returns=1)
def dag_to_circuit_qiskit(dag, num_qubits):
    """Reconstruct a Qiskit quantum circuit from a DAG representation.

    This function transforms a directed acyclic graph (DAG), where each node represents a quantum gate
    and contains metadata (e.g., gate type, acting qubits, parameters), back into a Qiskit
    :class:`qiskit.QuantumCircuit` object. The function handles standard gates, controlled gates,
    and optional Observable-I nodes used for measurement handling.

    Args:
        dag (networkx.DiGraph): The DAG representing a quantum circuit. Each node must include the following attributes:
            - 'gate': (str) Name of the gate.
            - 'qubits': (tuple) Qubit indices the gate acts upon.
            - 'parameters': (list or None) Parameters for parametric gates (e.g., rotations).

        num_qubits (int): Number of qubits in the original circuit.

    Returns:
        list: A list containing:
            - QuantumCircuit: The reconstructed Qiskit circuit object.
            - list or None: A list of qubit indices measured with an "Observable I" (identity) gate, or `None` if none are present.

    Notes:
        - Observable I nodes are removed from the DAG and recorded separately as classical measurement placeholders.
        - The function uses a topological sort of the DAG to ensure gate application respects dependencies.
        - The controlled-X gate is accepted under either name (Qibo's ``CNOT`` or Qiskit/CUDA-Q's ``CX``) and emitted as Qiskit's ``cx``.
        - If a `measure` gate is encountered, classical control flow (via `c_if`) is assumed.
        - The function assumes all remaining qubits are measured at the end unless otherwise specified.

    Example:
        .. code-block:: python

            import Qdislib as qd

            dag = qd.circuit_qiskit_to_dag(qiskit_circuit)
            reconstructed_circuit, _ = qd.dag_to_circuit_qiskit(dag, num_qubits=2)
            print(reconstructed_circuit.draw())
    """
    if dag is None:
        return None

    topo_order = list(networkx.topological_sort(dag))

    # Optionally handle measurements, assuming all qubits are measured at the end
    obs_I = []
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data["gate"] == "Observable I":
            # print(node)
            obs_I.append(node_data["qubits"][0])
            dag.remove_node(node)

    dag, highest_qubit, smalles_qubit = _update_qubits_serie(dag)

    # Create an empty Qibo circuit
    qreg_q = QuantumRegister(num_qubits, "q")
    creg_c = ClassicalRegister(num_qubits, "c")
    circuit = QuantumCircuit(qreg_q, creg_c)

    # Traverse the DAG in topological order
    topo_order = list(networkx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data["gate"]

        # Skip the observable I nodes (we'll handle them separately)
        if gate_name == "Observable I":
            continue

        # Two-qubit controlled-X: Qibo names it CNOT, Qiskit/CUDA-Q name it CX.
        # Emit Qiskit's native cx so all three backends keep this gate as-is
        # (reconstruction runs after cutting, so no CZ decomposition is needed).
        elif gate_name.upper() in ("CNOT", "CX"):
            circuit.cx(*node_data["qubits"])
            continue

        # Get the qubits this gate acts on

        qubits = node_data["qubits"]
        parameters = node_data["parameters"]

        # Qiskit method name for this gate (e.g. "cz", "rz", "measure").
        gate_method = gate_name.lower()

        if parameters:
            getattr(circuit, gate_method)(*parameters, *qubits)
        elif gate_method == "measure":
            # Gate-cut measure: measure, then force the qubit to |1> (matches the
            # Qibo/CUDA-Q gate-cut measure convention).
            clbit = qubits
            getattr(circuit, gate_method)(*qubits, *clbit)
            with circuit.if_test((*qubits, *clbit)):
                circuit.x(*qubits)
        else:
            getattr(circuit, gate_method)(*qubits)

    if obs_I:
        observables = ["Z"] * num_qubits
        for element in obs_I:
            observables[element] = "I"

        observables = "".join(observables)
        return [circuit, observables]

    observables = ["Z"] * num_qubits
    observables = "".join(observables)
    return [circuit, observables]
