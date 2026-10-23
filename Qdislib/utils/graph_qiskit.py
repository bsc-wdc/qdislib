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

from Qdislib.utils.graph_qibo import _final_measurements, _update_qubits_serie

from Qdislib.utils.runtime import task

import math


def _emit_qiskit_gate(circuit, gate_name, params, qubits):
    """Append one DAG gate to a Qiskit circuit. The DAG's canonical names are
    OpenQASM's, which are Qiskit's; other vocabularies (a DAG built by hand, the
    cut gadgets' ``CU3``) go through :func:`~Qdislib.utils.gates.canonical` first."""
    from Qdislib.utils.gates import canonical

    name, params = canonical(gate_name, list(params or []))
    getattr(circuit, name)(*params, *qubits)


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
    from Qdislib.utils.circuit_dag import ops_to_dag

    # A final measure (e.g. from measure_all) is left out: see _final_measurements.
    final = _final_measurements(
        [
            (
                gate.operation.name == "measure",
                (
                    []  # a barrier is no operation on its qubits
                    if gate.operation.name == "barrier"
                    else [circuit.find_bit(q).index for q in gate.qubits]
                ),
            )
            for gate in circuit.data
        ]
    )
    ops = [
        (
            f"{gate.operation.name}_{index}".upper(),
            gate.operation.name,
            # the global index, so qubits of different registers do not collide
            tuple(circuit.find_bit(q).index for q in gate.qubits),
            gate.operation.params,
        )
        for index, gate in enumerate(circuit.data, start=1)
        if gate.operation.name != "barrier" and index - 1 not in final
    ]
    return ops_to_dag(ops, obs_I, sdk="qiskit")


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
    cut_rounds, round_registers = {}, {1: creg_c}

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data["gate"]

        # Skip the observable I nodes (we'll handle them separately)
        if gate_name == "Observable I":
            continue

        qubits = node_data["qubits"]
        parameters = node_data["parameters"]

        if gate_name.lower() == "measure":
            # Gate-cut measure: a projective Z measurement that collapses the
            # qubit to its outcome |s> and lets it CONTINUE (the CZ cut's
            # (I +/- Z)/2 projector, Mitarai-Fujii arXiv:1909.07534). The outcome
            # is recorded in a cut-point register for the reconstruction sign;
            # the qubit is NOT reset. Resetting it (the old "force to |0>") threw
            # away the collapse and gave the wrong state to any gate after the cut.
            # Every cut keeps its own bit: a qubit's first cut goes to c[q], its
            # k-th to c{k}[q] (added on demand, so one cut per qubit -- the common
            # case -- still gives "<meas> <c>" keys). One shared c[q] kept only the
            # last of several cuts on a wire.
            (qubit,) = qubits
            k = cut_rounds[qubit] = cut_rounds.get(qubit, 0) + 1
            if k not in round_registers:
                round_registers[k] = ClassicalRegister(num_qubits, f"c{k}")
                circuit.add_register(round_registers[k])
            circuit.measure(qubit, round_registers[k][qubit])
        else:
            # CNOT->cx, U3->u, CU3->cu, ... (Qiskit renamed the OpenQASM-2 gates).
            _emit_qiskit_gate(circuit, gate_name, parameters, qubits)

    if obs_I:
        observables = ["Z"] * num_qubits
        for element in obs_I:
            observables[element] = "I"

        observables = "".join(observables)
        return [circuit, observables]

    observables = ["Z"] * num_qubits
    observables = "".join(observables)
    return [circuit, observables]


# --- SDK table adapters (see Qdislib.utils.sdks) ------------------------------


def qiskit_gates(circuit):
    """``[(label, qubits, is_measurement), ...]`` per instruction (a barrier acts
    on no qubit: it keeps its slot but is never a gate)."""
    return [
        (
            instruction.operation.name.upper(),
            (
                ()
                if instruction.operation.name == "barrier"
                else tuple(circuit.find_bit(q).index for q in instruction.qubits)
            ),
            instruction.operation.name == "measure",
        )
        for instruction in circuit.data
    ]


def qiskit_num_qubits(circuit):
    return circuit.num_qubits


def qiskit_statevector(circuit):
    """The exact state, qubit 0 the most significant bit."""
    import numpy as np
    from qiskit.quantum_info import Statevector

    n = circuit.num_qubits
    data = np.asarray(Statevector(circuit).data, dtype=complex)
    return data.reshape([2] * n).transpose(list(range(n))[::-1]).reshape(-1)
