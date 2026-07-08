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
Qdislib graph utils.

This file contains all auxiliary graph classes and functions.
"""

from __future__ import annotations  # defer annotation eval so qibo can be optional

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_qibo_to_dag",
    "dag_to_circuit_qibo",
    "plot_dag",
]

import networkx
import matplotlib.pyplot as plt
import numpy as np
import typing

try:
    from qibo import models, gates

    qibo_available = True
except ImportError:  # qibo is optional; functions that need it fail when called
    models = gates = None
    qibo_available = False

from Qdislib.utils.runtime import task

# Cache for the case-insensitive gate lookup (built lazily on first miss).
_QIBO_GATE_BY_UPPER = None


def _resolve_qibo_gate(name: str):
    """Resolve a gate class name to its ``qibo.gates`` class, case-insensitively.

    DAG nodes store the gate as ``gate.__class__.__name__`` (correct case, e.g.
    ``"iSWAP"``, ``"fSim"``), but different call sites upper/lower-case it before a
    lookup. qibo's class names are case-sensitive, so a plain
    ``getattr(gates, name.upper())`` silently misses every mixed-case gate. Try the
    exact attribute first (the common path), then fall back to a case-insensitive
    match so ``"ISWAP"``/``"iswap"`` both resolve to :class:`qibo.gates.iSWAP`.

    :raises AttributeError: if no gate class matches ``name`` in any case.
    """
    cls = getattr(gates, name, None)
    if isinstance(cls, type):
        return cls
    global _QIBO_GATE_BY_UPPER
    if _QIBO_GATE_BY_UPPER is None:
        _QIBO_GATE_BY_UPPER = {
            attr.upper(): obj
            for attr in dir(gates)
            if isinstance((obj := getattr(gates, attr, None)), type)
        }
    try:
        return _QIBO_GATE_BY_UPPER[name.upper()]
    except KeyError:
        raise AttributeError(f"qibo.gates has no gate class {name!r}")


def _unpack_parameters(parameters):
    """Normalise a DAG node's ``parameters`` into a flat tuple to splat into a gate
    constructor. qibo returns a tuple, but guard against a bare scalar or ``None`` so
    multi-parameter gates (``U2``, ``U3``, ``fSim``) get each angle as its own
    positional argument instead of one packed tuple."""
    if parameters is None:
        return ()
    if isinstance(parameters, (tuple, list, np.ndarray)):
        return tuple(parameters)
    return (parameters,)


def circuit_qibo_to_dag(circuit: models.Circuit, obs_I=None) -> networkx.DiGraph:
    """Convert a Qibo circuit into a directed acyclic graph (DAG) where each node represents a gate.

    This function constructs a DAG from a Qibo circuit. Each node in the graph corresponds to a gate in the circuit,
    storing information such as gate type, applied qubits, and parameters. Edges represent qubit-based dependencies
    between gates, and optionally, nodes representing identity observables ("I") can be appended to sink nodes in
    the graph.

    :param circuit: Qibo :class:`Circuit` object to convert into a DAG.
    :param obs_I: Optional iterable indicating which qubits have identity ("I") observables appended. Used for internal workflow when cutting. Default None.

    :return: A :class:`networkx.DiGraph` object where nodes represent gates and edges indicate gate dependencies.

    Node Attributes:
        - ``gate``: Name of the gate (e.g., "H", "CNOT").
        - ``qubits``: Tuple of qubit indices the gate acts on.
        - ``parameters``: Tuple of parameters for the gate (e.g., rotation angles).

    Edge Attributes:
        - ``color``: Connection type (e.g., "blue" or "red") used to distinguish dependency types.

    Example:
        .. code-block:: python

            from qibo import models, gates
            import networkx as nx
            import Qdislib as qd

            # Define a Qibo circuit
            c = models.Circuit(2)
            c.add(gates.H(0))
            c.add(gates.CNOT(0, 1))

            # Convert it to a DAG
            dag = qd.circuit_qibo_to_dag(c)

            # Visualize or analyze the DAG
            print(dag.nodes(data=True))
            print(dag.edges(data=True))
    """
    # Create a directed graph
    dag = networkx.DiGraph()

    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.queue, start=1):
        # Unique identifier for each gate instance
        gate_name = f"{gate.__class__.__name__}_{gate_idx}".upper()

        # Add the gate to the DAG, including the gate type, qubits, and parameters
        dag.add_node(
            gate_name,
            gate=gate.__class__.__name__,
            qubits=gate.qubits,
            parameters=gate.parameters,
        )

        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
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

        for qubit in gate.qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and gate_name != pred_gate
                    and not dag.has_edge(pred_gate, gate_name)
                ):
                    dag.add_edge(pred_gate, gate_name, color="red")

    if obs_I:
        for idx, pauli in enumerate(obs_I):
            if pauli != "I":
                continue
            # Attach the Observable-I marker after the LAST gate acting on qubit
            # idx. That gate may be multi-qubit (e.g. a CZ), so requiring a sink
            # whose qubits are exactly (idx,) silently dropped the marker whenever
            # a qubit ended on a two-qubit gate -- and the identity then wrongly
            # became a Z on the cutting path.
            last = None
            for node, data in dag.nodes(data=True):
                if (
                    idx in (data.get("qubits") or ())
                    and data.get("gate") != "Observable I"
                ):
                    last = node
            if last is not None:
                dag.add_node(
                    f"OBSI_{idx}",
                    gate="Observable I",
                    qubits=(idx,),
                    parameters=(),
                )
                dag.add_edge(last, f"OBSI_{idx}", color="blue")

    return dag


def plot_dag(dag: networkx.DiGraph) -> None:
    """
    Visualize a directed acyclic graph (DAG) representation of a quantum circuit.

    This function uses `matplotlib` and `networkx` to plot a DAG where each node represents
    a quantum gate and edges indicate qubit dependencies between gates. Edges can be colored
    differently based on the "color" attribute in the graph:

    - **Blue**: Standard gate dependencies.
    - **Red (dotted)**: Additional or alternative dependencies (e.g., same qubit paths).

    :param dag: A :class:`networkx.DiGraph` object representing a quantum circuit DAG.

    :return: None

    The plot includes:
        - Nodes representing gates, labeled with gate names.
        - Edges with color-coded dependencies.
        - A spring layout for visually appealing positioning.

    Example:
        .. code-block:: python

            import networkx as nx
            import Qdislib as qd

            # Define a Qibo circuit
            c = models.Circuit(2)
            c.add(gates.H(0))
            c.add(gates.CNOT(0, 1))

            # Convert it to a DAG
            dag = qd.circuit_qibo_to_dag(c)

            plot_dag(dag)
    """
    # Set up graph layout
    pos = networkx.spring_layout(dag)

    # Draw edges for the first group with blue color
    edges_first_group = [
        (edge[0], edge[1]) for edge in dag.edges.data("color") if edge[2] == "blue"
    ]
    networkx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_first_group,
        edge_color="blue",
        width=2.0,
        alpha=0.7,
    )

    edges_second_group = [
        (edge[0], edge[1]) for edge in dag.edges.data("color") if edge[2] == "red"
    ]
    networkx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_second_group,
        edge_color="red",
        width=2.0,
        alpha=0.7,
        style="dotted",
    )

    # Draw nodes
    node_labels = networkx.get_node_attributes(dag, "gate")
    networkx.draw_networkx_nodes(dag, pos, node_color="skyblue", node_size=1000)
    networkx.draw_networkx_labels(dag, pos, labels=node_labels, font_size=10)

    # Show plot
    plt.title("Circuit DAG")
    plt.show()


@task(returns=1)
def dag_to_circuit_qibo(dag: networkx.DiGraph, num_qubits: int) -> models.Circuit:
    """
    Create a Qibo quantum circuit from a DAG representation.

    This function takes a directed acyclic graph (DAG) where each node represents a quantum gate,
    and reconstructs a Qibo circuit by adding gates in topological order. Optionally, it detects
    "Observable I" nodes and separates them from the circuit for measurement processing.

    :param dag: A :class:`networkx.DiGraph` representing the quantum circuit DAG. Each node must
                contain `gate`, `qubits`, and `parameters` attributes.
    :param num_qubits: The number of qubits in the original circuit (used to initialize the circuit).

    :return: A list of two elements:
             - A :class:`qibo.models.Circuit` object reconstructed from the DAG.
             - A list of qubit indices where "Observable I" measurements are applied (or `None` if not present). Used for internal workflow when cutting.

    Gate nodes must have the following format:
        - `gate`: Name of the gate (must match a class in `qibo.gates`).
        - `qubits`: A tuple of qubit indices the gate acts on.
        - `parameters`: A tuple of gate parameters (or `None` if not applicable).

    Notes:
        - Measurement gates (e.g., `"Observable I"`) are excluded from reconstruction but their qubit
          indices are returned for postprocessing.
        - Gate classes are resolved case-insensitively (so mixed-case names like ``iSWAP``/``fSim`` round-trip)
          and each stored parameter is splatted as its own argument (so multi-parameter gates like ``U2``/``U3`` work).

    Example:
        .. code-block:: python

            import Qdislib as qd
            import networkx as nx

            dag = ...  # previously constructed DAG
            circuit, observables = qd.dag_to_circuit_qibo(dag, num_qubits=5)
            print(circuit)
            print(observables)
    """
    if dag is None:  # empty component (e.g. a qubit touched only by the cut gate)
        return None

    # Traverse the DAG in topological order
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
    circuit = models.Circuit(highest_qubit)

    topo_order = list(networkx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_raw = node_data["gate"]  # original class name (case preserved)
        gate_name = gate_raw.upper()  # for the case-insensitive marker checks below

        # Skip the observable I nodes (we'll handle them separately)
        if gate_name == "OBSERVABLE I":
            continue

        if gate_name == "MEASURE":
            # Gate-cut measure: a projective Z measurement that collapses the qubit
            # to its outcome |s> and lets it CONTINUE (the CZ cut's (I +/- Z)/2
            # projector, Mitarai-Fujii arXiv:1909.07534). The collapse gate's
            # per-shot samples give the cut-point outcome for the reconstruction
            # sign; the qubit is NOT reset (the old RX force discarded the collapse
            # and broke reconstruction for gates after the cut).
            qubits = node_data["qubits"]
            circuit.add(gates.M(*qubits, collapse=True))
            continue

        # Get the qubits this gate acts on
        qubits = node_data["qubits"]
        params = _unpack_parameters(node_data["parameters"])

        # Resolve the gate class case-insensitively (handles iSWAP, fSim, ...) and
        # splat each parameter as its own positional argument (handles U2, U3, ...).
        gate_class = _resolve_qibo_gate(gate_raw)
        if params:
            circuit.add(gate_class(*qubits, *params))
        else:
            circuit.add(gate_class(*qubits))

    if obs_I:
        observables = ["Z"] * circuit.nqubits
        for element in obs_I:
            observables[element] = "I"

        observables = "".join(observables)
        return [circuit, observables]

    observables = ["Z"] * circuit.nqubits
    observables = "".join(observables)
    return [circuit, observables]


def _max_qubit(graph: networkx.DiGraph) -> float:
    """Get the highest Qubit value.

    :param graph: Graph to explore.
    :return: The value of the highest Qubit.
    """
    # Initialize a variable to keep track of the highest Qubits value
    max_qubits = float("-inf")  # Start with the lowest possible number

    # Iterate over the nodes and check their 'Qubits' attribute
    for node, data in graph.nodes(data=True):
        qubits = data.get("qubits", 0)  # Default to 0 if 'Qubits' is not present
        for qubit in qubits:
            if qubit > max_qubits:
                max_qubits = qubit
                # Keep track of the node with the highest Qubit value
    return max_qubits


# Not a @task: it is a cheap in-place qubit reindex. As a plain function it runs
# inline inside the per-subcircuit task (all-in-one) and returns a concrete graph
# on the master (3-step), so its result can be wrapped in a serializable
# Subcircuit without a fragile compss_wait_on of an INOUT task output.
def _update_qubits(
    s: typing.List[typing.Any],
) -> typing.Tuple[typing.List[typing.Any], float]:
    """Update qubits task.

    :param s: Graph to explore.
    :return: The updated graph and the highest qubit value.
    """
    if not s.nodes():
        return None, None
    my_set = set()
    for node, _ in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    for node, _ in s.nodes(data=True):
        new_tuple = ()
        for qubit in s.nodes[node]["qubits"]:
            len_missing = _count_missing_up_to(my_set, qubit)
            new_qubit = qubit - len_missing
            new_tuple = new_tuple + (new_qubit,)
        s.nodes[node]["qubits"] = new_tuple

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    return s, highest_qubit


def _remove_red_edges(graph: networkx.DiGraph) -> networkx.DiGraph:
    """Remove red edges from the given graph.

    :param graph: Graph to process.
    :return: Updated input graph without red edges.
    """
    copy_dag = graph.copy()
    red_edges = []

    for ed in copy_dag.edges:
        # A few reconnection edges are added without a colour attribute; treat a
        # missing colour as "not red" (keep it) rather than raising.
        if copy_dag.get_edge_data(ed[0], ed[1]).get("color") == "red":
            red_edges.append(ed)

    copy_dag.remove_edges_from(red_edges)
    return copy_dag


def _count_missing_up_to(nums: typing.Set[int], max_num: int) -> int:
    """Retrieve the amount of missing numbers in the nums set.

    :param nums: Set of numbers.
    :param max_num: Highest number in nums.
    :return: Amount of missing numbers.
    """
    # Create a set of all numbers from 0 to max_num
    full_set = set(range(max_num + 1))

    # Subtract the given set from the full set to get the missing numbers
    missing_numbers = full_set - nums

    # Return the count of missing numbers
    return len(missing_numbers)


def _update_qubits_serie(
    s: typing.List[typing.Any],
) -> typing.Tuple[typing.List[typing.Any], int, int]:
    """Update the given serie of qubits.

    :param s: Input qubit serie.
    :return: Updated qubit serie, highest qubit value and lowest qubit value.
    """
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    for node, data in s.nodes(data=True):
        new_tuple = ()
        for qubit in s.nodes[node]["qubits"]:
            len_missing = _count_missing_up_to(my_set, qubit)
            new_qubit = qubit - len_missing
            new_tuple = new_tuple + (new_qubit,)
        s.nodes[node]["qubits"] = new_tuple

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    smallest_qubit = min(my_set) - _count_missing_up_to(my_set, min(my_set))
    return s, highest_qubit, smallest_qubit


def _max_qubits_graph(s: typing.Any) -> int:
    """Return the number of distinct qubits touched by the gates in a DAG.

    :param s: A circuit DAG whose nodes carry a ``"qubits"`` attribute.
    :return: The qubit count (highest contiguous qubit index after compaction).
    """
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    return highest_qubit
