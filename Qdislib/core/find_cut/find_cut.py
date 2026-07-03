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

"""Find cut algorithms."""

import logging

logger = logging.getLogger(__name__)

__all__ = ["find_cut", "draw_cut"]

import networkx
import typing

from Qdislib.utils.runtime import task, compss_wait_on

try:
    import qibo

    qibo_available = True
except ImportError:  # qibo is optional
    qibo = None
    qibo_available = False

try:
    import qiskit
    import qiskit.qasm2

    qiskit_available = True
except ImportError:  # qiskit is optional
    qiskit = None
    qiskit_available = False

import pymetis

from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
from Qdislib.utils.backends import iter_gates, num_qubits

try:
    from qiskit_addon_cutting.automated_cut_finding import (
        find_cuts,
        OptimizationParameters,
        DeviceConstraints,
    )
except ImportError:  # qiskit-addon-cutting optional (only the ibm-ckt path)
    find_cuts = OptimizationParameters = DeviceConstraints = None


def _find_nodes_with_qubit(
    graph: networkx.Graph,
    node: typing.Any,
    qubit: typing.Any,
    direction: str = "predecessor",
) -> typing.List[typing.Any]:
    """Find predecessor or successor nodes with specific qubit.

    :param graph: Graph.
    :param node: Node.
    :param qubit: Qubit.
    :param direction: Direction, defaults to "predecessor".
    :raises ValueError: Unexpected direction error. Must be "predecessor" or "successor".
    :return: List of nodes with qubit.
    """
    if direction == "predecessor":
        neighbors = graph.predecessors(node)
    elif direction == "successor":
        neighbors = graph.successors(node)
    else:
        raise ValueError("Direction must be either 'predecessor' or 'successor'")

    # Filter neighbors based on the qubit data
    nodes_with_qubit = [n for n in neighbors if qubit in graph.nodes[n]["qubits"]]
    return nodes_with_qubit


def _circuit_to_qubit_graph(circuit):
    """Build the qubit-connectivity graph used for *gate* cut finding.

    Nodes are qubits and each multi-qubit gate becomes an edge between the qubits
    it couples, labelled with the gate's ``NAME_index`` identifier. Partitioning
    this graph and cutting the edges that cross the partition is exactly gate
    cutting. Backend-agnostic (Qibo / Qiskit / CUDA-Q) via :func:`iter_gates`.

    Parallel edges are intentional here: two distinct gates acting on the same
    qubit pair are two distinct gate cuts, so both belong in the cut-set.

    :param circuit: A Qibo, Qiskit or CUDA-Q circuit.
    :return: A :class:`networkx.MultiGraph` of qubit connectivity.
    """
    graph = networkx.MultiGraph()
    for qubit in range(num_qubits(circuit)):
        graph.add_node(qubit, gate=qubit)

    for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
        if len(qubits) > 1:
            graph.add_edge(*qubits, gate=f"{name}_{index}")

    return graph


def _circuit_to_gate_graph(circuit):
    """Build the gate-dependency graph used for *wire* cut finding.

    Nodes are the circuit's multi-qubit gates; an edge links a gate to the most
    recent earlier gate that shares a qubit with it. Cutting an edge corresponds
    to a wire cut between those two gates. Backend-agnostic (Qibo / Qiskit /
    CUDA-Q) via :func:`iter_gates`.

    A pair of gates is connected by at most one edge. When two consecutive gates
    share more than one qubit, adding one edge per shared qubit would create
    duplicate parallel edges carrying the same ``(pred, gate)`` label; those
    duplicates inflate the cut-set size and yield spurious repeated wire cuts.
    The ``has_edge`` guard below keeps a single edge per gate pair.

    :param circuit: A Qibo, Qiskit or CUDA-Q circuit.
    :return: A :class:`networkx.MultiGraph` whose nodes carry ``gate`` (name) and
        ``qubits`` attributes, and whose edges carry a ``(pred_name, gate_name)``
        ``gate`` label.
    """
    graph = networkx.MultiGraph()

    node_id = 0
    for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
        if len(qubits) < 2:
            continue

        gate_name = f"{name}_{index}"
        graph.add_node(node_id, gate=gate_name, qubits=qubits)

        # Connect to the nearest earlier gate on each shared qubit, but never add
        # more than one edge between the same pair of gates (see docstring).
        for qubit in qubits:
            for pred in reversed(list(graph.nodes)):
                pred_qubits = graph.nodes[pred].get("qubits")
                if pred_qubits and qubit in pred_qubits and pred != node_id:
                    if not graph.has_edge(pred, node_id):
                        pred_name = graph.nodes[pred]["gate"]
                        graph.add_edge(pred, node_id, gate=(pred_name, gate_name))
                    break

        node_id += 1

    return graph


def _find_best_kernighan_lin_partition(
    graph, max_qubits, max_cuts=8, num_runs=10, seed=None
):
    best_cutset_size = float("inf")
    best_partition = []
    best_cutset = []

    for run in range(num_runs):
        # Perform balanced partition using Kernighan-Lin. On graphs too small to
        # bisect (e.g. a circuit with one or two gates) networkx raises
        # ValueError ("min() arg is an empty sequence"); treat that as a run that
        # produced no valid partition rather than crashing find_cut.
        #
        # A distinct per-run seed keeps the runs diverse (KL is randomised) while
        # making the whole search reproducible when a seed is given; seed=None
        # preserves the previous non-deterministic behaviour.
        run_seed = None if seed is None else seed + run
        try:
            partition = networkx.algorithms.community.kernighan_lin_bisection(
                graph.to_undirected(), seed=run_seed
            )
        except ValueError:
            continue

        # Extract the edges in the cut-set
        cutset = list(
            networkx.edge_boundary(graph, partition[0], partition[1], data=True)
        )
        cutset_size = len(cutset)
        flag = False
        for ele in sorted(list(c) for c in partition):
            if max_qubits and len(ele) > max_qubits:
                flag = True
        if max_cuts and cutset_size > max_cuts:
            flag = True
        # This run violates the qubit/cut budget: skip it and try the next
        # random partition rather than abandoning the whole search.
        if flag:
            continue

        # If this cut-set is smaller, update the best partition and cut-set
        if cutset_size < best_cutset_size:
            best_cutset_size = cutset_size
            best_partition = partition
            best_cutset = cutset

    if best_cutset:
        best_cutset = [el[2]["gate"] for el in best_cutset]
    else:
        return [], []
    return best_partition, best_cutset


def _find_best_girvan_newman_cutset(graph, max_qubits, max_cuts):
    from networkx.algorithms.community import girvan_newman

    # Perform Girvan-Newman community detection
    communities = girvan_newman(graph)

    # Defaults in case no partition is produced (e.g. a graph with < 2 nodes,
    # where girvan_newman yields nothing and the loop body never runs).
    partition_list, cut_edges = [], []
    # Get multiple partitions at different levels of edge removal
    for i, partition in enumerate(communities):
        partition_list = list(partition)  # Convert partition to list of sets
        # Find edges that are between different partitions
        cut_edges = []
        for u, v, gate in graph.edges(data=True):
            # Check if nodes u and v belong to different components
            component_u = [c for c in partition_list if u in c][0]
            component_v = [c for c in partition_list if v in c][0]
            if component_u != component_v:
                # data = graph.get_edge_data(u, v)
                # gates_list = tuple([value["gate"] for value in data.values()])
                cut_edges = cut_edges + [gate["gate"]]

        flag = True
        for ele in sorted(list(c) for c in partition_list):
            if max_qubits and len(ele) > max_qubits:
                flag = False
        if flag:
            break

    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return partition_list, cut_edges


def _find_best_spectral_clustering(graph, clusters, max_qubits, max_cuts, seed=None):
    from sklearn.cluster import SpectralClustering

    # Spectral clustering requires at least as many nodes as requested clusters;
    # clamp so an over-large target can never raise (needs n_samples >= n_clusters).
    num_clusters = min(clusters, graph.number_of_nodes())
    if num_clusters < 2:
        return [], []

    # Convert the graph to adjacency matrix form
    adjacency_matrix = networkx.to_numpy_array(graph)

    sc = SpectralClustering(
        n_clusters=num_clusters,
        affinity="precomputed",
        assign_labels="kmeans",
        random_state=seed,
    )
    labels = sc.fit_predict(adjacency_matrix)

    # Group nodes based on clusters
    clusters = {}
    for node, label in enumerate(labels):
        clusters.setdefault(label, []).append(node)

    flag = False
    for cluster_id, nodes in clusters.items():
        if max_qubits and len(nodes) > max_qubits:
            flag = True

    if flag:
        return [], []

    clusters = list(clusters.values())

    # Find the cut edges: edges between nodes in different clusters
    cut_edges = []
    for u, v, gate in graph.edges(data=True):
        if labels[u] != labels[v]:
            # data = graph.get_edge_data(u, v)
            # gates_list = [value["gate"] for value in data.values()]
            cut_edges = cut_edges + [gate["gate"]]

    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return clusters, cut_edges


def _metis_partition(graph, num_clusters, max_nodes_per_cluster, max_cuts):
    # Clamp the partition count to the number of nodes (METIS needs nparts <= n).
    num_clusters = min(num_clusters, graph.number_of_nodes())
    if num_clusters < 2:
        return [], []

    # Convert the NetworkX graph to a format suitable for METIS
    adjacency_list = []
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        adjacency_list.append(neighbors)
    _, parts = pymetis.part_graph(adjacency=adjacency_list, nparts=num_clusters)
    # Group nodes based on the partition labels
    clusters = {}
    for node, part in enumerate(parts):
        clusters.setdefault(part, []).append(node)

    # Check for cluster size constraint violations
    for cluster_id, nodes in clusters.items():
        if max_nodes_per_cluster and len(nodes) > max_nodes_per_cluster:
            logger.info(f"Cluster {cluster_id} exceeds the max size.")
            return [], []

    clusters = list(clusters.values())
    cut_edges = []
    for u, v, gate in graph.edges(data=True):
        if parts[u] != parts[v]:
            # data = graph.get_edge_data(u, v)
            # gates_list = [value["gate"] for value in data.values()]
            cut_edges = cut_edges + [gate["gate"]]
    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return clusters, cut_edges


def _score_candidates(candidates, weights=None, forbidden=None):
    """Pick the lowest-scoring partitioning candidate.

    Each candidate maps a method name to a ``(partition, cut_set)`` pair produced
    by one of the graph-partitioning helpers. The objective (lower is better)
    balances three terms::

        score = | w1 * cut_size  +  w2 * num_components  +  w3 * largest_component |

    rewarding few cuts, more components, and small (balanced) components.

    :param candidates: Ordered mapping ``{name: (partition, cut_set)}``.
    :param weights: Optional overrides for the ``w1``/``w2``/``w3`` coefficients.
    :param forbidden: Gate identifiers the user forbade cutting; any candidate
        whose cut-set touches one is discarded, so the winner is the best cut that
        respects the constraint (or ``[]`` if none does).
    :return: ``(best_cut_set, best_score, best_partition)`` for the winner.

    Candidate scores are logged at DEBUG level (``qd.set_log_level("DEBUG")``).
    """
    weights = weights or {}
    w1 = weights.get("w1", 100)  # cut size            (minimize)
    w2 = weights.get("w2", -1)  # number of components (maximize)
    w3 = weights.get("w3", 1)  # largest component   (minimize -> balance)
    forbidden = set(forbidden or ())
    logger.debug("Using weights: w1=%s, w2=%s, w3=%s", w1, w2, w3)

    best_score = float("inf")
    best_cut_set = []
    best_partition = []
    for name, (partition, cut_set) in candidates.items():
        if forbidden and forbidden.intersection(cut_set):
            logger.debug("%s: skipped (cuts a forbidden gate)", name)
            continue
        largest_component = max(
            [len(component) for component in partition] + [float("-inf")]
        )
        score = abs(w1 * len(cut_set) + w2 * len(partition) + w3 * largest_component)
        logger.debug(
            "%s: cut_set=%s components=%s largest=%s score=%s",
            name,
            cut_set,
            len(partition),
            largest_component,
            score,
        )
        if score < best_score:
            best_score = score
            best_cut_set = cut_set
            best_partition = partition
    return best_cut_set, best_score, best_partition


def _plot_partition(graph, partition):
    """Draw a partitioned circuit graph (manual debug aid; call it yourself)."""
    import matplotlib.pyplot as plt

    pos = networkx.spring_layout(graph)
    plt.figure(figsize=(12, 12))
    colors = ["orange", "purple", "green", "red", "blue"]
    for i, component in enumerate(partition):
        networkx.draw_networkx_nodes(
            graph,
            pos,
            nodelist=list(component),
            node_color=colors[i % len(colors)],
            node_size=200,
            alpha=0.5,
        )
    networkx.draw_networkx_edges(graph, pos, alpha=0.5)
    networkx.draw_networkx_labels(graph, pos)
    plt.title("Graph Partition Visualization")
    plt.show()


def draw_cut(circuit, cut, filename=None):
    """Visualise a circuit's DAG with the cut locations highlighted.

    Gate cuts (gate identifiers, e.g. ``["CZ_3"]``) highlight the cut gate nodes;
    wire cuts (endpoint pairs, e.g. ``[("H_1", "CZ_2")]``) also highlight the edge
    between the endpoints. Requires matplotlib.

    :param circuit: The circuit that was cut (Qibo/Qiskit/CUDA-Q, or a DAG).
    :param cut: A cut as returned by :func:`find_cut`.
    :param filename: If given, save the figure there; otherwise show it.
    """
    import matplotlib.pyplot as plt
    from Qdislib.utils.backends import convert_to_dag

    dag = convert_to_dag(circuit)
    if not isinstance(dag, networkx.DiGraph):
        raise TypeError("draw_cut requires a circuit convertible to a Qdislib DAG.")

    is_wire = bool(cut) and isinstance(cut[0], (tuple, list))
    cut_nodes = set()
    cut_edges = set()
    if is_wire:
        for source, target in cut:
            cut_edges.add((source, target))
            cut_nodes.update((source, target))
    else:
        cut_nodes.update(cut)

    pos = networkx.spring_layout(dag)
    plt.figure(figsize=(12, 8))
    node_colors = ["#d62728" if n in cut_nodes else "#aec7e8" for n in dag.nodes]
    edge_colors = [
        "#d62728" if (u, v) in cut_edges or (v, u) in cut_edges else "#bbbbbb"
        for u, v in dag.edges
    ]
    networkx.draw(
        dag,
        pos,
        with_labels=True,
        node_color=node_colors,
        edge_color=edge_colors,
        node_size=600,
        font_size=7,
    )
    plt.title(f"{'Wire' if is_wire else 'Gate'} cut: {cut}")
    if filename:
        plt.savefig(filename, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


@task(returns=2)
def _find_cut_gate(
    circuit,
    max_qubits=None,
    max_cuts=None,
    max_components=None,
    weights=None,
    seed=None,
    forbidden_gates=None,
):
    """Find the best *gate* cut for ``circuit`` by partitioning its qubit graph.

    Builds the qubit-connectivity graph, runs the four partitioning strategies
    (Kernighan-Lin, Girvan-Newman, spectral clustering, METIS) and keeps the
    lowest-scoring cut. Backend-agnostic (Qibo / Qiskit / CUDA-Q). ``seed`` makes
    the randomised strategies (Kernighan-Lin, spectral) reproducible.

    :return: ``(best_cut, best_score)`` where ``best_cut`` is a list of gate
        identifiers such as ``["CZ_3"]``.
    """
    logger.debug("Finding best gate cut")

    graph = _circuit_to_qubit_graph(circuit)
    candidates = {
        "Kernighan-Lin": _find_best_kernighan_lin_partition(
            graph=graph, max_qubits=max_qubits, max_cuts=max_cuts, seed=seed
        ),
        "Girvan-Newman": _find_best_girvan_newman_cutset(
            graph, max_qubits, max_cuts=max_cuts
        ),
        "Spectral clustering": _find_best_spectral_clustering(
            graph,
            clusters=max_components,
            max_qubits=max_qubits,
            max_cuts=max_cuts,
            seed=seed,
        ),
        "METIS": _metis_partition(
            graph,
            num_clusters=max_components,
            max_nodes_per_cluster=max_qubits,
            max_cuts=max_cuts,
        ),
    }

    best_cut, best_score, _ = _score_candidates(
        candidates, weights, forbidden=forbidden_gates
    )
    logger.debug("Best gate cut: %s (score %s)", best_cut, best_score)
    return best_cut, best_score


@task(returns=2)
def _find_cut_wire(
    circuit,
    max_qubits=None,
    max_cuts=None,
    max_components=None,
    weights=None,
    seed=None,
    forbidden_gates=None,
):
    """Find the best *wire* cut for ``circuit`` by partitioning its gate graph.

    Builds the gate-dependency graph, runs the four partitioning strategies and
    keeps the lowest-scoring cut. The winning cut-set (pairs of gate identifiers)
    is then refined into concrete wire cuts: for each pair, any gate sitting
    between the two endpoints that shares a qubit with both becomes the actual
    cut location. Backend-agnostic (Qibo / Qiskit / CUDA-Q).

    :return: ``(best_cut, best_score)`` where ``best_cut`` is a list of
        ``(source_gate, target_gate)`` identifier pairs.
    """
    logger.debug("Finding best wire cut")

    graph = _circuit_to_gate_graph(circuit)
    candidates = {
        "Kernighan-Lin": _find_best_kernighan_lin_partition(
            graph=graph, max_qubits=max_qubits, max_cuts=max_cuts, seed=seed
        ),
        "Girvan-Newman": _find_best_girvan_newman_cutset(
            graph, max_qubits=max_qubits, max_cuts=max_cuts
        ),
        "Spectral clustering": _find_best_spectral_clustering(
            graph,
            clusters=max_components,
            max_qubits=max_qubits,
            max_cuts=max_cuts,
            seed=seed,
        ),
        "METIS": _metis_partition(
            graph,
            num_clusters=max_components,
            max_nodes_per_cluster=max_qubits,
            max_cuts=max_cuts,
        ),
    }

    cut_set, best_score, _ = _score_candidates(
        candidates, weights, forbidden=forbidden_gates
    )

    # Refine each cut pair into a concrete wire cut: locate a gate between the
    # two endpoints that shares a qubit with both, and cut just after it.
    gates = iter_gates(circuit)
    best_cut = []
    for source, target in cut_set:
        source_idx = int(source.split("_")[1])
        target_idx = int(target.split("_")[1])
        _, source_qubits = gates[source_idx - 1]
        _, target_qubits = gates[target_idx - 1]

        refined = False
        for idx, (gate_name, gate_qubits) in enumerate(
            gates[source_idx : target_idx - 1], start=source_idx
        ):
            if set(gate_qubits).intersection(source_qubits, target_qubits):
                best_cut.append((source, f"{gate_name}_{idx + 1}"))
                refined = True
        if not refined:
            best_cut.append((source, target))

    logger.debug("Best wire cut: %s (score %s)", best_cut, best_score)
    return best_cut, best_score


def find_cut(
    circuit,
    max_qubits=None,
    max_cuts=None,
    wire_cut=True,
    gate_cut=True,
    implementation="qdislib",
    weights=None,
    seed=None,
    forbidden_gates=None,
    forbidden_qubits=None,
):
    """
    Automatically partitions a quantum circuit to optimize for parallel execution
    and hardware constraints, using either the Qdislib or IBM Qiskit-based backends.

    The FindCut algorithm from the Qdislib library is designed to find the best cut to
    decompose a quantum circuit into subcircuits that adhere to specific
    constraints such as the number of qubits per subcircuit, the number of allowed cuts,
    and the number of resulting components. It supports both wire cutting and gate cutting
    strategies, which are evaluated through multiple graph partitioning techniques
    (Kernighan-Lin, Girvan-Newman, Spectral Decomposition, METIS). It scores all possible options
    with a loss function and then decides which cut is better for that specific circuit and constraints.

    Parameters
    ----------
    circuit : QuantumCircuit, qibo.models.Circuit, or networkx.DiGraph
        The input quantum circuit to be partitioned.

    max_qubits : int, optional
        Maximum number of qubits allowed per subcircuit. Helps enforce hardware
        constraints (e.g., maximum qubits per QPU). When given, it also determines
        how many subcircuits the clustering strategies aim for
        (``ceil(num_qubits / max_qubits)``); otherwise a simple bisection is used.

    max_cuts : int, optional
        Maximum number of cuts permitted across the circuit. Controls the complexity
        and computational overhead of the partitioning.

    wire_cut : bool, default=True
        Whether to consider wire cuts (cutting dependencies between gates)
        during partitioning.

    gate_cut : bool, default=True
        Whether to consider gate cuts (cutting the circuit by splitting two-qubit gates)
        during partitioning.

    implementation : str, default="qdislib"
        Specifies the backend implementation to use for circuit partitioning. "qdislib": Uses graph-theoretical approaches (Kernighan-Lin, etc.) "ibm-ckt": Uses Qiskit's optimizer and hardware-aware constraints

    weights : dict, optional
        A dictionary specifying weights for the loss function used in Qdislib to evaluate
        partitioning

    seed : int, optional
        Seed for the randomised partitioning strategies (Kernighan-Lin and spectral
        clustering). Pass an integer for a reproducible cut; the default (``None``)
        keeps the previous non-deterministic behaviour. Only affects the
        ``implementation="qdislib"`` path.

    forbidden_gates : list of str, optional
        Gate identifiers (e.g. ``["CZ_3"]``) that must **not** be cut. Candidate
        cuts touching any of them are discarded, so the result is the best cut that
        leaves those gates intact (or ``[]`` if none does). Only affects the
        ``implementation="qdislib"`` path.

    forbidden_qubits : list of int, optional
        Qubit indices that must not be cut: every multi-qubit gate acting on one of
        them is added to ``forbidden_gates``, so no cut touches those qubits. Only
        affects the ``implementation="qdislib"`` path.

    Diagnostic details about the partitioning process are logged at DEBUG level;
    enable them with ``qd.set_log_level("DEBUG")``.

    Returns
    -------
    best_cut : list
        A list describing the best cut found. This is a list of gate names and their indices where cuts should be applied (e.g., ['CX_3']).

    Raises
    ------
    NameError
        If neither `wire_cut` nor `gate_cut` is enabled.

    TypeError
        If the chosen implementation fails to return a valid result due to configuration.

    Notes
    -----
    - In the Qdislib implementation:
        * Circuits are transformed into graphs and evaluated through multiple graph partitioning
          strategies. The best partition is chosen based on a custom loss function.
        * The function supports concurrent evaluation of both gate and wire cuts.
    - In the IBM implementation:
        * The circuit is compiled for compatibility with Qiskit.
        * Gate cuts are identified using an optimizer, subject to QPU constraints.
        * The returned result focuses on gate-level decomposition only.

    Example
    -------
    >>> import Qdislib as qd
    >>> cut = find_cut(my_circuit, max_qubits=5, gate_cut=True, wire_cut=False)
    >>> print(cut)  # Output: Best cut information (e.g., ['CX_3'])
    """

    if implementation == "qdislib":
        # Target number of subcircuits for the clustering-based strategies
        # (spectral, METIS). Derived from max_qubits when given (enough pieces so
        # each fits the qubit budget), otherwise a simple bisection. The
        # clustering helpers clamp this to the graph size, so it is only an upper
        # bound and can never request more clusters than there are nodes.
        if max_qubits:
            max_components = max(2, -(-num_qubits(circuit) // max_qubits))
        else:
            max_components = 2

        # Fold forbidden_qubits into the forbidden-gate set: any multi-qubit gate
        # touching a forbidden qubit becomes a forbidden cut. The gate id is built
        # exactly as the qubit/gate graphs label their cut edges (``NAME_index``
        # over iter_gates), so the names line up with the candidate cut-sets.
        forbidden = set(forbidden_gates or ())
        if forbidden_qubits:
            protected = set(forbidden_qubits)
            for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
                if len(qubits) > 1 and protected.intersection(qubits):
                    forbidden.add(f"{name}_{index}")

        if gate_cut:
            cut_gate, score_gate = _find_cut_gate(
                circuit,
                max_qubits=max_qubits,
                max_cuts=max_cuts,
                max_components=max_components,
                weights=weights,
                seed=seed,
                forbidden_gates=forbidden,
            )

        if wire_cut:
            cut_wire, score_wire = _find_cut_wire(
                circuit,
                max_qubits=max_qubits,
                max_cuts=max_cuts,
                max_components=max_components,
                weights=weights,
                seed=seed,
                forbidden_gates=forbidden,
            )

        if gate_cut and wire_cut:
            cut_gate = compss_wait_on(cut_gate)
            score_gate = compss_wait_on(score_gate)
            cut_wire = compss_wait_on(cut_wire)
            score_wire = compss_wait_on(score_wire)
        elif gate_cut:
            cut_gate = compss_wait_on(cut_gate)
            score_gate = compss_wait_on(score_gate)
        elif wire_cut:
            cut_wire = compss_wait_on(cut_wire)
            score_wire = compss_wait_on(score_wire)
        else:
            raise ValueError("Enable at least one of wire_cut or gate_cut.")

        if wire_cut and gate_cut:
            best = cut_gate if score_gate <= score_wire else cut_wire
        elif gate_cut:
            best = cut_gate
        elif wire_cut:
            best = cut_wire
        else:
            raise ValueError("Enable at least one of wire_cut or gate_cut.")

        if not best:
            logger.warning(
                "find_cut found no valid cut under the given constraints "
                "(max_qubits=%s, max_cuts=%s): the circuit may already fit on the "
                "target as-is, or the constraints may be too tight to satisfy. Try "
                "relaxing max_qubits / max_cuts.",
                max_qubits,
                max_cuts,
            )
        return best

    elif implementation == "ibm-ckt":
        if find_cuts is None:
            raise ImportError(
                "implementation='ibm-ckt' needs the qiskit-addon-cutting package, "
                "which is not installed. Install it with "
                "`pip install qiskit-addon-cutting`, or use the default "
                "implementation='qdislib'."
            )
        if qibo_available and isinstance(circuit, qibo.models.Circuit):
            qasm_str = circuit.to_qasm()
            circuit = qiskit.qasm2.loads(qasm_str)
        elif isinstance(circuit, networkx.DiGraph):
            circuit = dag_to_circuit_qiskit(circuit)

        # Specify settings for the cut-finding optimizer
        optimization_settings = OptimizationParameters(seed=50)

        # Specify the size of the QPUs available
        if max_qubits is None:
            max_qubits_ckt = circuit.num_qubits // 2 + 1
        else:
            max_qubits_ckt = max_qubits
        device_constraints = DeviceConstraints(qubits_per_subcircuit=max_qubits_ckt)

        cut_circuit, metadata = find_cuts(
            circuit, optimization_settings, device_constraints
        )
        logger.debug(
            "Found solution using %s cuts with a sampling overhead of %s. "
            "Lowest cost solution found: %s.",
            len(metadata["cuts"]),
            metadata["sampling_overhead"],
            metadata["minimum_reached"],
        )
        for cut in metadata["cuts"]:
            logger.debug("%s at circuit instruction index %s", cut[0], cut[1])

        list_gates_cut = []
        for cut in metadata["cuts"]:
            if cut[0] == "Gate Cut":
                gate_name = circuit.data[cut[1]].operation.name.upper()
                list_gates_cut.append(f"{gate_name}_{cut[1]+1}")
        logger.debug("Gate cuts: %s", list_gates_cut)
        return list_gates_cut

    else:
        raise ValueError(
            f"Unknown implementation={implementation!r}; "
            "use 'qdislib' (default) or 'ibm-ckt'."
        )
