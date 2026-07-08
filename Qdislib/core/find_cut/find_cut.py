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

    A pair of gates is connected by at most one edge, even when they share more
    than one qubit: parallel edges carrying the same ``(pred, gate)`` label would
    break the min-cut partitioners and yield spurious repeated wire cuts. Instead
    the single edge records how many qubits the pair shares as its ``weight`` --
    the number of wires a cut between them actually severs (``_expand_wire_cuts``
    cuts one per shared qubit at execution time). Cut-size scoring sums these
    weights (see :func:`_gate_edge_weights` / :func:`_rank_candidates`), so a pair
    sharing two qubits is correctly costed as two wire cuts, not one.

    :param circuit: A Qibo, Qiskit or CUDA-Q circuit.
    :return: A :class:`networkx.MultiGraph` whose nodes carry ``gate`` (name) and
        ``qubits`` attributes, and whose edges carry a ``(pred_name, gate_name)``
        ``gate`` label plus a ``weight`` (the shared-qubit count).
    """
    graph = networkx.MultiGraph()

    node_id = 0
    for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
        if len(qubits) < 2:
            continue

        gate_name = f"{name}_{index}"
        graph.add_node(node_id, gate=gate_name, qubits=qubits)

        # Connect to the nearest earlier gate on each shared qubit. Keep a single
        # edge per gate pair (parallel edges break the partitioners), but bump its
        # weight for every additional shared qubit -- that weight is how many wires
        # a cut there severs.
        for qubit in qubits:
            for pred in reversed(list(graph.nodes)):
                pred_qubits = graph.nodes[pred].get("qubits")
                if pred_qubits and qubit in pred_qubits and pred != node_id:
                    if graph.has_edge(pred, node_id):
                        key = next(iter(graph[pred][node_id]))
                        graph[pred][node_id][key]["weight"] += 1
                    else:
                        pred_name = graph.nodes[pred]["gate"]
                        graph.add_edge(
                            pred, node_id, gate=(pred_name, gate_name), weight=1
                        )
                    break

        node_id += 1

    return graph


def _gate_edge_weights(graph):
    """Map each gate-graph edge's ``(pred, gate)`` label to its wire multiplicity.

    Built from :func:`_circuit_to_gate_graph`'s edge ``weight`` (the number of
    qubits the two gates share). :func:`_rank_candidates` sums these over a cut-set
    so a cut between gates sharing *n* qubits is scored as *n* wire cuts. Gate cuts
    (qubit-graph edges, single-gate labels) are absent here and default to 1.
    """
    return {
        data["gate"]: data.get("weight", 1)
        for _, _, data in graph.edges(data=True)
        if "gate" in data
    }


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
        # Count wires, not edges: an edge between gates sharing two qubits is two
        # wire cuts (its weight), so both the budget check and best-cut selection
        # use the summed weight.
        cutset_size = sum(edge[2].get("weight", 1) for edge in cutset)
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
    partition_list, cut_edges, cut_wires = [], [], 0
    # Get multiple partitions at different levels of edge removal
    for i, partition in enumerate(communities):
        partition_list = list(partition)  # Convert partition to list of sets
        # Find edges that are between different partitions
        cut_edges = []
        cut_wires = 0  # wire count = summed edge weight (shared-qubit multiplicity)
        for u, v, gate in graph.edges(data=True):
            # Check if nodes u and v belong to different components
            component_u = [c for c in partition_list if u in c][0]
            component_v = [c for c in partition_list if v in c][0]
            if component_u != component_v:
                cut_edges = cut_edges + [gate["gate"]]
                cut_wires += gate.get("weight", 1)

        flag = True
        for ele in sorted(list(c) for c in partition_list):
            if max_qubits and len(ele) > max_qubits:
                flag = False
        if flag:
            break

    if max_cuts and cut_wires > max_cuts:
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


def _score_candidates(candidates, weights=None, forbidden=None, edge_weights=None):
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
    ranked = _rank_candidates(candidates, weights, forbidden, edge_weights)
    return ranked[0] if ranked else ([], float("inf"), [])


def _rank_candidates(candidates, weights=None, forbidden=None, edge_weights=None):
    """All valid partitioning candidates, lowest (best) score first.

    Uses the same objective as :func:`_score_candidates`; drops empty cut-sets
    and any candidate touching a ``forbidden`` gate. Powers both ``find_cut``'s
    single-best result and its ranked ``top_n`` list.

    :param edge_weights: Optional ``{cut-label: wire count}`` map (from
        :func:`_gate_edge_weights`) so a wire cut between gates sharing *n* qubits
        counts as *n* toward the cut-size term; labels missing from it default to 1
        (e.g. gate cuts, which are always a single cut).
    :return: ``[(cut_set, score, partition), ...]`` sorted best-first.
    """
    weights = weights or {}
    edge_weights = edge_weights or {}
    w1 = weights.get("w1", 100)  # cut size            (minimize)
    w2 = weights.get("w2", -1)  # number of components (maximize)
    w3 = weights.get("w3", 1)  # largest component   (minimize -> balance)
    forbidden = set(forbidden or ())
    logger.debug("Using weights: w1=%s, w2=%s, w3=%s", w1, w2, w3)

    scored = []
    for name, (partition, cut_set) in candidates.items():
        if not cut_set:
            continue
        if forbidden and forbidden.intersection(cut_set):
            logger.debug("%s: skipped (cuts a forbidden gate)", name)
            continue
        largest_component = max(
            [len(component) for component in partition] + [float("-inf")]
        )
        # Cut size counts severed WIRES, not graph edges: a cut between two gates
        # that share two qubits weighs 2 (see _gate_edge_weights).
        cut_size = sum(edge_weights.get(cut, 1) for cut in cut_set)
        score = abs(w1 * cut_size + w2 * len(partition) + w3 * largest_component)
        logger.debug(
            "%s: cut_set=%s wires=%s components=%s largest=%s score=%s",
            name,
            cut_set,
            cut_size,
            len(partition),
            largest_component,
            score,
        )
        scored.append((cut_set, score, partition))
    scored.sort(key=lambda item: item[1])
    return scored


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
    all_candidates=False,
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

    if all_candidates:
        ranked = _rank_candidates(candidates, weights, forbidden_gates)
        return [(cut, score) for cut, score, _ in ranked], None

    best_cut, best_score, _ = _score_candidates(
        candidates, weights, forbidden=forbidden_gates
    )
    logger.debug("Best gate cut: %s (score %s)", best_cut, best_score)
    return best_cut, best_score


def _refine_wire_cut(cut_set, gates):
    """Refine ``(source, target)`` gate pairs into concrete wire-cut locations.

    For each pair, any gate sitting between the two endpoints that shares a qubit
    with both becomes the actual cut location; otherwise the original pair is kept.
    """
    refined_cut = []
    for source, target in cut_set:
        source_idx = int(source.split("_")[1])
        target_idx = int(target.split("_")[1])
        _, source_qubits = gates[source_idx - 1]
        _, target_qubits = gates[target_idx - 1]

        found = False
        for idx, (gate_name, gate_qubits) in enumerate(
            gates[source_idx : target_idx - 1], start=source_idx
        ):
            if set(gate_qubits).intersection(source_qubits, target_qubits):
                refined_cut.append((source, f"{gate_name}_{idx + 1}"))
                found = True
        if not found:
            refined_cut.append((source, target))
    return refined_cut


@task(returns=2)
def _find_cut_wire(
    circuit,
    max_qubits=None,
    max_cuts=None,
    max_components=None,
    weights=None,
    seed=None,
    forbidden_gates=None,
    all_candidates=False,
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

    # Refine each cut pair into a concrete wire cut: locate a gate between the
    # two endpoints that shares a qubit with both, and cut just after it.
    gates = iter_gates(circuit)

    # Wire counts per candidate edge, so a cut between gates sharing two qubits is
    # scored as two wire cuts (it severs two wires) rather than one.
    edge_weights = _gate_edge_weights(graph)

    if all_candidates:
        ranked = _rank_candidates(candidates, weights, forbidden_gates, edge_weights)
        return [(_refine_wire_cut(cut, gates), score) for cut, score, _ in ranked], None

    cut_set, best_score, _ = _score_candidates(
        candidates, weights, forbidden=forbidden_gates, edge_weights=edge_weights
    )
    best_cut = _refine_wire_cut(cut_set, gates)
    logger.debug("Best wire cut: %s (score %s)", best_cut, best_score)
    return best_cut, best_score


def _sparsecut(circuit, coupling_map, layout, max_cuts, forbidden):
    """Hardware-aware gate-cut selection ("SparseCut", Guo et al., `arXiv:2511.05492
    <https://arxiv.org/abs/2511.05492>`_).

    Rank the circuit's two-qubit gates by the physical distance between their qubits on
    the QPU **coupling graph** and cut the longest-routing ones first -- removing
    exactly the gates that would otherwise need the most SWAPs (and so cause the most
    crosstalk / error) on a connectivity-limited device. A distance of 1 means the two
    physical qubits are already adjacent (no SWAP, nothing to gain from cutting), so
    only gates that actually route (distance > 1, or unreachable) are candidates.

    :param circuit: The circuit to cut.
    :param coupling_map: QPU connectivity as an iterable of physical-qubit edges
        ``[(0, 1), (1, 2), ...]`` (undirected).
    :param layout: Logical-qubit -> physical-qubit mapping (a dict or a list indexed by
        logical qubit); defaults to the identity mapping.
    :param max_cuts: Cut at most this many gates (the longest-routing); ``None`` cuts
        every routing gate.
    :param forbidden: Gate ids (``"NAME_index"``) that must not be cut.
    :return: ``[(gate_id, distance), ...]`` sorted by distance, longest first.
    """
    graph = networkx.Graph()
    graph.add_edges_from((int(a), int(b)) for a, b in coupling_map)
    distances = dict(networkx.all_pairs_shortest_path_length(graph))

    if layout is None:
        layout = {q: q for q in range(num_qubits(circuit))}
    elif not isinstance(layout, dict):
        layout = {q: p for q, p in enumerate(layout)}

    forbidden = set(forbidden or ())
    ranked = []
    for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
        if len(qubits) != 2:
            continue
        gate_id = f"{name}_{index}"
        if gate_id in forbidden:
            continue
        try:
            distance = distances[layout[qubits[0]]][layout[qubits[1]]]
        except KeyError:
            distance = float("inf")  # unreachable on this coupling map -> cut first
        if distance > 1:
            ranked.append((gate_id, distance))
    ranked.sort(key=lambda item: item[1], reverse=True)
    if max_cuts is not None:
        ranked = ranked[:max_cuts]
    return ranked


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
    top_n=None,
    coupling_map=None,
    layout=None,
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

    top_n : int, optional
        Return the ``top_n`` best candidate cuts (lowest score first) instead of the
        single winner, as a list of ``{"cut": [...], "kind": "gate"|"wire",
        "score": float}`` dicts -- e.g. to compare options and size them with
        :func:`~Qdislib.cut_cost`. ``None`` (default) returns the single best cut.
        Only affects the ``implementation="qdislib"`` path.

    coupling_map : list of (int, int), optional
        The target QPU connectivity, as physical-qubit edges ``[(0, 1), (1, 2), ...]``.
        When given, ``find_cut`` switches to **hardware-aware "SparseCut" mode**
        (Guo et al., `arXiv:2511.05492 <https://arxiv.org/abs/2511.05492>`_): it ranks
        the two-qubit gates by the physical distance between their qubits on the
        coupling graph and returns the longest-routing ones as a **gate cut** -- cutting
        exactly the gates that would otherwise need the most SWAPs (crosstalk / error)
        on a connectivity-limited device. Combine with ``max_cuts`` to bound how many
        gates are cut. This mode ignores ``max_qubits`` / ``wire_cut`` (it is a
        distance-driven gate cut). Only affects the ``implementation="qdislib"`` path.

    layout : dict or list, optional
        Logical-qubit -> physical-qubit mapping used with ``coupling_map`` (a dict, or a
        list indexed by logical qubit). Defaults to the identity mapping.

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

        # Hardware-aware "SparseCut" mode (arXiv:2511.05492): given the QPU coupling
        # map, cut the longest-routing two-qubit gates (a gate cut ranked by physical
        # qubit distance) so the fragments route with the fewest SWAPs / least crosstalk.
        if coupling_map is not None:
            ranked = _sparsecut(circuit, coupling_map, layout, max_cuts, forbidden)
            if top_n is not None:
                # Longest distance first -> best; score = -distance to keep the
                # "lowest score first" convention shared with the other modes.
                return [
                    {"cut": [gate_id], "kind": "gate", "score": -distance}
                    for gate_id, distance in ranked[:top_n]
                ]
            best = [gate_id for gate_id, _ in ranked]
            if not best:
                logger.warning(
                    "find_cut (SparseCut) found no routing two-qubit gate to cut: "
                    "every two-qubit gate is already adjacent on the coupling map, or "
                    "the routing ones are forbidden."
                )
            return best

        # Ranked mode: return the top_n candidate cuts (best score first) across
        # the enabled strategies, each as {"cut", "kind", "score"}, instead of the
        # single winner. Lets the caller compare options (and cost them with
        # cut_cost) rather than trust one pick.
        if top_n is not None:
            ranked = []
            if gate_cut:
                gate_cands, _ = _find_cut_gate(
                    circuit,
                    max_qubits=max_qubits,
                    max_cuts=max_cuts,
                    max_components=max_components,
                    weights=weights,
                    seed=seed,
                    forbidden_gates=forbidden,
                    all_candidates=True,
                )
                for cut, score in compss_wait_on(gate_cands):
                    ranked.append({"cut": cut, "kind": "gate", "score": score})
            if wire_cut:
                wire_cands, _ = _find_cut_wire(
                    circuit,
                    max_qubits=max_qubits,
                    max_cuts=max_cuts,
                    max_components=max_components,
                    weights=weights,
                    seed=seed,
                    forbidden_gates=forbidden,
                    all_candidates=True,
                )
                for cut, score in compss_wait_on(wire_cands):
                    ranked.append({"cut": cut, "kind": "wire", "score": score})
            ranked.sort(key=lambda candidate: candidate["score"])
            return ranked[:top_n]

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
