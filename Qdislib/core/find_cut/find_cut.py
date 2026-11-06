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

import itertools
import logging

logger = logging.getLogger(__name__)

__all__ = ["find_cut", "draw_cut"]

import networkx

from Qdislib.core.find_cut.draw import draw_cut  # noqa: F401  (public API)
import typing

from Qdislib.utils.runtime import task, compss_wait_on

import pymetis

from Qdislib.utils.backends import convert_to_dag, iter_gates, num_qubits


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
        if max_cuts is not None and cutset_size > max_cuts:
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

    if max_cuts is not None and cut_wires > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return partition_list, cut_edges


_warned_no_sklearn = False


def _find_best_spectral_clustering(graph, clusters, max_qubits, max_cuts, seed=None):
    global _warned_no_sklearn
    try:
        from sklearn.cluster import SpectralClustering
    except ImportError as exc:
        # scikit-learn missing, or its bundled OpenMP runtime failing to load
        # (glibc < 2.32: "cannot allocate memory in static TLS block"). The other
        # strategies still run, so find_cut degrades instead of crashing.
        if not _warned_no_sklearn:
            logger.warning(
                "find_cut: spectral clustering skipped, scikit-learn could not be "
                "imported (%s). If this is a static TLS error, start Python with "
                "LD_PRELOAD pointing at scikit_learn.libs/libgomp*.so.",
                exc,
            )
            _warned_no_sklearn = True
        return [], []

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
    cut_edges, cut_size = [], 0
    for u, v, gate in graph.edges(data=True):
        if labels[u] != labels[v]:
            cut_edges.append(gate["gate"])
            cut_size += gate.get("weight", 1)  # wires a gate-graph edge severs

    if max_cuts is not None and cut_size > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return clusters, cut_edges


def _metis_partition(graph, num_clusters, max_nodes_per_cluster, max_cuts):
    """METIS k-way partition minimising the *weighted* edge cut.

    The graphs here are multigraphs: two gates on the same qubit pair are two
    edges, and a gate-graph edge carries the number of wires it severs as its
    ``weight``. METIS takes a simple graph, so each neighbour pair becomes one edge
    whose weight is the sum over its parallel edges -- otherwise two gates between
    the same qubits would cost METIS the same as one.
    """
    num_clusters = min(num_clusters, graph.number_of_nodes())
    if num_clusters < 2:
        return [], []

    nodes = list(graph.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    weights = {}
    for u, v, data in graph.edges(data=True):
        if u == v:
            continue
        key = (min(index[u], index[v]), max(index[u], index[v]))
        weights[key] = weights.get(key, 0) + int(data.get("weight", 1))
    neighbours = [[] for _ in nodes]
    for (a, b), w in weights.items():
        neighbours[a].append((b, w))
        neighbours[b].append((a, w))
    xadj, adjncy, eweights = [0], [], []
    for row in neighbours:
        for b, w in sorted(row):
            adjncy.append(b)
            eweights.append(w)
        xadj.append(len(adjncy))
    if not adjncy:
        return [], []
    # pymetis >= 2025 takes a CSRAdjacency; plain xadj/adjncy are going away in 2027
    if hasattr(pymetis, "CSRAdjacency"):
        graph_in = {"adjacency": pymetis.CSRAdjacency(xadj, adjncy)}
    else:
        graph_in = {"xadj": xadj, "adjncy": adjncy}
    _, parts = pymetis.part_graph(num_clusters, eweights=eweights, **graph_in)
    clusters = {}
    for node, part in zip(nodes, parts):
        clusters.setdefault(part, []).append(node)

    for cluster_id, members in clusters.items():
        if max_nodes_per_cluster and len(members) > max_nodes_per_cluster:
            logger.info(f"Cluster {cluster_id} exceeds the max size.")
            return [], []

    clusters = list(clusters.values())
    part_of = dict(zip(nodes, parts))
    cut_edges, cut_size = [], 0
    for u, v, data in graph.edges(data=True):
        if part_of[u] != part_of[v]:
            cut_edges.append(data["gate"])
            cut_size += data.get("weight", 1)
    if max_cuts is not None and cut_size > max_cuts:
        return [], []
    elif not cut_edges:
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


def _wire_cut_fits(dag, cut, max_qubits, max_cuts, forbidden_qubits):
    """Whether wire cut ``cut`` is usable and within the user's limits.

    Checks the cut as ``wire_cutting`` will apply it, on a copy of ``dag``: no cut
    on a ``forbidden_qubits`` wire, at most ``max_cuts`` single-wire cuts, an actual
    split into more pieces, and no piece wider than ``max_qubits`` qubits. (The
    gate graph's partitions hold gates, so their size says nothing about qubits.)
    """
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expand_wire_cuts,
        _insert_wire_cut_markers,
    )
    from Qdislib.utils.graph_qibo import _max_qubit, _remove_red_edges

    try:
        wire_cuts = _expand_wire_cuts(dag, cut)
    except KeyError:  # an endpoint the DAG doesn't have
        return False
    if not wire_cuts:
        return False
    if forbidden_qubits and set(forbidden_qubits) & {q for _, _, q in wire_cuts}:
        return False
    if max_cuts is not None and len(wire_cuts) > max_cuts:
        return False
    marked = _insert_wire_cut_markers(dag.copy(), wire_cuts, _max_qubit(dag))
    split = _remove_red_edges(marked).to_undirected()
    pieces = list(networkx.connected_components(split))
    if len(pieces) <= networkx.number_connected_components(dag.to_undirected()):
        return False
    if max_qubits is not None:
        for piece in pieces:
            width = {q for node in piece for q in marked.nodes[node].get("qubits", ())}
            if len(width) > max_qubits:
                return False
    return True


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
    forbidden_qubits=None,
):
    """Find the best *wire* cut for ``circuit`` by partitioning its gate graph.

    Builds the gate-dependency graph, runs the four partitioning strategies and
    keeps the lowest-scoring cut. The winning cut-set (pairs of gate identifiers)
    is then refined into concrete wire cuts: for each pair, any gate sitting
    between the two endpoints that shares a qubit with both becomes the actual
    cut location. Backend-agnostic (Qibo / Qiskit / CUDA-Q).

    Every candidate is checked as it will be applied (:func:`_wire_cut_fits`):
    ``max_qubits`` bounds each resulting piece's qubits and ``forbidden_qubits``
    the wires that may be cut. The partitioners themselves get no ``max_qubits``:
    their parts are sets of gates, so a gate count is no qubit bound.

    :return: ``(best_cut, best_score)`` where ``best_cut`` is a list of
        ``(source_gate, target_gate)`` identifier pairs.
    """
    logger.debug("Finding best wire cut")

    graph = _circuit_to_gate_graph(circuit)
    # The partitioners never see forbidden_qubits; make severing a protected wire
    # very expensive in the graph they partition (Kernighan-Lin and spectral
    # clustering weigh edges) so they look for splits elsewhere. Scoring and the
    # final check still use the real graph.
    part_graph = graph
    if forbidden_qubits:
        protected = set(forbidden_qubits)
        part_graph = graph.copy()
        for u, v, data in part_graph.edges(data=True):
            shared = set(graph.nodes[u]["qubits"]) & set(graph.nodes[v]["qubits"])
            if shared & protected:
                data["weight"] = data.get("weight", 1) * 1000
    candidates = {
        "Kernighan-Lin": _find_best_kernighan_lin_partition(
            graph=part_graph, max_qubits=None, max_cuts=None, seed=seed
        ),
        "Girvan-Newman": _find_best_girvan_newman_cutset(
            part_graph, max_qubits=None, max_cuts=max_cuts
        ),
        "Spectral clustering": _find_best_spectral_clustering(
            part_graph,
            clusters=max_components,
            max_qubits=None,
            max_cuts=max_cuts,
            seed=seed,
        ),
        "METIS": _metis_partition(
            part_graph,
            num_clusters=max_components,
            max_nodes_per_cluster=None,
            max_cuts=max_cuts,
        ),
    }

    # Refine each cut pair into a concrete wire cut: locate a gate between the
    # two endpoints that shares a qubit with both, and cut just after it.
    gates = iter_gates(circuit)

    # Wire counts per candidate edge, so a cut between gates sharing two qubits is
    # scored as two wire cuts (it severs two wires) rather than one.
    edge_weights = _gate_edge_weights(graph)

    dag = convert_to_dag(circuit)

    def fitting(candidates):
        ranked = []
        for cut, score, _ in _rank_candidates(
            candidates, weights, forbidden_gates, edge_weights
        ):
            refined = _refine_wire_cut(cut, gates)
            if _wire_cut_fits(dag, refined, max_qubits, max_cuts, forbidden_qubits):
                ranked.append((refined, score))
            else:
                logger.debug("Wire cut %s: skipped (breaks a limit)", refined)
        return ranked

    ranked = fitting(candidates)
    # A qubit bound or a protected wire can need more pieces than the target count
    # derived from max_qubits (a chain of n qubits needs about n / (max_qubits - 1)).
    # Ask the clustering strategies for more pieces until one fits.
    clusters = max_components or 2
    most = min(graph.number_of_nodes(), clusters + 8)
    while not ranked and clusters < most:
        clusters += 1
        ranked = fitting(
            {
                "Spectral clustering": _find_best_spectral_clustering(
                    part_graph,
                    clusters=clusters,
                    max_qubits=None,
                    max_cuts=max_cuts,
                    seed=seed,
                ),
                "METIS": _metis_partition(
                    part_graph,
                    num_clusters=clusters,
                    max_nodes_per_cluster=None,
                    max_cuts=max_cuts,
                ),
            }
        )

    if all_candidates:
        return ranked, None

    best_cut, best_score = ranked[0] if ranked else ([], float("inf"))
    logger.debug("Best wire cut: %s (score %s)", best_cut, best_score)
    return best_cut, best_score


def _routing_distances(circuit, coupling_map, layout):
    """``{gate_id: physical distance}`` of every two-qubit gate on the coupling map
    (``inf`` when its qubits are not connected at all)."""
    graph = networkx.Graph()
    graph.add_edges_from((int(a), int(b)) for a, b in coupling_map)
    distances = dict(networkx.all_pairs_shortest_path_length(graph))
    if layout is None:
        layout = {q: q for q in range(num_qubits(circuit))}
    elif not isinstance(layout, dict):
        layout = {q: p for q, p in enumerate(layout)}
    routing = {}
    for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
        if len(qubits) != 2:
            continue
        try:
            distance = distances[layout[qubits[0]]][layout[qubits[1]]]
        except KeyError:
            distance = float("inf")
        routing[f"{name}_{index}"] = distance
    return routing


def _sparsecut(
    circuit, coupling_map, layout, max_qubits, max_cuts, forbidden, seed, time_limit
):
    """Hardware-aware gate cut ("SparseCut", Guo et al., `arXiv:2511.05492
    <https://arxiv.org/abs/2511.05492>`_): cut the gates that would need the most
    SWAPs on the device.

    Only gates whose two qubits are not adjacent on the coupling map are cut, and a
    gate's cost is divided by its routing distance, so the cheapest cut uses the
    longest-routing gates. Unlike a plain ranking, the cut returned separates the
    circuit (every cut separates its two sides), within ``max_qubits`` when given.

    :return: ``(cut, distances)``: gate ids, longest-routing first, and the routing
        distance of every two-qubit gate.
    """
    from Qdislib.core.find_cut.cost_search import cost_aware_cut

    distances = _routing_distances(circuit, coupling_map, layout)
    adjacent = {gid for gid, d in distances.items() if d <= 1}
    scale = {
        gid: (0.0 if d == float("inf") else 1.0 / d)
        for gid, d in distances.items()
        if d > 1
    }
    cuts, _ = cost_aware_cut(
        circuit,
        max_qubits,
        max_cuts,
        hadamard=False,
        wires=False,
        joint=False,
        time_limit=time_limit,
        seed=seed,
        forbidden_gates=set(forbidden or ()) | adjacent,
        balanced=False,
        gate_scale=scale,
    )
    cut = sorted(cuts["gate"], key=lambda gid: -distances[gid])
    return cut, distances


def separable_components(circuit, max_qubits=None):
    """Groups of qubits the circuit already splits into, needing no cuts.

    Every strategy :func:`find_cut` runs -- Kernighan-Lin, Girvan-Newman,
    spectral, METIS -- searches for the *cheapest cut*. None of them asks
    whether a cut is needed at all, and a circuit whose qubit interaction graph
    is already disconnected does not need one: the state is a product across
    that boundary, so the fragments are independent circuits whose results
    multiply. Cutting one anyway costs ``gamma = 3`` and nine times the shots
    for nothing.

    It is not a corner case. Bernstein-Vazirani only entangles a qubit where its
    secret bit is 1, so a 24-qubit instance leaves twelve qubits touched by no
    two-qubit gate at all, and splits 12/12 for free -- where ``find_cut`` used
    to propose a gate cut and hand back lopsided 1q and 11q fragments.

    Components are packed into as few groups as possible, largest first, so the
    answer is a partition that respects ``max_qubits`` rather than merely a list
    of components.

    :param max_qubits: Qubit budget per fragment. ``None`` means only that the
        circuit must split into at least two parts, so any grouping will do.
    :return: List of qubit-index lists, one per fragment, or ``None`` when the
        circuit cannot be separated without cutting.
    """
    graph = networkx.Graph()
    total = num_qubits(circuit)
    graph.add_nodes_from(range(total))
    for _, qubits in iter_gates(circuit):
        for left, right in itertools.combinations(sorted(qubits), 2):
            graph.add_edge(left, right)

    components = sorted(
        (sorted(c) for c in networkx.connected_components(graph)), key=len, reverse=True
    )
    if len(components) < 2:
        return None
    if max_qubits is None:
        # Any split will do, so take the cheapest balanced one: largest
        # component alone against everything else.
        return [components[0], sorted(q for c in components[1:] for q in c)]
    if len(components[0]) > max_qubits:
        # One component is already too wide; it has to be cut, and that is the
        # normal search's job.
        return None

    # First-fit-decreasing: enough fragments that each fits the budget, and no
    # more.
    bins: typing.List[typing.List[int]] = []
    for component in components:
        for group in bins:
            if len(group) + len(component) <= max_qubits:
                group.extend(component)
                break
        else:
            bins.append(list(component))
    if len(bins) < 2:
        return None
    return [sorted(group) for group in bins]


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
    mixed=False,
    hadamard=True,
    joint=True,
    exact="auto",
    time_limit=10.0,
):
    """
    Automatically partitions a quantum circuit to optimize for parallel execution
    and hardware constraints.

    The FindCut algorithm from the Qdislib library is designed to find the best cut to
    decompose a quantum circuit into subcircuits that adhere to specific
    constraints such as the number of qubits per subcircuit, the number of allowed cuts,
    and the number of resulting components. It supports both wire cutting and gate cutting
    strategies, which are evaluated through multiple graph partitioning techniques
    (Kernighan-Lin, Girvan-Newman, Spectral Decomposition, METIS). It scores all possible options
    with a loss function and then decides which cut is better for that specific circuit and constraints.

    Parameters
    ----------
    circuit : a circuit of any SDK
        The input quantum circuit to be partitioned: Qiskit, qibo, CUDA-Q, PennyLane
        or QiliSDK.

    max_qubits : int, optional
        Maximum number of qubits allowed per subcircuit. Helps enforce hardware
        constraints (e.g., maximum qubits per QPU). When given, it also determines
        how many subcircuits the clustering strategies aim for
        (``ceil(num_qubits / max_qubits)``); otherwise a simple bisection is used.

    max_cuts : int, optional
        Maximum number of cuts permitted across the circuit. Controls the complexity
        and computational overhead of the partitioning. ``0`` allows no cut, so the
        result is ``[]`` (``{"gate": [], "wire": []}`` with ``mixed=True``).

    wire_cut : bool, default=True
        Whether to consider wire cuts (cutting dependencies between gates)
        during partitioning.

    gate_cut : bool, default=True
        Whether to consider gate cuts (cutting the circuit by splitting two-qubit gates)
        during partitioning.

    implementation : str, default="qdislib"
        Kept for compatibility; ``"qdislib"`` is the only value. (``"ibm-ckt"``,
        delegating to qiskit-addon-cutting, was removed in 2.0.0: it was slow on
        larger circuits and dropped the wire cuts it found.)

    weights : dict, optional
        Overrides for the score of a candidate cut,
        ``|w1 * cuts + w2 * fragments + w3 * largest fragment|``, the lowest wins:
        ``{"w1": 100, "w2": -1, "w3": 1}`` by default. Ignored with ``mixed=True`` and
        ``coupling_map``.

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
        :func:`~Qdislib.cutting_overhead`. ``None`` (default) returns the single best
        cut.
        Only affects the ``implementation="qdislib"`` path.

    coupling_map : list of (int, int), optional
        The target QPU connectivity, as physical-qubit edges ``[(0, 1), (1, 2), ...]``.
        When given, ``find_cut`` switches to **hardware-aware "SparseCut" mode**
        (Guo et al., `arXiv:2511.05492 <https://arxiv.org/abs/2511.05492>`_): it cuts
        only gates whose qubits are not adjacent on the device, preferring the
        longest-routing ones (the gates that would need the most SWAPs), and returns
        the cheapest such gate cut that **separates** the circuit, within
        ``max_qubits`` when given and ``max_cuts``; ``[]`` if no such cut exists.

    layout : dict or list, optional
        Logical-qubit -> physical-qubit mapping used with ``coupling_map`` (a dict, or a
        list indexed by logical qubit). Defaults to the identity mapping.

    mixed : bool, default=False
        Search every kind of cut at once and price it the way
        :func:`~Qdislib.circuit_cutting` runs it: each gate cut by the cheaper of
        measure-and-prepare and the Hadamard test (angle-aware), the gates across one
        boundary cut jointly, wire cuts where a qubit's timeline is better split. The
        objective is the sampling overhead gamma; every fragment, its Hadamard ancillas
        included, fits ``max_qubits`` (``None``: halves), and every cut separates.
        Returns ``{"gate": [...], "hadamard": [...], "wire": [...]}`` for
        :func:`~Qdislib.circuit_cutting`. ``gate_cut`` / ``wire_cut`` switch kinds off.

    hadamard : bool, default=True
        With ``mixed=True``, allow Hadamard-test gate cuts.

    joint : bool, default=True
        With ``mixed=True``, price the Hadamard cuts across one boundary jointly (as
        ``circuit_cutting(..., joint=True)`` runs them).

    exact : bool or "auto", default="auto"
        With ``mixed=True`` (and SparseCut), seed the search with the exact optimal
        gate partition (a MILP, SciPy/HiGHS); ``"auto"`` does so up to 40 qubits.

    time_limit : float, default=10
        Seconds the ``mixed=True`` / SparseCut search may take; it returns the best
        cut found by then.

    Diagnostic details about the partitioning process are logged at DEBUG level;
    enable them with ``qd.set_log_level("DEBUG")``.

    Returns
    -------
    best_cut : list or dict
        The best cut found: gate names (e.g. ``["CX_3"]``) or ``(source, target)`` wire
        edges, ``[]`` when none is needed or none fits. With ``mixed=True``
        ``{"gate": [...], "hadamard": [...], "wire": [...]}``. With ``top_n``, a list of
        ``{"cut", "kind", "score"}`` dicts, best first.

    Raises
    ------
    ValueError
        If ``max_qubits`` is less than 1, ``max_cuts`` is negative, ``implementation``
        is unknown, or neither ``wire_cut`` nor ``gate_cut`` is enabled.

    Notes
    -----
    - Without ``mixed`` or ``coupling_map``, circuits are transformed into graphs and
      evaluated through multiple graph partitioning strategies; the best partition is
      chosen based on a custom loss function.
    - A circuit that already splits into independent parts within ``max_qubits`` needs
      no cut, and gets none.

    Example
    -------
    >>> import Qdislib as qd
    >>> cut = qd.find_cut(my_circuit, max_qubits=5, gate_cut=True, wire_cut=False)
    >>> print(cut)  # the best cut, e.g. ['CX_3']
    """
    # 0 is a budget, not "no limit": check against None, never truthiness.
    if max_qubits is not None and max_qubits < 1:
        raise ValueError(f"max_qubits must be at least 1, got {max_qubits}")
    if max_cuts is not None and max_cuts < 0:
        raise ValueError(f"max_cuts must be 0 or more, got {max_cuts}")

    if implementation == "ibm-ckt":
        raise ValueError(
            "implementation='ibm-ckt' was removed in 2.0.0 (it was slow on larger "
            "circuits and dropped its wire cuts); use find_cut(..., mixed=True), which "
            "searches gate, Hadamard and wire cuts together."
        )
    if implementation != "qdislib":
        raise ValueError(f"Unknown implementation={implementation!r}; use 'qdislib'.")

    # Before searching for a cut, check whether one is needed. A circuit
    # whose interaction graph is already disconnected splits for free, and
    # every strategy below would instead return a cut costing gamma = 3.
    separable = separable_components(circuit, max_qubits)
    if separable is not None:
        logger.info(
            "find_cut: the circuit already separates into %d fragments of "
            "%s qubits with no cut at all; returning none.",
            len(separable),
            [len(g) for g in separable],
        )
        if mixed:
            return {"gate": [], "hadamard": [], "wire": []}
        if top_n is not None:
            return [{"cut": [], "kind": "none", "score": 0.0}]
        return []

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

    # Mixed mode: one cost-weighted partition of the unified qubit-timeline graph,
    # classifying each severed connection as a gate cut or a wire cut -> a combined
    # {"gate": [...], "wire": [...]} for gate_wire_cutting.
    if mixed:
        from Qdislib.core.find_cut.cost_search import cost_aware_cut

        cuts, _ = cost_aware_cut(
            circuit,
            max_qubits,
            max_cuts,
            gates=gate_cut,
            hadamard=hadamard and gate_cut,
            joint=joint,
            wires=wire_cut,
            exact=exact,
            time_limit=time_limit,
            seed=seed,
            forbidden_gates=forbidden,
            forbidden_qubits=forbidden_qubits,
        )
        return cuts

    # Hardware-aware "SparseCut" mode (arXiv:2511.05492): given the QPU coupling
    # map, cut the longest-routing two-qubit gates (a gate cut ranked by physical
    # qubit distance) so the fragments route with the fewest SWAPs / least crosstalk.
    if coupling_map is not None:
        best, distances = _sparsecut(
            circuit,
            coupling_map,
            layout,
            max_qubits,
            max_cuts,
            forbidden,
            seed,
            time_limit,
        )
        if not best and max_cuts != 0:
            logger.warning(
                "find_cut (SparseCut) found no separating cut made of routing "
                "gates: the non-adjacent two-qubit gates (not forbidden) do not "
                "split the circuit within the limits."
            )
        if top_n is not None:
            # One separating cut; score = -total routing distance, so the
            # "lowest score first" convention of the other modes holds.
            if not best:
                return []
            score = -sum(distances[g] for g in best)
            return [{"cut": best, "kind": "gate", "score": score}][:top_n]
        return best

    if not (gate_cut or wire_cut):
        raise ValueError("Enable at least one of wire_cut or gate_cut.")

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
                forbidden_qubits=forbidden_qubits,
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
            forbidden_qubits=forbidden_qubits,
        )

    # Waiting only now: under PyCOMPSs the two searches run side by side.
    if gate_cut:
        best, best_score = compss_wait_on(cut_gate), compss_wait_on(score_gate)
    if wire_cut:
        cut_wire, score_wire = compss_wait_on(cut_wire), compss_wait_on(score_wire)
        if not gate_cut or score_wire < best_score:
            best = cut_wire

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
