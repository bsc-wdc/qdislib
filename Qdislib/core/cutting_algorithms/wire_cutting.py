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

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "wire_cutting",
    "wire_cutting_subcircuits",
    "wire_cutting_subcircuit_reconstruction",
]

import math
import networkx
import numpy as np
import typing

import re

try:
    import qibo
    from qibo import gates

    qibo_available = True
except ImportError:  # qibo is optional; functions that need it fail when called
    qibo = gates = None
    qibo_available = False

try:
    import cudaq

    cudaq_available = True
except ImportError:
    cudaq = None
    cudaq_available = False

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN

from Qdislib.utils.graph_qibo import (
    _max_qubit,
    _update_qubits,
    _remove_red_edges,
)
from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel


from Qdislib.utils.backends import (
    convert_to_dag,
    change_basis_dag,
    normalize_observables,
    make_subcircuit,
    CutSubcircuits,
)
from Qdislib.utils.exceptions import UnsupportedBackendError, InvalidCutError


def wire_cutting(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    sync: bool = True,
    observables: str = None,
    shots: int = 1024,
    backend: str = "numpy",
    qpu: str = None,
    software: str = "qibo",
    cache=None,
    cache_mode="expectation",
):
    """Apply wire cutting to a quantum circuit to simplify its structure for distributed or hardware-constrained execution of subcircuits.

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
    software : str, default="qibo"
        Backend used to build and execute the subcircuits: "qibo", "qiskit" or "cudaq".
    cache : CircuitCache, optional
        Enable the semantic subcircuit cache (see
        :func:`~Qdislib.utils.circuit_cache.default_circuit_cache`). When given,
        evaluation runs the CPU-simulation subcircuit workflow (GPU/QPU dispatch
        is bypassed). ``None`` (default) disables caching.
    cache_mode : str or CacheStrategy, default="expectation"
        Cache granularity; ``"expectation"``, ``"counts"``, or a custom strategy.

    Returns
    -------
    final_recons : float
        The reconstructed expectation value of the original quantum circuit after wire cutting.

    Notes
    -----
    - If the input circuit has multiple disconnected components, they are processed independently.
    - Each cut doubles the number of circuit evaluations required.
    - The final expectation value is scaled by 1 / 2^(number of cuts).
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
    # A weighted sum of Pauli observables (a Hamiltonian) is evaluated term by
    # term: <sum_i c_i P_i> = sum_i c_i <P_i>. Each term reuses the single
    # Pauli-string path below. A plain string (or None) skips this and runs once.
    if observables is not None and not isinstance(observables, str):
        total = 0.0
        for coeff, pauli in normalize_observables(observables):
            base = rand_qc.copy() if isinstance(rand_qc, networkx.Graph) else rand_qc
            total += coeff * wire_cutting(
                base,
                cut,
                sync=sync,
                observables=pauli,
                shots=shots,
                backend=backend,
                qpu=qpu,
                software=software,
                cache=cache,
                cache_mode=cache_mode,
            )
        return total

    # Semantic cache: reuse the recommended subcircuit workflow, which routes
    # each subcircuit through expectation_value(cache=...). This runs the CPU
    # simulation path, so GPU/QPU dispatch is bypassed when caching is enabled.
    if cache is not None:
        from Qdislib.utils.backends import expectation_value

        subcircuits = wire_cutting_subcircuits(
            rand_qc, cut, observables=observables, software=software
        )
        values = [
            expectation_value(s, shots=shots, cache=cache, cache_mode=cache_mode)
            for s in subcircuits
        ]
        return subcircuits.reconstruct(values)

    dag = convert_to_dag(rand_qc, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    if networkx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected())
        ]
        results = []
        total_wire_cuts = 0
        for s in S:
            num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in cut:
                if s.has_edge(*c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                cut_wires = _expand_wire_cuts(s, tmp_cuts)
                total_wire_cuts += len(cut_wires)
                graphs = _generate_wire_cutting(
                    s, cut_wires, num_qubits=num_qubits, qpu=qpu, software=software
                )
                graphs = _sum_results(graphs)
                results.append(graphs)
            else:
                results.append(
                    _build_and_run_subcircuit(s, software, shots, backend, qpu)
                )
        if sync:
            results = compss_wait_on(results)
        final_recons = 1 / (2**total_wire_cuts) * math.prod(results)
        return final_recons
    else:
        num_qubits = _max_qubit(dag)
        if cut:
            cut_wires = _expand_wire_cuts(dag, cut)
            results = _generate_wire_cutting(
                dag, cut_wires, num_qubits=num_qubits, qpu=qpu, software=software
            )

            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(cut_wires)) * sum(results)
        else:
            final_recons = _build_and_run_subcircuit(dag, software, shots, backend, qpu)
        return final_recons


def _expand_wire_cuts(dag, edges_to_replace):
    """Expand each ``(source, target)`` cut into one entry per shared qubit.

    A wire cut severs a single qubit line. When the two gates of a cut share more
    than one qubit, separating them requires cutting *each* shared wire, so this
    expands the pair into one ``(source, target, qubit)`` single-wire cut per
    shared qubit (the common case of a single shared qubit yields one entry,
    unchanged).

    :param dag: The circuit DAG (nodes carry a ``qubits`` attribute).
    :param edges_to_replace: List of ``(source, target)`` gate-name pairs.
    :return: List of ``(source, target, qubit)`` single-wire cuts.
    """
    wire_cuts = []
    for source, target in edges_to_replace:
        shared = sorted(
            set(dag.nodes[source].get("qubits", ()))
            & set(dag.nodes[target].get("qubits", ()))
        )
        for qubit in shared:
            wire_cuts.append((source, target, qubit))
    return wire_cuts


@task(returns=1)
def _build_and_run_subcircuit(component, software, shots, backend, qpu):
    """Build one subcircuit in the chosen backend and return its expectation value.

    Made a single task so that, under PyCOMPSs, the native circuit is built **and**
    executed in one worker process: the picklable component DAG crosses the task
    boundary, never the fragile native circuit (Qibo collapse symbol / CUDA-Q
    kernel). The inner ``dag_to_circuit_*`` / ``_expec_value_*`` calls run inline.

    Backend-agnostic execution dispatch: the subcircuit is built with the
    matching ``dag_to_circuit_*`` and evaluated with the matching
    ``_expec_value_*``. A QPU target ("MN_Ona") overrides the software choice.

    :param component: A connected-component DAG of the (already cut) circuit.
    :param software: Execution backend: "qibo", "qiskit" or "cudaq".
    :param shots: Number of measurement shots.
    :param backend: Qibo simulation backend (only used when ``software == "qibo"``).
    :param qpu: Optional QPU target (e.g. "MN_Ona").
    :return: The subcircuit's expectation value (a PyCOMPSs future under COMPSs).
    """
    s_new, highest_qubit = _update_qubits(component)

    # Backend-agnostic execution: the selected backend evaluates the subcircuit.
    # Wire cutting has no measure-and-prepare sign, so there are no prep_qubits;
    # ``backend`` is the Qibo platform when software == "qibo".
    from Qdislib.core.backend.registry import build_execution_backend, build_backend
    from Qdislib.core.backend.base_backend import BackendConfig

    if qpu == "MN_Ona":
        qpu_backend = build_backend("BSC_QPU_BACKEND", BackendConfig(num_qubits=highest_qubit))
        return qpu_backend.evaluate_expectation(
            s_new, highest_qubit, (), shots, backend or "numpy"
        )

    extra = {"platform": backend} if software == "qibo" and backend else {}
    exec_backend = build_execution_backend(software, num_qubits=highest_qubit, **extra)
    return exec_backend.evaluate_expectation(s_new, highest_qubit, (), shots)


def _generate_wire_cutting(
    dag: typing.Any,
    wire_cuts: typing.List[typing.Tuple[str, str, int]],
    num_qubits: int,
    shots: int = 10000,
    backend: str = "numpy",
    qpu=None,
    software: str = "qibo",
) -> typing.List[float]:
    """Apply wire cuts to the DAG and evaluate the reconstruction terms.

    Each wire cut ``(source, target, qubit)`` replaces the wire carrying ``qubit``
    from ``source`` to ``target`` with a measure node ``O_i`` (on the original
    qubit) and a prepare node ``PS_i`` (on a fresh qubit line), so the two sides
    become independent. The wire-cutting identity is then expanded over the
    ``8`` measure/prepare basis terms per cut, giving ``8 ** n_cuts`` sign-adjusted
    products to sum.

    :param dag: The circuit DAG to modify (mutated in place).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples
        (see :func:`_expand_wire_cuts`).
    :param num_qubits: Highest qubit index in the circuit; fresh qubit lines are
        numbered above it.
    :return: List of ``8 ** len(wire_cuts)`` reconstruction terms.
    """
    for index, (source, target, qubit) in enumerate(wire_cuts, start=1):
        # The shared edge is removed once; further wires between the same pair
        # are already disconnected.
        if dag.has_edge(source, target):
            dag.remove_edge(source, target)

        # Relabel this qubit on the target and its successors onto a fresh line so
        # the prepared (PS) side is disconnected from the measured (O) side.
        successors = [node for node in dag.nodes if dag.has_edge(target, node)]
        for node in [target] + successors:
            node_qubits = dag.nodes[node].get("qubits")
            if node_qubits and qubit in node_qubits:
                dag.nodes[node]["qubits"] = tuple(
                    num_qubits + index if q == qubit else q for q in node_qubits
                )

        dag.add_node(f"O_{index}", gate="S", qubits=(qubit,), parameters=())
        dag.add_node(
            f"PS_{index}", gate="T", qubits=(num_qubits + index,), parameters=()
        )
        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")

    reconstruction = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag.copy())

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )
        graph_components = [networkx.DiGraph().copy() for _ in range(num_components)]

        _generate_subcircuits_wire_cutting(
            copy_graph,
            num_qubits + len(wire_cuts),
            index,
            wire_cuts,
            graph_components,
        )

        exp_value = [
            _build_and_run_subcircuit(s, software, shots, backend, qpu)
            for s in graph_components
        ]
        exp_value = _change_sign(exp_value, index)
        reconstruction.append(exp_value)

    return reconstruction


# Plain function (not a @task): the wire-cut substitution is cheap graph surgery,
# run on the master, populating graph_components in place with concrete DiGraphs.
# Keeping it task-free (no COLLECTION_OUT futures) lets the generated subcircuits
# be wrapped as serializable Subcircuit objects under PyCOMPSs.
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
            updated_dag.remove_node(f"PS_{idx2}")

        # I 1
        elif index == 1:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "Observable I"
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
            raise ValueError(f"Unexpected wire-cut substitution index: {index}")

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
def _counts_qibo(qc, shots=1024, backend="numpy"):
    """Simulate a Qibo circuit and return its measurement counts (frequencies).

    Appends a measurement on every qubit in place -- so the circuit is single-use
    (evaluate once); deep-copying is not viable here as it severs the symbol
    binding that gate-cut collapse measurements rely on -- then runs it.

    :param qc: A Qibo ``Circuit``.
    :param shots: Number of measurement shots.
    :param backend: Qibo backend name.
    :return: A ``{bitstring: count}`` dict (qubit-0-first keys).
    """
    qibo.set_backend(backend)
    qc.add(gates.M(*range(qc.nqubits)))
    result = qc(nshots=shots)
    return dict(result.frequencies(binary=True))


def _expectation_from_counts_qibo(freq, observables, shots=1024):
    """Compute ``<O>`` for a Z/I Pauli string from Qibo measurement counts.

    Qibo frequency keys are already qubit-0-first bit strings, so they align
    directly with the observable string.

    :param freq: A ``{bitstring: count}`` dict as returned by :func:`_counts_qibo`.
    :param observables: Per-qubit Pauli string over ``{I, Z}``.
    :param shots: Number of shots the counts were drawn from.
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
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


def _expectation_from_samples_qibo(qc, observables, shots, backend, prep_qubits):
    """Decoupled ``<O>`` for gate-cut subcircuits with identity on a cut qubit.

    Mirrors :func:`~Qdislib.core.cutting_algorithms.gate_cutting._expectation_from_counts_qiskit`
    for Qibo. A cut qubit measured for the decomposition (``gates.M(collapse=True)``)
    carries a sign normally read off its forced final measurement. When the user
    observes it with identity, we recover the sign from the *cut-point* outcome
    ``m`` (the collapse gate's per-shot samples) as ``-(-1)**m`` instead. This
    needs per-shot samples -- not aggregated frequencies -- because ``m`` can be
    correlated with the other qubits' outcomes.

    :param qc: The subcircuit (its collapse ``M`` gates hold the cut-point outcomes).
    :param observables: Per-qubit Pauli string over ``{I, Z}`` (qubit-0-first).
    :param shots: Number of measurement shots.
    :param backend: Qibo backend name.
    :param prep_qubits: Renumbered qubits measured for the cut.
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    qibo.set_backend(backend)
    collapse = {
        gate.qubits[0]: gate
        for gate in qc.queue
        if isinstance(gate, gates.M) and getattr(gate, "collapse", False)
    }
    qc.add(gates.M(*range(qc.nqubits)))
    final = np.asarray(qc(nshots=shots).samples())  # (shots, nqubits), qubit-0-first
    contribution = np.ones(len(final))
    for q, obs in enumerate(observables):
        if q in prep_qubits and obs == "I" and q in collapse:
            m = np.asarray(collapse[q].result.samples()).ravel()
            contribution *= -((-1.0) ** m)
        elif obs == "Z":
            contribution *= (-1.0) ** final[:, q]
        elif obs == "I":
            continue
        else:
            raise ValueError(f"Unsupported observable {obs}")
    return float(contribution.mean())


def _expec_value_qibo(
    subcirc: typing.Any, shots=1024, backend="numpy", prep_qubits=None
) -> float:
    """Execute the given circuit.

    :param subcirc: Circuit to execute.
    :param prep_qubits: Renumbered cut qubits measured for a gate cut; when any is
        observed with identity, its decomposition sign is recovered from the
        cut-point measurement (see :func:`_expectation_from_samples_qibo`).
    :raises ValueError: If there is an unsupported observable.
    :return: The circuit expected value.
    """
    if subcirc is None:
        return None

    if type(subcirc) != qibo.models.circuit.Circuit:
        observables = subcirc[1]
        qc = subcirc[0]
    else:
        observables = "Z" * subcirc.nqubits
        qc = subcirc

    if prep_qubits:
        return _expectation_from_samples_qibo(
            qc, observables, shots, backend, set(prep_qubits)
        )

    freq = _counts_qibo(qc, shots, backend)
    return _expectation_from_counts_qibo(freq, observables, shots)


def _counts_cudaq(kernel, shots=1024):
    """Sample a CUDA-Q kernel and return its measurement counts.

    Selects the GPU (``nvidia``) target when available, falling back to CPU
    (``qpp-cpu``).

    :param kernel: A built CUDA-Q kernel.
    :param shots: Number of measurement shots.
    :return: A ``{bitstring: count}`` dict.
    """
    try:
        cudaq.set_target("nvidia")
        logger.info("Using GPU backend.")
    except Exception:
        cudaq.set_target("qpp-cpu")
        logger.info("GPU unavailable. Falling back to CPU backend.")

    freq = cudaq.sample(kernel, shots_count=shots)
    return {str(key): int(value) for key, value in freq.items()}


def _expectation_from_counts_cudaq(freq, observables, shots=1024):
    """Compute ``<O>`` for a Z/I Pauli string from CUDA-Q measurement counts.

    :param freq: A ``{bitstring: count}`` dict as returned by :func:`_counts_cudaq`.
    :param observables: Per-qubit Pauli string over ``{I, Z}``.
    :param shots: Number of shots the counts were drawn from.
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
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


def _expectation_cudaq_decoupled(
    kernel, observables, shots, prep_qubits, measured_order, num_qubits
):
    """Decoupled ``<O>`` for CUDA-Q gate-cut subcircuits with identity on a cut qubit.

    Mirrors the Qiskit/Qibo decoupling. With ``explicit_measurements=True`` each
    per-shot key is ``[cut-point outcomes in mz-call order] + [final measurement of
    q0..q(n-1)]``. For a measured cut qubit observed with identity, the
    measure-and-prepare sign is recovered from its (last) cut-point outcome ``m``
    as ``-(-1)**m`` rather than the traced-out final bit.

    :param kernel: Kernel built with ``explicit_final=True``.
    :param observables: Per-qubit Pauli string over ``{I, Z}`` (qubit-0-first).
    :param shots: Number of shots.
    :param prep_qubits: Renumbered qubits measured for the cut.
    :param measured_order: Cut-measured qubits in ``mz``-call order.
    :param num_qubits: Qubit count (width of the final-measurement segment).
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    result = cudaq.sample(kernel, shots_count=shots, explicit_measurements=True)
    n_mid = len(measured_order)
    last_pos = {qb: i for i, qb in enumerate(measured_order)}
    total = 0.0
    drawn = 0
    for key in result:
        count = result.count(key)
        drawn += count
        contribution = 1
        for q, obs in enumerate(observables):
            if q in prep_qubits and obs == "I" and q in last_pos:
                contribution *= -((-1) ** int(key[last_pos[q]]))
            elif obs == "Z":
                contribution *= (-1) ** int(key[n_mid + q])
            elif obs == "I":
                continue
            else:
                raise ValueError(f"Unsupported observable {obs}")
        total += contribution * count
    return total / (drawn or 1)


@task(returns=1)
def _expec_value_cudaq(subcirc: typing.Any, shots=1024, prep_qubits=None) -> float:
    """Execute the given circuit.

    :param subcirc: Circuit to execute.
    :param prep_qubits: Renumbered cut qubits measured for a gate cut; when any is
        observed with identity, its decomposition sign is recovered from the
        cut-point measurement (see :func:`_expectation_cudaq_decoupled`).
    :raises ValueError: If there is an unsupported observable.
    :return: The circuit expected value.
    """
    # Accept either an already-built kernel (``[kernel, observables]`` -- the
    # materialised/serial form) or a legacy ``[(dag, num_qubits), observables]``
    # payload whose kernel must still be built here.
    if (
        isinstance(subcirc, (list, tuple))
        and len(subcirc) == 2
        and isinstance(subcirc[1], str)
    ):
        payload, observables = subcirc[0], subcirc[1]
    else:
        payload, observables = subcirc, None

    is_legacy_dag = (
        isinstance(payload, tuple)
        and len(payload) == 2
        and hasattr(payload[0], "nodes")
    )

    # Identity on a measured cut qubit: build with explicit final measurements and
    # recover the sign from the cut-point outcome (decoupled from the observable).
    if is_legacy_dag and prep_qubits:
        kernel, kernel_observables, measured_order = dag_to_circuit_cudaq_kernel(
            payload[0], payload[1], explicit_final=True
        )
        if observables is None:
            observables = kernel_observables
        return _expectation_cudaq_decoupled(
            kernel, observables, shots, set(prep_qubits), measured_order, payload[1]
        )

    if is_legacy_dag:
        # Legacy payload: build the kernel here. Keep the observable string that
        # came with the payload (its DAG is already Observable-I-stripped, so the
        # kernel builder would otherwise recompute an all-Z observable).
        subcirc, kernel_observables = dag_to_circuit_cudaq_kernel(
            payload[0], payload[1]
        )
        if observables is None:
            observables = kernel_observables
    else:
        subcirc = payload

    freq = _counts_cudaq(subcirc, shots)
    return _expectation_from_counts_cudaq(freq, observables, shots)


@task(returns=1)
def _expec_value_qibo_qpu(subcirc, shots=1024, method="numpy"):
    """Evaluate a subcircuit on a BSC QPU via the shared-filesystem bridge.

    Ships the *unmeasured* subcircuit to the QPU listener
    (:func:`~Qdislib.utils.qpu_bridge.run_qpu_listener`, which adds the
    measurements and runs it on the hardware), then computes the expectation from
    the returned counts. The wait is bounded (see
    :func:`~Qdislib.utils.qpu_bridge.submit_and_wait`). ``method`` is unused here,
    kept for signature parity with the CPU/GPU evaluators.
    """
    if subcirc is None:
        return None

    if type(subcirc) != qibo.models.circuit.Circuit:
        observables = subcirc[1]
        subcirc = subcirc[0]
    else:
        observables = "Z" * subcirc.nqubits

    from Qdislib.utils.qpu_bridge import submit_and_wait

    result = submit_and_wait(subcirc)
    freq = result.counts()
    # Normalise by the actual shot count returned by the QPU (its num_bins may
    # differ from the classical-side ``shots`` argument).
    total_shots = sum(freq.values()) or 1
    return _expectation_from_counts_qibo(freq, observables, total_shots)


@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign(expectation_value, index):
    """Get the product of all expected values and apply change of sign if necessary.

    :param expectation_value: List of expected values.
    :param index: Number to be used to determine if the sign has to be changed.
    :return: The expected value.
    """
    expectation_value = [x for x in expectation_value if x is not None]
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


def wire_cutting_subcircuits(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    observables=None,
    software="qiskit",
):
    """
    Partition a quantum circuit into subcircuits via wire (edge) cutting.

    This is the entry point of the **recommended three-step workflow**: generate
    the subcircuits, evaluate each one (in your own loop, on whatever backend or
    hardware you like), then recombine the results. It gives you full control over
    *how* each subcircuit is executed, unlike :func:`wire_cutting`, which runs the
    whole pipeline in one call::

        import Qdislib as qd

        # 1. Cut the circuit into independently executable subcircuits.
        subcircuits = qd.wire_cutting_subcircuits(circuit, [("H_1", "CZ_2")], software="qibo")

        # 2. Evaluate each subcircuit (here with the built-in helper, but this
        #    loop is yours: dispatch to real QPUs, batch, cache, etc.).
        values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

        # 3. Recombine into the original circuit's expectation value.
        result = subcircuits.reconstruct(values)

    The cut metadata (number of cuts/components) is carried by the returned object,
    so :meth:`~Qdislib.utils.backends.CutSubcircuits.reconstruct` needs nothing else.

    Parameters
    ----------
    rand_qc : Union[qiskit.QuantumCircuit, qibo.models.Circuit, cudaq.PyKernelDecorator, nx.Graph]
        The quantum circuit to be cut. This can be a Qiskit ``QuantumCircuit``, a
        Qibo ``Circuit``, a CUDA-Q kernel, or an already-built DAG/networkx graph.

    cut : List[Tuple[str, str]]
        Wire cuts as tuples of the two gate identifiers (``"<GATE>_<INDEX>"``) that
        the cut wire runs between, e.g. ``("H_1", "CZ_2")``. Use
        :func:`~Qdislib.find_cut` to discover good candidates.

    observables : str, optional
        Per-qubit observable string (e.g. ``"ZZ"``). When given, the measurement
        basis is rotated accordingly before cutting. Defaults to the all-``Z``
        observable.

    software : str, optional
        The backend the subcircuits are emitted for, which must match the backend
        used to evaluate them:
        - ``"qiskit"`` (default): Qiskit ``QuantumCircuit`` objects
        - ``"qibo"``: Qibo ``Circuit`` objects
        - ``"cudaq"``: CUDA-Q ``(dag, num_qubits)`` payloads

    Returns
    -------
    CutSubcircuits
        A ``list`` subclass holding every subcircuit needed to reconstruct the
        original expectation value (iterate/index it like a list), plus the
        ``number_cuts``/``number_components`` metadata and a
        :meth:`~Qdislib.utils.backends.CutSubcircuits.reconstruct` method.

    Raises
    ------
    ValueError
        If ``cut`` is malformed or an unsupported ``software`` option is provided.

    Notes
    -----
    - If the circuit has multiple disconnected components, each is processed separately.
    - When ``cut`` is empty the whole circuit is returned as a single subcircuit with
      updated qubit indices.
    - Edges in ``cut`` are assumed to exist in the graph; otherwise they are ignored.

    Examples
    --------
    >>> import Qdislib as qd
    >>> from qibo import models, gates
    >>> circuit = models.Circuit(2)
    >>> circuit.add(gates.H(0)); circuit.add(gates.CZ(0, 1)); circuit.add(gates.RY(0, 0.8))
    >>> subcircuits = qd.wire_cutting_subcircuits(circuit, [("H_1", "CZ_2")], software="qibo")
    >>> values = [qd.expectation_value(s, shots=8192) for s in subcircuits]
    >>> subcircuits.reconstruct(values)  # doctest: +SKIP
    -0.717...
    """
    pattern = re.compile(r"^[A-Za-z]+_\d+$")  # e.g. CX_1, H_0, RZ_3

    if not isinstance(cut, list):
        raise InvalidCutError(
            "Invalid format of wire cuts: expected a list of tuples with two strings each."
        )
    for t in cut:
        if not (
            isinstance(t, tuple) and len(t) == 2 and all(isinstance(x, str) for x in t)
        ):
            raise InvalidCutError(
                "Invalid format of wire cuts: expected a list of tuples with two strings each."
            )
        if not all(pattern.match(x) for x in t):
            raise InvalidCutError(
                "Invalid wire name format: each string in wire_cut tuples must match '<GATE>_<NUMBER>', e.g. 'CX_1'."
            )

    dag = convert_to_dag(rand_qc, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    # Serial -> native backend circuit; PyCOMPSs -> serializable Subcircuit.
    dispatch_map = {
        "qibo": lambda d, n: make_subcircuit(d, n, "qibo"),
        "qiskit": lambda d, n: make_subcircuit(d, n, "qiskit"),
        "cudaq": lambda d, n: make_subcircuit(d, n, "cudaq"),
    }

    if networkx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected())
        ]
        subcircuits = []
        number_cuts = 0
        for s in S:
            num_qubits = _max_qubit(s)
            tmp_cuts = [c for c in cut if s.has_edge(*c)]
            if tmp_cuts:
                number_cuts += len(_expand_wire_cuts(s, tmp_cuts))
                subcircuits += _generate_wire_cutting_subcircuits(
                    s, tmp_cuts, num_qubits=num_qubits, software=software
                )
            else:
                s_new, highest_qubit = _update_qubits(s)
                try:
                    subcircuits.append(dispatch_map[software](s_new, highest_qubit))
                except KeyError:
                    raise UnsupportedBackendError(
                        f"Unsupported software backend: '{software}'. "
                        f"Expected one of {list(dispatch_map.keys())}."
                    )
    else:
        num_qubits = _max_qubit(dag)
        if cut:
            number_cuts = len(_expand_wire_cuts(dag, cut))
            subcircuits = _generate_wire_cutting_subcircuits(
                dag, cut, num_qubits=num_qubits, software=software
            )
        else:
            number_cuts = 0
            s_new, highest_qubit = _update_qubits(dag)
            try:
                subcircuits = [dispatch_map[software](s_new, highest_qubit)]
            except KeyError:
                raise UnsupportedBackendError(
                    f"Unsupported software backend: '{software}'. "
                    f"Expected one of {list(dispatch_map.keys())}."
                )

    number_components = max(1, len(subcircuits) // (8**number_cuts))
    return CutSubcircuits(
        subcircuits,
        number_cuts=number_cuts,
        number_components=number_components,
        reconstructor=wire_cutting_subcircuit_reconstruction,
    )


def _generate_wire_cutting_subcircuits(
    dag: typing.Any,
    edges_to_replace: typing.List[typing.Tuple[str, str]],
    num_qubits: int,
    software,
) -> typing.List[typing.Any]:
    """Apply wire cuts to the DAG and build the resulting subcircuits.

    Mirrors :func:`_generate_wire_cutting` but returns subcircuits (in the chosen
    ``software`` backend) instead of evaluating them. Each ``(source, target)``
    cut is expanded to one cut per shared qubit (see :func:`_expand_wire_cuts`),
    so two gates sharing several qubits are separated correctly.

    :param dag: The circuit DAG to modify (mutated in place).
    :param edges_to_replace: List of ``(source, target)`` gate-name cut pairs.
    :param num_qubits: Highest qubit index; fresh qubit lines are numbered above it.
    :param software: Target backend ("qibo", "qiskit" or "cudaq").
    :return: Flat list of subcircuits across all ``8 ** n_cuts`` basis terms.
    """
    wire_cuts = _expand_wire_cuts(dag, edges_to_replace)

    for index, (source, target, qubit) in enumerate(wire_cuts, start=1):
        if dag.has_edge(source, target):
            dag.remove_edge(source, target)

        successors = [node for node in dag.nodes if dag.has_edge(target, node)]
        for node in [target] + successors:
            node_qubits = dag.nodes[node].get("qubits")
            if node_qubits and qubit in node_qubits:
                dag.nodes[node]["qubits"] = tuple(
                    num_qubits + index if q == qubit else q for q in node_qubits
                )

        dag.add_node(f"O_{index}", gate="S", qubits=(qubit,), parameters=())
        dag.add_node(
            f"PS_{index}", gate="T", qubits=(num_qubits + index,), parameters=()
        )
        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")

    # Serial -> native backend circuit; PyCOMPSs -> serializable Subcircuit.
    dispatch_map = {
        "qibo": lambda d, n: make_subcircuit(d, n, "qibo"),
        "qiskit": lambda d, n: make_subcircuit(d, n, "qiskit"),
        "cudaq": lambda d, n: make_subcircuit(d, n, "cudaq"),
    }

    subcircuits = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag.copy())

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )
        graph_components = [networkx.DiGraph().copy() for _ in range(num_components)]

        _generate_subcircuits_wire_cutting(
            copy_graph,
            num_qubits + len(wire_cuts),
            index,
            wire_cuts,
            graph_components,
        )

        for s in graph_components:
            s_new, highest_qubit = _update_qubits(s)
            try:
                subcirc = dispatch_map[software](s_new, highest_qubit)
            except KeyError:
                raise UnsupportedBackendError(
                    f"Unsupported software backend: '{software}'. "
                    f"Expected one of {list(dispatch_map.keys())}."
                )
            subcircuits.append(subcirc)
    return subcircuits


def wire_cutting_subcircuit_reconstruction(
    expected_values, number_cuts, number_components=2
):
    """
    Reconstruct the expectation value of a quantum circuit from the individual
    expectation values of its subcircuits, based on wire cutting techniques.

    This function assumes that the circuit was cut along a specified number of wires,
    and each resulting subcircuit was simulated separately. The reconstruction process
    involves:

    - Grouping subcircuit results into pairs.
    - Applying sign changes based on index parity.
    - Averaging the results with a normalization factor dependent on the number of cuts.

    In the recommended workflow you don't call this directly: the object returned
    by :func:`wire_cutting_subcircuits` already carries ``number_cuts`` and
    ``number_components``, so ``subcircuits.reconstruct(values)`` forwards here for
    you. Call this function explicitly only when reconstructing from values whose
    cut metadata you are tracking yourself.

    :param expected_values: A flat list of expectation values, where each group of
        ``number_components`` corresponds to one wire-cut configuration.
    :param number_cuts: The number of wire cuts made in the original circuit; this
        sets the ``1 / 2^(number_cuts)`` normalization factor.
    :param number_components: The number of independent parts the subcircuits split into.
    :return: The reconstructed expectation value of the original (uncut) circuit.

    Example:
        >>> values = [0.8, -0.6, 0.7, -0.5] * 4
        >>> wire_cutting_subcircuit_reconstruction(values, number_cuts=1, number_components=2)
    """
    # Resolve any deferred (cache placeholder-hit) results once all subcircuit
    # tasks have finished: synchronise to the master, then refresh each
    # PendingResult from the now-populated cache. Doing this after the barrier
    # (rather than inside the per-group sign tasks) avoids worker starvation when
    # a deferred key's owner is still running in another group.
    from Qdislib.utils.circuit_cache import _resolve_value

    expected_values = compss_wait_on(expected_values)
    expected_values = [_resolve_value(v) for v in expected_values]

    if number_cuts == 0:
        return math.prod([v for v in expected_values if v is not None])

    expected_values = [
        expected_values[i : i + number_components]
        for i in range(0, len(expected_values), number_components)
    ]
    results = []
    for index, exp_values in enumerate(expected_values):
        result = _change_sign(exp_values, index)
        results.append(result)

    results = compss_wait_on(results)
    final_recons = 1 / (2**number_cuts) * sum(results)
    return final_recons
