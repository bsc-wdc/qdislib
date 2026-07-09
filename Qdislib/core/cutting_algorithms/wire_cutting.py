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
    "wire_cutting_sampling",
    "wire_cutting_sampling_subcircuits",
    "WireSampledCutSubcircuits",
    "wire_cutting_subcircuits",
    "wire_cutting_subcircuit_reconstruction",
]

import itertools
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

try:
    import pennylane as qml

    pennylane_available = True
except ImportError:  # pennylane is optional; functions that need it fail when called
    qml = None
    pennylane_available = False

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN

from Qdislib.utils.graph_qibo import (
    _max_qubit,
    _update_qubits,
    _remove_red_edges,
)
from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel
from Qdislib.utils.graph_pennylane import (
    _apply_pennylane_gate,
    _build_pennylane_qnode,
)


from Qdislib.utils.backends import (
    convert_to_dag,
    change_basis_dag,
    normalize_observables,
    make_subcircuit,
    _make_subcircuit_task,
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
    optimized: bool = False,
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
    optimized : bool, default=False
        Use the fragment-separated reconstruction, which runs ``3**n + 4**n``
        subcircuits for ``n`` cuts instead of ``8**n`` by evaluating the measure and
        prepare fragments independently and combining their cached values
        classically (identical result, far fewer circuits). Falls back to the
        standard method for cut topologies it cannot separate (e.g. two cuts in
        series on one wire).

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

    if optimized and cut:
        # Fragment-separated reconstruction: 3**n + 4**n subcircuits instead of
        # 8**n. Falls back to the standard method for cut topologies it cannot
        # separate (a fragment carrying both a measure and a prepare marker) or if
        # anything goes wrong, so optimized=True never changes the result.
        try:
            optimized_value = _wire_cutting_optimized(
                dag.copy(),
                _expand_wire_cuts(dag, cut),
                _max_qubit(dag),
                shots,
                backend,
                software,
            )
            # The optimized path returns a single reconstruction future (task ->
            # task, no master barrier); wait on it once here.
            if optimized_value is not None:
                return compss_wait_on(optimized_value)
            reason = "not applicable to this cut topology"
        except Exception as exc:  # never let the optimization break the computation
            reason = f"failed ({exc})"
        logger.info("Optimized wire cutting %s; using the standard method.", reason)

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
                    s,
                    cut_wires,
                    num_qubits=num_qubits,
                    shots=shots,
                    backend=backend,
                    qpu=qpu,
                    software=software,
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
                dag,
                cut_wires,
                num_qubits=num_qubits,
                shots=shots,
                backend=backend,
                qpu=qpu,
                software=software,
            )

            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(cut_wires)) * sum(results)
        else:
            final_recons = _build_and_run_subcircuit(dag, software, shots, backend, qpu)
        return final_recons


def _wire_cut_sampling_components(dag, wire_cuts, num_qubits):
    """Build the wire-cut fragments for every 8**n_cuts configuration (one pass).

    Mirrors :func:`_generate_wire_cutting`'s cut surgery (measure node ``O_i`` on the
    original qubit, prepare node ``PS_i`` on a fresh line ``num_qubits+i``, downstream
    relabelled onto that fresh line) and its base-8 configuration loop, but returns
    the raw fragment DAGs to be *sampled* rather than evaluating expectations.

    Returns ``(components, per_config_counts, cut_source_index, prepared_to_orig)``:
    a flat, configuration-major list of fragment DAGs; the number of fragments per
    configuration; the map ``cut-source qubit -> 1-based cut index`` (its base-8 digit
    decides whether it is a ``Z`` decomposition factor -- digits 2-7 -- or observed
    with ``I`` -- digits 0/1); and the map ``fresh line -> original qubit`` (a cut
    wire's output lives on its fresh line).
    """
    dag = dag.copy()
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

    # Cut-source qubit -> its 1-based cut index. The base-8 digit for that cut sets
    # the measure basis: identity for digits 0/1 (observe the cut-source with I), else
    # X/Y/Z (a Z decomposition factor). See :func:`_reconstruct_wire_cut_distribution`.
    cut_source_index = {
        qubit: index for index, (_, _, qubit) in enumerate(wire_cuts, start=1)
    }
    prepared_to_orig = {
        num_qubits + index: qubit
        for index, (_, _, qubit) in enumerate(wire_cuts, start=1)
    }

    components = []
    per_config_counts = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag.copy())
        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )
        graph_components = [networkx.DiGraph() for _ in range(num_components)]
        _generate_subcircuits_wire_cutting(
            copy_graph, num_qubits + len(wire_cuts), index, wire_cuts, graph_components
        )
        per_config_counts.append(len(graph_components))
        components.extend(graph_components)
    return components, per_config_counts, cut_source_index, prepared_to_orig


def _reconstruct_wire_cut_distribution(
    configs, num_cuts, num_qubits, cut_source_index, prepared_to_orig
):
    """Reconstruct the output distribution from sampled wire-cut fragments.

    Wire counterpart of ``_reconstruct_gate_cut_distribution``: each fragment supplies
    Pauli-``Z`` correlators ``<Z^b>`` (from one sampling); per configuration these
    multiply (signed by the base-8 rule -- flip on octal digits {3,5,7} -- over
    ``2**num_cuts``) into the global ``<Z^A>``, and ``P(x)=2**-N * sum_A (-1)**(A.x)
    <Z^A>`` (inverse Walsh-Hadamard). A measured cut-source qubit is a decomposition
    factor: for the X/Y/Z terms (base-8 digit >= 2) it is in the Z-string of its
    fragment, for the identity terms (digits 0/1) it is observed with ``I``. A cut
    wire's output lives on its fresh line, so ``prepared_to_orig`` maps that line back
    to the original output qubit. The signed result is abs-normalised. This matches
    the Perlin et al. (arXiv:2005.12702) reconstruction in one
    fragment pass instead of ``2**N`` expectation runs.
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard

    dim = 2**num_qubits
    z_correlators = np.zeros(dim)
    for index, per_component in enumerate(configs):
        # Base-8 digits, MSB-first: octal[i-1] is cut i's measure/prepare term.
        octal = oct(index)[2:].zfill(num_cuts)
        sign, digits = 1, index
        while digits:
            if digits % 8 in (3, 5, 7):
                sign = -sign
            digits //= 8

        contribution = np.full(dim, sign / (2**num_cuts), dtype=float)
        for corr, original in per_component:
            component_values = np.empty(dim)
            for global_mask in range(dim):
                local_mask = 0
                for local_index, qubit in enumerate(original):
                    cut = cut_source_index.get(qubit)
                    if cut is not None:
                        # Measured cut-source qubit: a Z decomposition factor for the
                        # X/Y/Z terms (digit >= 2); the identity terms (0, 1) observe
                        # it with I, so it stays out of the Z-string.
                        if int(octal[cut - 1]) >= 2:
                            local_mask |= 1 << local_index
                    elif (global_mask >> prepared_to_orig.get(qubit, qubit)) & 1:
                        # Output qubit (its own line, or a cut wire's fresh line).
                        local_mask |= 1 << local_index
                component_values[global_mask] = corr.get(local_mask, 0.0)
            contribution *= component_values
        z_correlators += contribution

    probabilities = _fast_walsh_hadamard(z_correlators) / dim
    keys = ["".join(str((x >> q) & 1) for q in range(num_qubits)) for x in range(dim)]
    norm = float(np.abs(probabilities).sum()) or 1.0
    return {key: abs(float(p)) / norm for key, p in zip(keys, probabilities)}


@task(returns=1, sampled=COLLECTION_IN)
def _reconstruct_wire_cut_distribution_task(
    sampled, num_cuts, num_qubits, per_config_counts, cut_source_index, prepared_to_orig
):
    """Task wrapper: regroup the flat per-fragment samples by configuration and
    reconstruct the distribution. ``COLLECTION_IN`` feeds the fragment futures
    straight in (task -> task, no master barrier before reconstruction)."""
    configs = []
    start = 0
    for count in per_config_counts:
        configs.append(sampled[start : start + count])
        start += count
    return _reconstruct_wire_cut_distribution(
        configs, num_cuts, num_qubits, cut_source_index, prepared_to_orig
    )


@task(
    returns=1,
    measure_flat=COLLECTION_IN,
    prepare_flat=COLLECTION_IN,
    plain_flat=COLLECTION_IN,
)
def _reconstruct_optimized_wire_distribution_task(
    measure_flat,
    measure_keys,
    measure_cuts,
    prepare_flat,
    prepare_keys,
    prepare_cuts,
    plain_flat,
    num_cuts,
    num_qubits,
    cut_source_index,
    prepared_to_orig,
):
    """Assemble the ``8**n`` per-config fragment correlators from the deduped samples,
    then reconstruct the distribution with the standard combiner.

    The fragment futures (each ``(corr, original)``) feed in through
    ``COLLECTION_IN`` (task -> task, no master barrier). The deduped measure corrs
    are looked up by measurement basis -- a Z-basis correlator serves both the
    identity (digits 0/1) and Z (digits 6/7) terms through the reconstruction's local
    mask -- and the prepare corrs by prepared state, with |-> / |-i> expanded as
    linear combinations of the base-state correlators (correlators are linear in the
    input state). Reuses :func:`_reconstruct_wire_cut_distribution` unchanged.
    """
    basis_of = {0: "Z", 1: "Z", 2: "X", 3: "X", 4: "Y", 5: "Y", 6: "Z", 7: "Z"}
    measure_corr = {}
    for (frag_idx, basis_combo), value in zip(measure_keys, measure_flat):
        measure_corr.setdefault(frag_idx, {})[basis_combo] = value
    prepare_corr = {}
    for (frag_idx, prep_combo), value in zip(prepare_keys, prepare_flat):
        prepare_corr.setdefault(frag_idx, {})[prep_combo] = value

    configs = []
    for index in range(8**num_cuts):
        settings = [(index // (8 ** (num_cuts - 1 - k))) % 8 for k in range(num_cuts)]
        per_component = []
        for frag_idx, cuts in enumerate(measure_cuts):
            basis_combo = tuple(basis_of[settings[cut - 1]] for cut in cuts)
            per_component.append(measure_corr[frag_idx][basis_combo])
        for frag_idx, cuts in enumerate(prepare_cuts):
            expansions = [_WIRE_PREP[settings[cut - 1]] for cut in cuts]
            combined, original = {}, None
            for choice in itertools.product(*expansions):
                coefficient = 1.0
                base_combo = []
                for weight, base_state in choice:
                    coefficient *= weight
                    base_combo.append(base_state)
                corr, original = prepare_corr[frag_idx][tuple(base_combo)]
                for mask, value in corr.items():
                    combined[mask] = combined.get(mask, 0.0) + coefficient * value
            per_component.append((combined, original))
        per_component.extend(plain_flat)
        configs.append(per_component)

    return _reconstruct_wire_cut_distribution(
        configs, num_cuts, num_qubits, cut_source_index, prepared_to_orig
    )


def _wire_sampling_optimized(
    dag, wire_cuts, num_qubits, shots, software, qpu, qpu_dict, seed
):
    """Fragment-separated sampling reconstruction (``3**n + 4**n`` fragment samplings).

    The sampling counterpart of :func:`_wire_cutting_optimized`: sample each measure
    fragment once per measurement basis and each prepare fragment once per base
    preparation, then reconstruct the distribution over the ``8**n`` configurations
    from those cached correlators. Returns the reconstruction future, or ``None`` if
    the cut topology is not fragment-separable (mixed fragment).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _sample_gate_subcircuit

    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, plain_frags = fragments

    cut_source_index = {
        qubit: index for index, (_, _, qubit) in enumerate(wire_cuts, start=1)
    }
    prepared_to_orig = {
        num_qubits + index: qubit
        for index, (_, _, qubit) in enumerate(wire_cuts, start=1)
    }

    counter = itertools.count()

    def _sample(component):
        derived = None if seed is None else seed + next(counter)
        return _sample_gate_subcircuit(
            component, software, shots, qpu, qpu_dict, seed=derived
        )

    measure_flat, measure_keys = [], []
    for frag_idx, (sub, cuts, _q) in enumerate(measure_frags):
        for basis_combo in itertools.product(_MEASURE_BASES, repeat=len(cuts)):
            graph = sub.copy()
            for cut, basis in zip(cuts, basis_combo):
                _apply_wire_measure_option(graph, cut, basis + "P")
            measure_flat.append(_sample(graph))
            measure_keys.append((frag_idx, basis_combo))
    prepare_flat, prepare_keys = [], []
    for frag_idx, (sub, cuts) in enumerate(prepare_frags):
        for prep_combo in itertools.product(_PREP_STATES, repeat=len(cuts)):
            graph = sub.copy()
            for cut, state in zip(cuts, prep_combo):
                _apply_wire_prep_state(graph, cut, state)
            prepare_flat.append(_sample(graph))
            prepare_keys.append((frag_idx, prep_combo))
    plain_flat = [_sample(sub.copy()) for sub in plain_frags]

    return _reconstruct_optimized_wire_distribution_task(
        measure_flat,
        measure_keys,
        [cuts for _s, cuts, _q in measure_frags],
        prepare_flat,
        prepare_keys,
        [cuts for _s, cuts in prepare_frags],
        plain_flat,
        len(wire_cuts),
        num_qubits,
        cut_source_index,
        prepared_to_orig,
    )


def wire_cutting_sampling(
    circuit,
    cut,
    shots=1024,
    software="qiskit",
    sync=True,
    qpu=False,
    qpu_dict=None,
    seed=None,
    optimized: bool = False,
):
    """Reconstruct the output **probability distribution** of a wire-cut circuit.

    The wire-cutting counterpart of
    :func:`~Qdislib.core.cutting_algorithms.gate_cutting.gate_cutting_sampling`:
    instead of a single expectation value, reconstruct the full distribution over
    output bitstrings. Each wire-cut fragment is sampled **once** (so under PyCOMPSs
    the fragments run in parallel), its Pauli-``Z`` correlators are combined via the
    wire-cut quasiprobability decomposition (base-8; Perlin et al., arXiv:2005.12702)
    into ``<Z^A>``, and ``P(x)=2**-N * sum_A (-1)**(A.x) <Z^A>`` is the inverse
    Walsh-Hadamard transform, abs-normalised into a valid distribution.

    :param circuit: A Qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param cut: The wire cut(s), as returned by :func:`~Qdislib.find_cut`.
    :param shots: Measurement shots per fragment.
    :param software: Execution backend, or the name of any registered backend.
    :param sync: Under PyCOMPSs, wait for the tasks and return the distribution; set
        ``False`` to keep the future unresolved (advanced).
    :param qpu: Route each fragment's sampling to a QPU (via
        ``build_execution_backend``) instead of a local simulator.
    :param qpu_dict: QPU-name -> max-qubit map selecting the QPU (e.g. ``{"MN_Ona": 4}``).
    :param seed: Optional seed for a reproducible distribution. Each fragment gets a
        distinct derived seed (``seed + i``); ``None`` keeps the previous
        non-deterministic sampling. Ignored on a QPU.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _sample_gate_subcircuit,
        _dag_num_qubits,
    )

    dag = convert_to_dag(circuit)
    num_qubits = _dag_num_qubits(dag)
    wire_cuts = _expand_wire_cuts(dag, cut)

    if not wire_cuts:
        from Qdislib.core.backend.registry import build_execution_backend

        backend = build_execution_backend(
            software, num_qubits=num_qubits, qpu=bool(qpu), qpu_dict=qpu_dict
        )
        counts = backend.evaluate_samples(dag, num_qubits, shots, seed=seed)
        total = sum(counts.values()) or 1
        return {k: v / total for k, v in counts.items()}

    if optimized:
        # Fragment-separated sampling: 3**n + 4**n fragment samplings instead of
        # 8**n. Falls back to the standard method for non-separable topologies or
        # on error, so optimized=True never changes the (shot-noise) result.
        try:
            optimized_distribution = _wire_sampling_optimized(
                dag.copy(), wire_cuts, num_qubits, shots, software, qpu, qpu_dict, seed
            )
            if optimized_distribution is not None:
                return (
                    compss_wait_on(optimized_distribution)
                    if sync
                    else optimized_distribution
                )
            reason = "not applicable to this cut topology"
        except Exception as exc:  # never let the optimization break the computation
            reason = f"failed ({exc})"
        logger.info(
            "Optimized wire-cut sampling %s; using the standard method.", reason
        )

    components, per_config_counts, cut_source_index, prepared_to_orig = (
        _wire_cut_sampling_components(dag, wire_cuts, num_qubits)
    )
    sampled = [
        _sample_gate_subcircuit(
            comp,
            software,
            shots,
            qpu,
            qpu_dict,
            seed=None if seed is None else seed + i,
        )
        for i, comp in enumerate(components)
    ]
    distribution = _reconstruct_wire_cut_distribution_task(
        sampled,
        len(wire_cuts),
        num_qubits,
        per_config_counts,
        cut_source_index,
        prepared_to_orig,
    )
    return compss_wait_on(distribution) if sync else distribution


class WireSampledCutSubcircuits(list):
    """Wire-cut sampling fragments that also know how to be recombined.

    The wire, sampling counterpart of
    :class:`~Qdislib.core.cutting_algorithms.gate_cutting.SampledCutSubcircuits`:
    behaves like the plain ``list`` of fragments (iterate/index/``len``) and carries
    the metadata to reconstruct the output distribution::

        subs = qd.wire_cutting_sampling_subcircuits(circuit, cut, software="qiskit")
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(sampled)
    """

    def __init__(
        self,
        subcircuits,
        *,
        number_cuts,
        num_qubits,
        per_config_counts,
        cut_source_index,
        prepared_to_orig,
    ):
        super().__init__(subcircuits)
        #: Number of wire cuts applied.
        self.number_cuts = number_cuts
        #: Qubit count of the original circuit (distribution width).
        self.num_qubits = num_qubits
        self._per_config_counts = per_config_counts
        self._cut_source_index = cut_source_index
        self._prepared_to_orig = prepared_to_orig

    def reconstruct(self, sampled_results):
        """Recombine per-fragment samples into the output distribution.

        :param sampled_results: One :func:`sample_subcircuit` result per fragment,
            in the same order as the fragments in this object.
        :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
        """
        distribution = _reconstruct_wire_cut_distribution_task(
            sampled_results,
            self.number_cuts,
            self.num_qubits,
            self._per_config_counts,
            self._cut_source_index,
            self._prepared_to_orig,
        )
        return compss_wait_on(distribution)


def wire_cutting_sampling_subcircuits(
    circuit, cut, software="qiskit", qpu=False, qpu_dict=None
):
    """Partition a circuit into **wire-cut sampling** fragments (three-step workflow).

    The wire counterpart of
    :func:`~Qdislib.core.cutting_algorithms.gate_cutting.gate_cutting_sampling_subcircuits`:
    generate the fragments, sample each yourself, then reconstruct the output
    distribution::

        subs = qd.wire_cutting_sampling_subcircuits(circuit, cut, software="qiskit")
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(sampled)

    :param circuit: A Qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param cut: The wire cut(s), as returned by :func:`~Qdislib.find_cut`.
    :param software: Backend the fragments are sampled on.
    :param qpu: Route fragment sampling to a QPU when each fragment is later sampled.
    :param qpu_dict: QPU-name -> max-qubit map selecting the QPU (e.g. ``{"MN_Ona": 4}``).
    :return: A :class:`WireSampledCutSubcircuits` (a ``list`` of fragments plus the
        cut metadata and a ``reconstruct`` method).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _dag_num_qubits,
        _SampleUnit,
    )

    dag = convert_to_dag(circuit)
    num_qubits = _dag_num_qubits(dag)
    wire_cuts = _expand_wire_cuts(dag, cut)
    components, per_config_counts, cut_source_index, prepared_to_orig = (
        _wire_cut_sampling_components(dag, wire_cuts, num_qubits)
    )
    units = [
        _SampleUnit(component, software, qpu, qpu_dict) for component in components
    ]
    return WireSampledCutSubcircuits(
        units,
        number_cuts=len(wire_cuts),
        num_qubits=num_qubits,
        per_config_counts=per_config_counts,
        cut_source_index=cut_source_index,
        prepared_to_orig=prepared_to_orig,
    )


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
def _create_subcircuit(component):
    """Prepare one wire-cut subcircuit for execution (COMPSs task).

    The *creation* half, split from execution so the two are independently
    schedulable/routable. Renumbers the component to a standalone ``0..n``
    subcircuit DAG. Returns a serializable spec ``(dag, num_qubits)``; the native
    circuit is built later, in :func:`_execute_subcircuit` (a Qibo collapse circuit
    / CUDA-Q kernel cannot cross a task boundary).
    """
    s_new, highest_qubit = _update_qubits(component)
    return (s_new, highest_qubit)


@task(returns=1)
def _execute_subcircuit(prepared, software, shots, backend, qpu):
    """Build + run a prepared wire-cut subcircuit; return its expectation value.

    The *execution* half (COMPSs task) -- this is where the native circuit is built
    and run, so it is the task to route (e.g. ``@constraint`` or a QPU backend) to
    control where subcircuits execute. Backend-agnostic dispatch; a QPU target
    ("MN_Ona") overrides the software choice. Wire cutting has no measure-and-prepare
    sign, so there are no prep_qubits; ``backend`` names the Qibo platform (a generic
    ``extra`` config key that non-Qibo backends ignore). The native circuit never
    leaves this task.
    """
    s_new, highest_qubit = prepared

    from Qdislib.core.backend.registry import build_execution_backend, build_backend
    from Qdislib.core.backend.base_backend import BackendConfig

    if qpu == "MN_Ona":
        qpu_backend = build_backend(
            "BSC_QPU_BACKEND", BackendConfig(num_qubits=highest_qubit)
        )
        return qpu_backend.evaluate_expectation(
            s_new, highest_qubit, (), shots, backend or "numpy"
        )

    extra = {"platform": backend} if backend else {}
    exec_backend = build_execution_backend(software, num_qubits=highest_qubit, **extra)
    return exec_backend.evaluate_expectation(s_new, highest_qubit, (), shots)


def _build_and_run_subcircuit(component, software, shots, backend, qpu):
    """Prepare then execute one wire-cut subcircuit as two COMPSs tasks.

    Orchestrates :func:`_create_subcircuit` (prepare) then :func:`_execute_subcircuit`
    (build + run), so subcircuit creation and execution are separate, independently
    routable tasks (the prepared-spec future feeds create -> execute directly,
    task -> task). Called from the master, so it submits both without nesting.
    Returns the execution future (the subcircuit's expectation value).
    """
    prepared = _create_subcircuit(component)
    return _execute_subcircuit(prepared, software, shots, backend, qpu)


def _insert_wire_cut_markers(dag, wire_cuts, num_qubits):
    """Insert the measure (``O_i``) and prepare (``PS_i``) marker nodes for each cut.

    For every single-wire cut ``(source, target, qubit)``, removes the ``source ->
    target`` edge, relabels ``qubit`` on ``target`` and its successors onto a fresh
    line ``num_qubits + i`` (so the prepared side is disconnected from the measured
    side), and adds ``O_i`` (on the original ``qubit``, downstream of ``source``) and
    ``PS_i`` (on the fresh line, upstream of ``target``). The marker gates are
    placeholders (``S`` / ``T``); the standard and optimized paths overwrite them per
    measure/prepare basis. Mutates ``dag`` in place and returns it.

    :param dag: The circuit DAG (mutated in place).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :return: The mutated DAG.
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
    return dag


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
    _insert_wire_cut_markers(dag, wire_cuts, num_qubits)

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
def _apply_wire_substitution(
    updated_dag: typing.Any,
    num_qubits: int,
    idx: int,
    edges_to_replace: typing.List[typing.Tuple[str, str]],
) -> typing.Any:
    """Apply the base-8 measure/prepare substitution for wire-cut configuration ``idx``
    to the O_i/PS_i marker nodes (mutates ``updated_dag`` in place, no component split).
    Factored out of :func:`_generate_subcircuits_wire_cutting` so the mixed wire+gate
    driver can compose the wire substitution with the gate substitution on one DAG."""
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
    _apply_wire_substitution(updated_dag, num_qubits, idx, edges_to_replace)
    updated_dag = _remove_red_edges(updated_dag)
    for i, c in enumerate(networkx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")
    return updated_dag


# ---------------------------------------------------------------------------
# Optimized wire cutting (fragment-separated: 3**n measure + 4**n prepare
# subcircuits instead of 8**n).
#
# The standard base-8 decomposition runs every (measure, prepare) configuration
# as a fresh circuit -> 8**n executions, re-running each fragment redundantly.
# The measure ("top") and prepare ("bottom") fragments are independent, so we run
# them ONCE each and combine their cached values in the 8**n sum classically. Two
# algebraic reductions keep the per-cut circuit count minimal: a Z-basis
# measurement yields BOTH <Z> and <I> from one circuit (4 measure observables ->
# 3 bases), and |-> / |-i> are linear combinations of prepared states (6 prepare
# states -> 4 circuits). Result: sum over fragments of 3**(measures) + 4**(preps).
# ---------------------------------------------------------------------------

_MEASURE_BASES = ("X", "Y", "Z")
_PREP_STATES = ("0", "1", "+", "+i")

#: Base-8 cut setting -> (measurement basis, extraction) for the measure fragment.
#: extraction "P" reads the measured Pauli (its eigenvalue), "I" the identity
#: (drop the cut qubit) -- both come from the same Z-basis circuit for settings 0/1/6/7.
_WIRE_MEASURE = {
    0: ("Z", "I"),
    1: ("Z", "I"),  # measure I (from the Z circuit), prepare |0>/|1>
    2: ("X", "P"),
    3: ("X", "P"),  # measure X, prepare |+>/|->
    4: ("Y", "P"),
    5: ("Y", "P"),  # measure Y, prepare |+i>/|-i>
    6: ("Z", "P"),
    7: ("Z", "P"),  # measure Z, prepare |0>/|1>
}

#: Base-8 cut setting -> prepared state as a linear combination of the 4 base
#: preparations (|-> = |0> + |1> - |+>, |-i> = |0> + |1> - |+i>).
_WIRE_PREP = {
    0: [(1.0, "0")],
    1: [(1.0, "1")],
    2: [(1.0, "+")],
    3: [(1.0, "0"), (1.0, "1"), (-1.0, "+")],
    4: [(1.0, "+i")],
    5: [(1.0, "0"), (1.0, "1"), (-1.0, "+i")],
    6: [(1.0, "0")],
    7: [(1.0, "1")],
}

#: Cut settings whose quasiprobability sign is -1 (matches _reconstruction base-8).
_WIRE_SIGN_FLIP = {3, 5, 7}


def _apply_wire_prep_state(graph, cut_index, state):
    """Set a prepare marker ``PS_i`` to prepare one of the 4 base states, in place.

    |0> drops the marker (the default input), |1> -> X, |+> -> H, |+i> -> H then S.
    """
    node = f"PS_{cut_index}"
    if node not in graph.nodes:
        return
    if state == "0":
        graph.remove_node(node)
    elif state == "1":
        graph.nodes[node]["gate"] = "X"
    elif state == "+":
        graph.nodes[node]["gate"] = "H"
    elif state == "+i":
        graph.nodes[node]["gate"] = "S"
        rotate = f"{node}_Hi"
        graph.add_node(
            rotate, gate="H", qubits=graph.nodes[node]["qubits"], parameters=()
        )
        graph.add_edge(rotate, node, color="blue")  # H then S = |+i>


def _apply_wire_measure_option(graph, cut_index, option):
    """Set a measure marker ``O_i`` to observe one materialized option, in place.

    "ZI" -> Observable I (identity on the cut qubit); "XP" -> H (measure X);
    "YP" -> SDG then H (measure Y); "ZP" -> drop the marker (measure Z directly).
    Mirrors the standard base-8 measure side but per (basis, observable), so a single
    ``expectation_value`` call yields the term the reconstruction needs.
    """
    node = f"O_{cut_index}"
    if node not in graph.nodes:
        return
    if option == "ZI":
        graph.nodes[node]["gate"] = "Observable I"
    elif option == "XP":
        graph.nodes[node]["gate"] = "H"
    elif option == "YP":
        graph.nodes[node]["gate"] = "SDG"
        rotate = f"{node}_HY"
        graph.add_node(
            rotate, gate="H", qubits=graph.nodes[node]["qubits"], parameters=()
        )
        graph.add_edge(node, rotate, color="blue")
    elif option == "ZP":
        graph.remove_node(node)


def _sample_wire_fragment(graph, shots, backend, software):
    """Sample a renumbered fragment DAG on the chosen backend; return qubit-0-first counts.

    :param graph: A renumbered ``0..k-1`` subcircuit DAG.
    :param shots: Number of measurement shots.
    :param backend: Qibo platform name (ignored by non-Qibo backends).
    :param software: ``"qibo"``/``"qiskit"``/``"cudaq"``/``"pennylane"`` or a
        registered backend name.
    :return: ``{bitstring: count}`` with qubit-0-first keys, or ``{}`` for an empty
        fragment.
    """
    from Qdislib.core.backend.registry import build_execution_backend

    if graph is None or graph.number_of_nodes() == 0:
        return {}
    _, highest = _update_qubits(graph)
    extra = {"platform": backend} if backend else {}
    exec_backend = build_execution_backend(software, num_qubits=highest, **extra)
    return exec_backend.evaluate_samples(graph, highest, shots)


@task(returns=1)
def _run_wire_measure_fragment(
    sub, cuts, cut_qubits, basis_combo, shots, backend, software
):
    """Evaluate one measure ("top") fragment for a measurement-basis configuration.

    Rotates each cut qubit's ``O_i`` marker into the requested basis (X -> H,
    Y -> SDG + H, Z -> measure directly), samples the fragment, then extracts, for
    every ``P``/``I`` combination of its cut qubits, the joint expectation of those
    Pauli/identity factors with ``Z`` on the fragment's other (observable) qubits.

    :param sub: The measure fragment DAG (carries the ``O_i`` markers).
    :param cuts: The cut indices this fragment owns (sorted).
    :param cut_qubits: ``{cut_index: original qubit}`` for each owned cut.
    :param basis_combo: One of ``("X"|"Y"|"Z", ...)`` aligned with ``cuts``.
    :param shots: Number of measurement shots.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :return: ``{extraction_combo: expectation}`` where ``extraction_combo`` is a
        tuple over ``cuts`` of ``"P"``/``"I"`` (``"I"`` only when the basis is Z).
    """
    graph = sub.copy()
    used = sorted(
        {q for node in graph.nodes for q in graph.nodes[node].get("qubits", ())}
    )
    new_position = {q: i for i, q in enumerate(used)}
    cut_positions = [new_position[cut_qubits[cut]] for cut in cuts]

    for cut, basis in zip(cuts, basis_combo):
        # X/Y/Z basis == the XP/YP/ZP measure option (I is recovered from the Z
        # circuit's counts below, not a separate circuit).
        _apply_wire_measure_option(graph, cut, basis + "P")

    counts = _sample_wire_fragment(graph, shots, backend, software)
    total = sum(counts.values()) or 1
    width = len(next(iter(counts), "")) if counts else 0
    observable_positions = [p for p in range(width) if p not in cut_positions]

    options = [(["P", "I"] if basis == "Z" else ["P"]) for basis in basis_combo]
    results = {}
    for extraction in itertools.product(*options):
        selected = [pos for pos, ext in zip(cut_positions, extraction) if ext == "P"]
        value = 0.0
        for bitstring, count in counts.items():
            bits = [int(b) for b in bitstring]
            parity = 0
            for pos in observable_positions:
                parity ^= bits[pos]
            for pos in selected:
                parity ^= bits[pos]
            value += (-1 if parity else 1) * (count / total)
        results[extraction] = value
    return results


class _OptimizedWireMeasure:
    """A measure ("top") fragment subcircuit for the optimized three-step workflow.

    Unlike a plain subcircuit (one circuit -> one expectation value), evaluating
    this returns the whole ``{extraction_combo: value}`` table for its measurement
    basis -- i.e. one circuit yields every P/I observable, so a Z-basis measurement
    gives both ``<Z>`` and ``<I>``. That is what keeps the measure side at ``3**n``
    (not ``4**n``): :func:`expectation_value` detects the ``_is_optimized_wire_measure``
    marker and calls :meth:`evaluate` instead of the single-value path. Picklable
    (holds only the fragment DAG + metadata), so it survives a PyCOMPSs task boundary.
    """

    _is_optimized_wire_measure = True

    def __init__(self, sub, cuts, cut_qubits, basis_combo, software):
        self.sub = sub
        self.cuts = cuts
        self.cut_qubits = cut_qubits
        self.basis_combo = basis_combo
        self.software = software

    def evaluate(self, shots):
        """Return the ``{extraction_combo: value}`` table (a future under PyCOMPSs)."""
        return _run_wire_measure_fragment(
            self.sub,
            self.cuts,
            self.cut_qubits,
            self.basis_combo,
            shots,
            None,
            self.software,
        )

    def __repr__(self):
        return f"_OptimizedWireMeasure(cuts={self.cuts}, basis={self.basis_combo})"


@task(returns=1)
def _run_wire_prepare_fragment(sub, cuts, prep_combo, shots, backend, software):
    """Evaluate one prepare ("bottom") fragment for a preparation configuration.

    Sets each cut qubit's ``PS_i`` marker to prepare the requested state
    (|0> -> drop, |1> -> X, |+> -> H, |+i> -> H then S), samples the fragment and
    returns its global ``<Z...Z>`` parity. Also used (with empty ``cuts``) to fold
    in a fragment untouched by any cut.

    :param sub: The prepare fragment DAG (carries the ``PS_i`` markers).
    :param cuts: The cut indices this fragment owns (sorted); empty for a plain fragment.
    :param prep_combo: One of ``("0"|"1"|"+"|"+i", ...)`` aligned with ``cuts``.
    :param shots: Number of measurement shots.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :return: The fragment's ``<Z...Z>`` expectation.
    """
    graph = sub.copy()
    for cut, prep in zip(cuts, prep_combo):
        _apply_wire_prep_state(graph, cut, prep)

    counts = _sample_wire_fragment(graph, shots, backend, software)
    total = sum(counts.values()) or 1
    value = 0.0
    for bitstring, count in counts.items():
        parity = 1 if bitstring.count("1") % 2 == 0 else -1
        value += parity * (count / total)
    return value


def _optimized_wire_fragments(dag, wire_cuts, num_qubits):
    """Insert cut markers and split the circuit into measure / prepare / plain fragments.

    Shared by the all-in-one and three-step optimized paths.

    :param dag: The circuit DAG (mutated: cut markers inserted).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :return: ``(measure_frags, prepare_frags, plain_frags)`` where ``measure_frags``
        is ``[(sub, owned_cuts, {cut: qubit}), ...]``, ``prepare_frags`` is
        ``[(sub, owned_cuts), ...]`` and ``plain_frags`` is ``[sub, ...]``; or
        ``None`` when a fragment carries BOTH a measure and a prepare marker (two
        cuts in series on one wire -> not fragment-separable) or a cut lacks a marker.
    """
    n_cuts = len(wire_cuts)
    _insert_wire_cut_markers(dag, wire_cuts, num_qubits)
    graph = _remove_red_edges(dag)

    measure_frags, prepare_frags, plain_frags = [], [], []
    measured, prepared = set(), set()
    for component in networkx.connected_components(graph.to_undirected()):
        sub = graph.subgraph(component).copy()
        measures = sorted(i for i in range(1, n_cuts + 1) if f"O_{i}" in sub.nodes)
        preps = sorted(i for i in range(1, n_cuts + 1) if f"PS_{i}" in sub.nodes)
        if measures and preps:
            return None  # mixed fragment -> not separable
        if measures:
            cut_qubits = {i: sub.nodes[f"O_{i}"]["qubits"][0] for i in measures}
            measure_frags.append((sub, measures, cut_qubits))
            measured.update(measures)
        elif preps:
            prepare_frags.append((sub, preps))
            prepared.update(preps)
        else:
            plain_frags.append(sub)

    if len(measured) != n_cuts or len(prepared) != n_cuts:
        return None
    return measure_frags, prepare_frags, plain_frags


@task(
    returns=1,
    measure_flat=COLLECTION_IN,
    prepare_flat=COLLECTION_IN,
    plain_flat=COLLECTION_IN,
)
def _reconstruct_optimized_wire_cut(
    measure_flat,
    measure_keys,
    measure_cuts,
    prepare_flat,
    prepare_keys,
    prepare_cuts,
    plain_flat,
    n_cuts,
):
    """Classical ``8**n`` signed sum over the cached fragment values (a COMPSs task).

    The fragment futures feed straight in through ``COLLECTION_IN`` (task -> task,
    no master barrier before reconstruction), so the whole optimized call is waited
    on exactly once -- like the standard method's ``_change_sign`` reconstruction.
    Rebuilds the per-fragment value tables from the flat futures + their keys, then
    evaluates the base-8 sum (sign flip on ``{3, 5, 7}``, prefactor ``1/2**n``),
    expanding |-> / |-i> as linear combinations of the prepared states.

    :param measure_flat: Measure-fragment results, one ``{extraction: value}`` dict
        per (fragment, basis-combo), in ``measure_keys`` order.
    :param measure_keys: ``(fragment_index, basis_combo)`` for each ``measure_flat`` entry.
    :param measure_cuts: The owned cut indices of each measure fragment.
    :param prepare_flat: Prepare-fragment ``<Z...Z>`` values, in ``prepare_keys`` order.
    :param prepare_keys: ``(fragment_index, prep_combo)`` for each ``prepare_flat`` entry.
    :param prepare_cuts: The owned cut indices of each prepare fragment.
    :param plain_flat: ``<Z...Z>`` of each cut-free fragment (a constant factor).
    :param n_cuts: Number of wire cuts.
    :return: The reconstructed expectation value.
    """
    measure_values = {}
    for (frag_idx, basis_combo), value in zip(measure_keys, measure_flat):
        measure_values.setdefault(frag_idx, {})[basis_combo] = value
    prepare_values = {}
    for (frag_idx, prep_combo), value in zip(prepare_keys, prepare_flat):
        prepare_values.setdefault(frag_idx, {})[prep_combo] = value
    plain_factor = 1.0
    for value in plain_flat:
        plain_factor *= value

    total = 0.0
    for index in range(8**n_cuts):
        settings = [(index // (8 ** (n_cuts - 1 - k))) % 8 for k in range(n_cuts)]
        sign = -1.0 if sum(s in _WIRE_SIGN_FLIP for s in settings) % 2 else 1.0
        term = sign

        for frag_idx, cuts in enumerate(measure_cuts):
            basis_combo = tuple(_WIRE_MEASURE[settings[cut - 1]][0] for cut in cuts)
            extraction = tuple(_WIRE_MEASURE[settings[cut - 1]][1] for cut in cuts)
            term *= measure_values[frag_idx][basis_combo][extraction]

        for frag_idx, cuts in enumerate(prepare_cuts):
            expansions = [_WIRE_PREP[settings[cut - 1]] for cut in cuts]
            fragment_value = 0.0
            for choice in itertools.product(*expansions):
                coefficient = 1.0
                base_combo = []
                for weight, base_state in choice:
                    coefficient *= weight
                    base_combo.append(base_state)
                fragment_value += (
                    coefficient * prepare_values[frag_idx][tuple(base_combo)]
                )
            term *= fragment_value

        total += term * plain_factor

    return total / (2**n_cuts)


def _wire_cutting_optimized(dag, wire_cuts, num_qubits, shots, backend, software):
    """Fragment-separated wire-cut reconstruction (``3**n + 4**n`` subcircuits).

    Splits the cut circuit into measure fragments (own cuts via their ``O_i``
    markers) and prepare fragments (via ``PS_i``); a fragment carrying both kinds
    of marker (e.g. two cuts in series on one wire) is not separable this way, so
    this returns ``None`` and the caller uses the standard base-8 method. Each
    measure fragment is run for its ``3**(#cuts)`` basis configurations and each
    prepare fragment for its ``4**(#cuts)`` preparations; the ``8**n`` signed sum
    (base-8, sign flip on ``{3, 5, 7}``, prefactor ``1/2**n``) is then evaluated
    classically over the cached fragment values, with |-> / |-i> expanded as
    linear combinations. Backend-agnostic; each fragment run is a PyCOMPSs task.

    :param dag: The circuit DAG (mutated: cut markers inserted).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :param shots: Number of measurement shots per fragment.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :return: The reconstructed expectation value, or ``None`` if the cut topology
        is not fragment-separable (mixed measure/prepare fragment).
    """
    n_cuts = len(wire_cuts)
    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, plain_frags = fragments

    # Fan out every fragment execution (PyCOMPSs tasks), flattened with their keys.
    measure_flat, measure_keys = [], []
    for frag_idx, (sub, cuts, cut_qubits) in enumerate(measure_frags):
        for combo in itertools.product(_MEASURE_BASES, repeat=len(cuts)):
            measure_flat.append(
                _run_wire_measure_fragment(
                    sub, cuts, cut_qubits, combo, shots, backend, software
                )
            )
            measure_keys.append((frag_idx, combo))
    prepare_flat, prepare_keys = [], []
    for frag_idx, (sub, cuts) in enumerate(prepare_frags):
        for combo in itertools.product(_PREP_STATES, repeat=len(cuts)):
            prepare_flat.append(
                _run_wire_prepare_fragment(sub, cuts, combo, shots, backend, software)
            )
            prepare_keys.append((frag_idx, combo))
    plain_flat = [
        _run_wire_prepare_fragment(sub, [], (), shots, backend, software)
        for sub in plain_frags
    ]

    # Feed the fragment futures straight into the reconstruction task (task -> task,
    # no master barrier); the caller waits once on the returned value.
    return _reconstruct_optimized_wire_cut(
        measure_flat,
        measure_keys,
        [cuts for _sub, cuts, _q in measure_frags],
        prepare_flat,
        prepare_keys,
        [cuts for _sub, cuts in prepare_frags],
        plain_flat,
        n_cuts,
    )


@task(returns=1, lst=COLLECTION_IN)
def _sum_results(lst: typing.List[int]) -> int:
    """Calculate the sum of all results.

    :param lst: List of partial results.
    :return: Sum of all partial results.
    """
    return sum(lst)


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
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
    # A qubit can be cut-measured MORE THAN ONCE (several gate cuts on the same wire,
    # e.g. a general gate decomposed into 3 CZs), so collect every collapse gate per
    # qubit -- a plain {qubit: gate} dict would keep only the last and drop the others'
    # signs.
    collapse = {}
    for gate in qc.queue:
        if isinstance(gate, gates.M) and getattr(gate, "collapse", False):
            collapse.setdefault(gate.qubits[0], []).append(gate)
    qc.add(gates.M(*range(qc.nqubits)))
    final = np.asarray(qc(nshots=shots).samples())  # (shots, nqubits), qubit-0-first
    contribution = np.ones(len(final))
    for q, obs in enumerate(observables):
        if q in prep_qubits and q in collapse:
            # Measured cut qubit: each (I +/- Z)/2 projector contributes the
            # decomposition sign -(-1)**m from its cut-point outcome m for *any*
            # observable; Z then multiplies by the final measurement.
            for cut_gate in collapse[q]:
                m = np.asarray(cut_gate.result.samples()).ravel()
                contribution *= -((-1.0) ** m)
            if obs == "Z":
                contribution *= (-1.0) ** final[:, q]
            elif obs != "I":
                raise ValueError(f"Unsupported observable {obs}")
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

    if prep_qubits is None:
        # 3-step subcircuit workflow: prep qubits aren't threaded in, so recover
        # them from the subcircuit's collapse (gate-cut) measurements. Lets the
        # workflow use the decoupled reconstruction like the all-in-one path.
        prep_qubits = [
            g.qubits[0]
            for g in qc.queue
            if isinstance(g, gates.M) and getattr(g, "collapse", False)
        ]

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

    Mirrors the Qiskit/Qibo decoupling. The kernel is executed with ``cudaq.run``
    (not ``cudaq.sample``, which does not collapse mid-circuit measurements) and
    returns, per shot, ``[cut-point outcomes in mz-call order] + [final measurement
    of q0..q(n-1)]``. For a measured cut qubit the measure-and-prepare sign is
    recovered from its cut-point outcome ``m`` as ``-(-1)**m`` for any observable
    (``Z`` then multiplies by the final measurement).

    :param kernel: Kernel built with ``explicit_final=True`` (a ``cudaq.run`` kernel).
    :param observables: Per-qubit Pauli string over ``{I, Z}`` (qubit-0-first).
    :param shots: Number of shots.
    :param prep_qubits: Renumbered qubits measured for the cut.
    :param measured_order: Cut-measured qubits in ``mz``-call order.
    :param num_qubits: Qubit count (width of the final-measurement segment).
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    records = cudaq.run(kernel, shots_count=shots)
    n_mid = len(measured_order)
    last_pos = {qb: i for i, qb in enumerate(measured_order)}
    total = 0.0
    drawn = 0
    for record in records:
        drawn += 1
        contribution = 1
        # record is an int; bit k = outcome k of [cut-points..., finals...].
        for q, obs in enumerate(observables):
            if q in prep_qubits and q in last_pos:
                # Measured cut qubit: decomposition sign -(-1)**m from the
                # cut-point outcome m for *any* observable; Z then multiplies by
                # the final measurement.
                contribution *= -((-1) ** ((record >> last_pos[q]) & 1))
                if obs == "Z":
                    contribution *= (-1) ** ((record >> (n_mid + q)) & 1)
                elif obs != "I":
                    raise ValueError(f"Unsupported observable {obs}")
            elif obs == "Z":
                contribution *= (-1) ** ((record >> (n_mid + q)) & 1)
            elif obs == "I":
                continue
            else:
                raise ValueError(f"Unsupported observable {obs}")
        total += contribution
    return total / (drawn or 1)


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
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

    if prep_qubits is None and is_legacy_dag:
        # 3-step subcircuit workflow: recover the gate-cut measured qubits from the
        # DAG so the workflow uses the decoupled reconstruction like the all-in-one
        # path (its measure nodes collapse-and-continue, so a plain read is wrong).
        dag_component = payload[0]
        prep_qubits = sorted(
            {
                (dag_component.nodes[n].get("qubits") or (None,))[0]
                for n in dag_component
                if dag_component.nodes[n].get("gate") == "measure"
            }
            - {None}
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


def _counts_pennylane(circuit, shots=1024, device="default.qubit", seed=None):
    """Sample a PennyLane subcircuit and return its measurement counts.

    Builds a QNode on ``device`` that replays the subcircuit's gates and returns
    ``qml.counts`` over all wires -- PennyLane's own execution path.

    :param circuit: A :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit`.
    :param shots: Number of measurement shots.
    :param device: PennyLane device short name.
    :param seed: Optional RNG seed.
    :return: A ``{bitstring: count}`` dict with **qubit-0-first** keys.
    """
    ops, num_qubits = circuit.ops, circuit.num_qubits

    def qfunc():
        for gate, qubits, params in ops:
            _apply_pennylane_gate(gate, qubits, params)
        return qml.counts(wires=range(num_qubits))

    qnode = _build_pennylane_qnode(qfunc, num_qubits, shots, device, seed)
    raw = qnode()
    return {str(key): int(value) for key, value in raw.items()}


def _expectation_from_counts_pennylane(freq, observables, shots=1024):
    """Compute ``<O>`` for a Z/I Pauli string from PennyLane measurement counts.

    PennyLane ``counts`` keys are qubit-0-first bit strings, so they align
    directly with the observable string (like Qibo/CUDA-Q).

    :param freq: A ``{bitstring: count}`` dict as returned by :func:`_counts_pennylane`.
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
        expectation_value += contribution * (value / shots)
    return expectation_value


def _expectation_from_samples_pennylane(
    circuit, observables, shots, prep_qubits, device="default.qubit", seed=None
):
    """Decoupled ``<O>`` for PennyLane gate-cut subcircuits with identity on a cut qubit.

    Mirrors :func:`_expectation_from_samples_qibo`. A cut qubit measured for the
    decomposition (a ``MEASURE`` node -> ``qml.measure``, collapse-and-continue)
    carries a sign normally read off its final measurement; when the user observes
    it with identity, the sign is recovered from the *cut-point* outcome ``m`` as
    ``-(-1)**m`` instead. This needs per-shot samples (final measurements and the
    mid-circuit ``MeasurementValue``, kept correlated by one QNode execution).

    :param circuit: A :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit`.
    :param observables: Per-qubit Pauli string over ``{I, Z}`` (qubit-0-first).
    :param shots: Number of measurement shots.
    :param prep_qubits: Renumbered qubits measured for the cut.
    :param device: PennyLane device short name.
    :param seed: Optional RNG seed.
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    ops, num_qubits = circuit.ops, circuit.num_qubits
    # Collapse-measured qubits, in circuit order (matches the sample-return order).
    collapse_order = []
    for gate, qubits, _params in ops:
        if gate.upper() == "MEASURE":
            for qubit in qubits:
                if qubit not in collapse_order:
                    collapse_order.append(qubit)

    def qfunc():
        mids = {}
        for gate, qubits, params in ops:
            value = _apply_pennylane_gate(gate, qubits, params)
            if value is not None:  # a MEASURE returns its MeasurementValue
                mids[qubits[0]] = value
        # Finals for every qubit, then each cut-point outcome -- one execution keeps
        # them correlated per shot (built here so the MeasurementValues are in scope).
        return [qml.sample(wires=range(num_qubits))] + [
            qml.sample(mids[qubit]) for qubit in collapse_order
        ]

    qnode = _build_pennylane_qnode(qfunc, num_qubits, shots, device, seed)
    out = qnode()
    final = np.asarray(out[0]).reshape(shots, num_qubits)  # (shots, n), qubit-0-first
    mid_samples = {
        qubit: np.asarray(out[1 + i]).ravel() for i, qubit in enumerate(collapse_order)
    }

    contribution = np.ones(shots)
    for qubit, obs in enumerate(observables):
        if qubit in prep_qubits and qubit in mid_samples:
            # Measured cut qubit: the (I +/- Z)/2 projector contributes the
            # decomposition sign -(-1)**m from the cut-point outcome m for *any*
            # observable; Z then multiplies by the final measurement.
            contribution *= -((-1.0) ** mid_samples[qubit])
            if obs == "Z":
                contribution *= (-1.0) ** final[:, qubit]
            elif obs != "I":
                raise ValueError(f"Unsupported observable {obs}")
        elif obs == "Z":
            contribution *= (-1.0) ** final[:, qubit]
        elif obs == "I":
            continue
        else:
            raise ValueError(f"Unsupported observable {obs}")
    return float(contribution.mean())


def _expec_value_pennylane(
    subcirc, shots=1024, prep_qubits=None, device="default.qubit"
):
    """Execute a PennyLane subcircuit and return its expectation value.

    :param subcirc: A :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit`, or
        a ``[PennylaneCircuit, observables]`` pair.
    :param shots: Number of measurement shots.
    :param prep_qubits: Renumbered cut qubits measured for a gate cut; when any is
        observed with identity, its decomposition sign is recovered from the
        cut-point measurement (see :func:`_expectation_from_samples_pennylane`).
        Auto-detected from the subcircuit's ``MEASURE`` nodes when ``None`` (the
        three-step workflow), matching the Qibo/CUDA-Q evaluators.
    :param device: PennyLane device short name.
    :return: The circuit expectation value.
    """
    if subcirc is None:
        return None

    if (
        isinstance(subcirc, (list, tuple))
        and len(subcirc) == 2
        and isinstance(subcirc[1], str)
    ):
        circuit, observables = subcirc[0], subcirc[1]
    else:
        circuit, observables = subcirc, "Z" * subcirc.num_qubits

    if prep_qubits is None:
        prep_qubits = sorted(
            {
                qubit
                for gate, qubits, _params in circuit.ops
                if gate.upper() == "MEASURE"
                for qubit in qubits
            }
        )

    if prep_qubits:
        return _expectation_from_samples_pennylane(
            circuit, observables, shots, set(prep_qubits), device
        )

    freq = _counts_pennylane(circuit, shots, device)
    return _expectation_from_counts_pennylane(freq, observables, shots)


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
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
    """Signed product of subcircuit values for one wire-cut term (base 8, {3, 5, 7}).

    See :func:`Qdislib.core.cutting_algorithms._reconstruction.signed_product`.
    """
    # Resolved and reduced on the worker (not the master), so the subcircuit
    # futures feed straight into this task -- no barrier before the sign step.
    from Qdislib.core.cutting_algorithms._reconstruction import signed_product

    return signed_product(expectation_value, index, 8, {3, 5, 7})


def _make_optimized_wire_reconstructor(
    n_cuts, measure_index, prepare_index, plain_positions, measure_cuts, prepare_cuts
):
    """Build the ``reconstruct`` callback for the optimized three-step subcircuits.

    Returns a function with the :class:`~Qdislib.utils.backends.CutSubcircuits`
    reconstructor signature that maps the flat list of per-subcircuit expectation
    values back to the fragment value tables (via the captured index maps) and
    evaluates the base-8 signed sum, expanding |-> / |-i> as linear combinations.
    ``number_cuts`` / ``number_components`` are ignored -- the layout is captured here.
    """

    def reconstruct(expected_values, number_cuts=None, number_components=None):
        from Qdislib.utils.circuit_cache import _resolve_value

        # Under PyCOMPSs the per-subcircuit values are futures; resolve them before
        # the master-side classical sum (the standard reconstructor feeds them into
        # @tasks instead, but the optimized sum is cheap and done here).
        expected_values = compss_wait_on(expected_values)
        values = [_resolve_value(v) for v in expected_values]
        measure_values = {}
        for (frag_idx, combo), position in measure_index.items():
            measure_values.setdefault(frag_idx, {})[combo] = values[position]
        prepare_values = {}
        for (frag_idx, combo), position in prepare_index.items():
            prepare_values.setdefault(frag_idx, {})[combo] = values[position]
        plain_factor = 1.0
        for position in plain_positions:
            plain_factor *= values[position]

        total = 0.0
        for index in range(8**n_cuts):
            settings = [(index // (8 ** (n_cuts - 1 - k))) % 8 for k in range(n_cuts)]
            sign = -1.0 if sum(s in _WIRE_SIGN_FLIP for s in settings) % 2 else 1.0
            term = sign
            for frag_idx, cuts in enumerate(measure_cuts):
                # One Z-basis measure circuit yields both <Z> and <I> (dict of
                # extractions), so the measure side is 3**n, matching the all-in-one.
                basis_combo = tuple(_WIRE_MEASURE[settings[cut - 1]][0] for cut in cuts)
                extraction = tuple(_WIRE_MEASURE[settings[cut - 1]][1] for cut in cuts)
                term *= measure_values[frag_idx][basis_combo][extraction]
            for frag_idx, cuts in enumerate(prepare_cuts):
                expansions = [_WIRE_PREP[settings[cut - 1]] for cut in cuts]
                fragment_value = 0.0
                for choice in itertools.product(*expansions):
                    coefficient = 1.0
                    base_combo = []
                    for weight, base_state in choice:
                        coefficient *= weight
                        base_combo.append(base_state)
                    fragment_value += (
                        coefficient * prepare_values[frag_idx][tuple(base_combo)]
                    )
                term *= fragment_value
            total += term * plain_factor
        return total / (2**n_cuts)

    return reconstruct


def _optimized_wire_cutting_subcircuits(dag, cut, num_qubits, software):
    """Three-step counterpart of :func:`_wire_cutting_optimized`.

    Returns a :class:`~Qdislib.utils.backends.CutSubcircuits` holding the reduced
    fragment subcircuits -- ``3**n`` measure probes (each an
    :class:`_OptimizedWireMeasure` returning all its P/I observables from one
    circuit) and ``4**n`` base preparations, i.e. ``3**n + 4**n`` instead of
    ``8**n * components`` -- whose ``reconstruct`` does the fragment-separated base-8
    sum. Every subcircuit is a normal ``expectation_value`` target, so the standard
    three-step loop is unchanged. Returns ``None`` when the cut topology is not
    fragment-separable, so the caller falls back to the standard generator.
    """
    if software not in ("qibo", "qiskit", "cudaq", "pennylane"):
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            "Expected one of ['qibo', 'qiskit', 'cudaq', 'pennylane']."
        )
    wire_cuts = _expand_wire_cuts(dag, cut)
    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, plain_frags = fragments

    subcircuits = []
    measure_index, prepare_index, plain_positions = {}, {}, []
    for frag_idx, (sub, cuts, cut_qubits) in enumerate(measure_frags):
        # One _OptimizedWireMeasure per basis combo (3**n) -- evaluating it returns
        # every P/I observable of that circuit, so <Z> and <I> share one circuit.
        for basis_combo in itertools.product(_MEASURE_BASES, repeat=len(cuts)):
            measure_index[(frag_idx, basis_combo)] = len(subcircuits)
            subcircuits.append(
                _OptimizedWireMeasure(
                    sub.copy(), cuts, cut_qubits, basis_combo, software
                )
            )
    for frag_idx, (sub, cuts) in enumerate(prepare_frags):
        for prep_combo in itertools.product(_PREP_STATES, repeat=len(cuts)):
            graph = sub.copy()
            for cut_index, state in zip(cuts, prep_combo):
                _apply_wire_prep_state(graph, cut_index, state)
            prepare_index[(frag_idx, prep_combo)] = len(subcircuits)
            subcircuits.append(_make_subcircuit_task(graph, software))
    for sub in plain_frags:
        plain_positions.append(len(subcircuits))
        subcircuits.append(_make_subcircuit_task(sub.copy(), software))

    reconstructor = _make_optimized_wire_reconstructor(
        len(wire_cuts),
        measure_index,
        prepare_index,
        plain_positions,
        [cuts for _s, cuts, _q in measure_frags],
        [cuts for _s, cuts in prepare_frags],
    )
    return CutSubcircuits(
        subcircuits,
        number_cuts=len(wire_cuts),
        number_components=len(subcircuits),
        reconstructor=reconstructor,
    )


def wire_cutting_subcircuits(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    observables=None,
    software="qiskit",
    optimized: bool = False,
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

    if optimized and cut:
        # Fragment-separated subcircuits (~4**n + 4**n instead of 8**n * components);
        # the returned object's reconstruct() does the classical base-8 sum. Falls
        # back to the standard generator for non-separable topologies or on error.
        try:
            optimized_subcircuits = _optimized_wire_cutting_subcircuits(
                dag.copy(), cut, _max_qubit(dag), software
            )
            if optimized_subcircuits is not None:
                return optimized_subcircuits
            reason = "not applicable to this cut topology"
        except UnsupportedBackendError:
            raise
        except (
            Exception
        ) as exc:  # never let the optimization break subcircuit generation
            reason = f"failed ({exc})"
        logger.info(
            "Optimized wire-cut subcircuits %s; using the standard generator.", reason
        )

    # Serial -> native backend circuit; PyCOMPSs -> serializable Subcircuit.
    dispatch_map = {
        "qibo": lambda d, n: make_subcircuit(d, n, "qibo"),
        "qiskit": lambda d, n: make_subcircuit(d, n, "qiskit"),
        "cudaq": lambda d, n: make_subcircuit(d, n, "cudaq"),
        "pennylane": lambda d, n: make_subcircuit(d, n, "pennylane"),
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

    if software not in ("qibo", "qiskit", "cudaq", "pennylane"):
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            "Expected one of ['qibo', 'qiskit', 'cudaq', 'pennylane']."
        )

    # Each subcircuit is *created* in its own task (_make_subcircuit_task), separate
    # from the per-subcircuit evaluation task -- so under PyCOMPSs these are
    # Subcircuit futures fed straight into expectation_value's eval task (task -> task).
    subcircuits = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag.copy())

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected()
        )
        graph_components = [networkx.DiGraph() for _ in range(num_components)]

        _generate_subcircuits_wire_cutting(
            copy_graph,
            num_qubits + len(wire_cuts),
            index,
            wire_cuts,
            graph_components,
        )

        for s in graph_components:
            subcircuits.append(_make_subcircuit_task(s, software))
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
    if number_cuts == 0:
        # No cuts: the reconstruction is just the product of the (single group of)
        # subcircuit values, computed on the master -- so sync + resolve here.
        from Qdislib.utils.circuit_cache import _resolve_value

        expected_values = compss_wait_on(expected_values)
        expected_values = [_resolve_value(v) for v in expected_values]
        return math.prod([v for v in expected_values if v is not None])

    # Feed each configuration's subcircuit futures straight into the sign task
    # (task -> task, no master barrier between _build_and_run_subcircuit /
    # _evaluate_subcircuit and _change_sign); the sign task resolves cache placeholders.
    groups = [
        expected_values[i : i + number_components]
        for i in range(0, len(expected_values), number_components)
    ]
    results = [_change_sign(group, index) for index, group in enumerate(groups)]
    results = compss_wait_on(results)
    return 1 / (2**number_cuts) * sum(results)
