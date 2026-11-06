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
    "wire_cutting_subcircuits",
    "wire_cutting_subcircuit_reconstruction",
]

import functools
import itertools
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

from Qdislib.utils.sdks import sdk_names as _sdk_names
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
    CombinedCutSubcircuits,
    Subcircuit,
)
from Qdislib.utils.exceptions import UnsupportedBackendError, InvalidCutError
from Qdislib.core.cutting_algorithms._reconstruction import next_seed
from Qdislib.core.cutting_algorithms._api import (
    standard_arguments,
    resolve_backend,
    checked_observable,
)


def _warn_not_optimized(what, reason):
    """``optimized=True`` was asked for and the standard method runs instead: the
    same result, but ``8**n`` configurations per fragment set."""
    from Qdislib.utils.logging_config import warn_once

    warn_once(
        f"optimized=True: optimized wire {what} {reason}, so this runs the standard "
        "method (the same result, at up to 8**n configurations for n cuts)."
    )


def _wire_cutting(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    sync: bool = True,
    observables: str = None,
    shots: int = 1024,
    backend: str = "numpy",
    qpu: str = None,
    software: str = None,
    cache=None,
    cache_mode="expectation",
    optimized: bool = False,
    joint: bool = True,
    seed=None,
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
    software : str, optional
        Backend used to build and execute the subcircuits: "qibo", "qiskit" or "cudaq".
        By default, the SDK the circuit was written in ("qibo" for a DAG).
    cache : CircuitCache, optional
        Enable the semantic subcircuit cache (see
        :func:`~Qdislib.utils.circuit_cache.default_circuit_cache`). When given,
        evaluation runs the CPU-simulation subcircuit workflow (GPU/QPU dispatch
        is bypassed). ``None`` (default) disables caching.
    cache_mode : str or CacheStrategy, default="expectation"
        Cache granularity; ``"expectation"``, ``"counts"``, or a custom strategy.
    optimized : bool, default=False
        Use the fragment-separated reconstruction: splits the cut circuit into
        measure fragments (``3**(#cuts)`` subcircuits each), prepare fragments
        (``4**(#cuts)`` each), and -- for chain topologies like two cuts in series
        on one wire, where a fragment owns both kinds of marker -- mixed fragments
        (``3**(#measure-cuts) * 4**(#prepare-cuts)`` each), combining their cached
        values classically (identical result, always fewer circuits than the
        standard method's ``8**n``, though never as few as a fully separable cut
        of the same size). Falls back to the standard method only if a cut's
        marker is missing entirely (a bookkeeping inconsistency, not a topology
        issue).

    joint : bool, default=True
        With ``optimized=True``, measure the cut wires each fragment sends together
        instead of one by one: ``2**k + 1`` measurement circuits for ``k`` wires
        leaving a fragment instead of ``3**k`` (5 instead of 9 at ``k = 2``, 9
        instead of 27 at ``k = 3``), reading every Pauli of the wires exactly once.
        Same value; usually fewer shots for the same error, not always for
        ``k = 2``. Some of those circuits add CZ gates between the fragment's cut
        qubits: pass ``False`` to measure each wire on its own (no entangling gates),
        e.g. on hardware with costly two-qubit gates. Ignored without ``optimized``.

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

        qc = QuantumCircuit(2)
        qc.h(0)     # Gate "H_1"
        qc.cx(0, 1) # Gate "CX_2"
        qc.measure_all()

        cuts = [("H_1","CX_2")]  # List of gate tuples
        reconstruction = _wire_cutting(qc, cuts, shots=2048, software="qiskit")
    """
    if software is None:
        # Run the subcircuits in the circuit's own SDK: a Qiskit circuit's gate
        # names (cx, ...) are not Qibo gate classes.
        from Qdislib.utils.circuit_cache import _detect_software

        software = _detect_software(rand_qc) or "qibo"

    # A weighted sum of Pauli observables (a Hamiltonian) is evaluated term by
    # term: <sum_i c_i P_i> = sum_i c_i <P_i>. Each term reuses the single
    # Pauli-string path below. A plain string (or None) skips this and runs once.
    if observables is not None and not isinstance(observables, str):
        # Every term's cut is submitted without waiting, so all of them run
        # together, and one task folds the term futures.
        from Qdislib.core.cutting_algorithms.gate_cutting import _weighted_sum

        coeffs, values = [], []
        for k, (coeff, pauli) in enumerate(normalize_observables(observables)):
            base = rand_qc.copy() if isinstance(rand_qc, networkx.Graph) else rand_qc
            coeffs.append(coeff)
            values.append(
                _wire_cutting(
                    base,
                    cut,
                    sync=False,
                    observables=pauli,
                    shots=shots,
                    backend=backend,
                    qpu=qpu,
                    software=software,
                    cache=cache,
                    cache_mode=cache_mode,
                    optimized=optimized,
                    joint=joint,
                    seed=next_seed(seed, "term", k),
                )
            )
        total = _weighted_sum(coeffs, values)
        return compss_wait_on(total) if sync else total

    # Semantic cache: reuse the recommended subcircuit workflow, which routes
    # each subcircuit through expectation_value(cache=...). This runs the CPU
    # simulation path, so GPU/QPU dispatch is bypassed when caching is enabled.
    if cache is not None:
        from Qdislib.utils.backends import expectation_value

        subcircuits = _wire_cutting_subcircuits(
            rand_qc, cut, observables=observables, software=software
        )
        values = [
            expectation_value(
                s,
                shots=shots,
                cache=cache,
                cache_mode=cache_mode,
                seed=next_seed(seed, "cached", i),
            )
            for i, s in enumerate(subcircuits)
        ]
        return subcircuits.reconstruct(values)

    dag = convert_to_dag(rand_qc, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    _validate_wire_cuts_partition(dag, _expand_wire_cuts(dag, cut), _max_qubit(dag))

    if optimized and cut:
        # Fragment-separated reconstruction: sum_f 3**i_f * 4**j_f subcircuits
        # instead of 8**n (measure, prepare, and -- for chain topologies -- mixed
        # fragments owning both kinds of marker all handled). Falls back to the
        # standard method only if anything goes wrong, so optimized=True never
        # changes the result.
        try:
            optimized_value = _wire_cutting_optimized(
                dag.copy(),
                _expand_wire_cuts(dag, cut),
                _max_qubit(dag),
                shots,
                backend,
                software,
                joint=joint,
                seed=seed,
            )
            # The optimized path returns a single reconstruction future (task ->
            # task, no master barrier); wait on it once here.
            if optimized_value is not None:
                return compss_wait_on(optimized_value) if sync else optimized_value
            reason = "not applicable to this cut topology"
        except Exception as exc:  # never let the optimization break the computation
            reason = f"failed ({exc})"
        _warn_not_optimized("cutting", reason)

    if networkx.number_connected_components(dag.to_undirected(as_view=True)) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected(as_view=True))
        ]
        results = []
        total_wire_cuts = 0
        for piece, s in enumerate(S):
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
                    seed=next_seed(seed, "piece", piece),
                )
                graphs = _sum_results(graphs)
                results.append(graphs)
            else:
                results.append(
                    _build_and_run_subcircuit(
                        s,
                        software,
                        shots,
                        backend,
                        qpu,
                        next_seed(seed, "piece", piece),
                    )
                )
        from Qdislib.core.cutting_algorithms._reconstruction import scaled_product

        final_recons = scaled_product(results, 1 / (2**total_wire_cuts))
        return compss_wait_on(final_recons) if sync else final_recons
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
                seed=seed,
            )

            from Qdislib.core.cutting_algorithms._reconstruction import scaled_sum

            final_recons = scaled_sum(results, 1 / (2 ** len(cut_wires)))
            return compss_wait_on(final_recons) if sync else final_recons
        else:
            final_recons = _build_and_run_subcircuit(
                dag, software, shots, backend, qpu, seed
            )
        return compss_wait_on(final_recons) if sync else final_recons


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
    placed = _place_wire_cuts(dag, wire_cuts, num_qubits)

    # Cut-source line -> its 1-based cut index. The base-8 digit for that cut sets
    # the measure basis: identity for digits 0/1 (observe the cut-source with I), else
    # X/Y/Z (a Z decomposition factor). See :func:`_reconstruct_wire_cut_distribution`.
    cut_source_index, prepared_to_orig = _wire_cut_line_maps(placed, num_qubits)

    components = []
    per_config_counts = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag)
        num_components = networkx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
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
    mixed_flat=COLLECTION_IN,
    plain_flat=COLLECTION_IN,
)
def _reconstruct_optimized_wire_distribution_task(
    measure_flat,
    measure_keys,
    measure_cuts,
    prepare_flat,
    prepare_keys,
    prepare_cuts,
    mixed_flat,
    mixed_keys,
    mixed_cuts,
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
    input state). Mixed fragments (own both kinds of marker, e.g. a chain of cuts
    in series) are looked up jointly by ``(basis_combo, prep_base_combo)`` -- one
    physical circuit per pair -- and still get the |-> / |-i> expansion on their
    prepare-type cuts. Reuses :func:`_reconstruct_wire_cut_distribution` unchanged.
    """
    # The wire Pauli a cut bit stands for, per base-8 digit (0/1 never set it).
    basis_of = {0: "Z", 1: "Z", 2: "X", 3: "X", 4: "Y", 5: "Y", 6: "Z", 7: "Z"}

    # Per fragment, one {(output_mask, pauli): value} table merged over its
    # measurement settings (per-wire or joint), and its qubit order.
    measure_tables, measure_original = {}, {}
    for (frag_idx, setting), (corr, original) in zip(measure_keys, measure_flat):
        cuts, cut_qubits = measure_cuts[frag_idx]
        readout = _setting_readout(setting, len(cuts))
        measure_tables.setdefault(frag_idx, []).append(
            _corr_pauli_table(corr, original, cuts, cut_qubits, readout)
        )
        measure_original[frag_idx] = original
    measure_tables = {f: _merge_tables(t) for f, t in measure_tables.items()}
    prepare_corr = {}
    for (frag_idx, prep_combo), value in zip(prepare_keys, prepare_flat):
        prepare_corr.setdefault(frag_idx, {})[prep_combo] = value
    mixed_tables, mixed_original = {}, {}
    for (frag_idx, setting, prep_combo), (corr, original) in zip(
        mixed_keys, mixed_flat
    ):
        m_cuts, cut_qubits, _p = mixed_cuts[frag_idx]
        readout = _setting_readout(setting, len(m_cuts))
        mixed_tables.setdefault(frag_idx, {}).setdefault(prep_combo, []).append(
            _corr_pauli_table(corr, original, m_cuts, cut_qubits, readout)
        )
        mixed_original[frag_idx] = original
    mixed_tables = {
        f: {prep: _merge_tables(t) for prep, t in by_prep.items()}
        for f, by_prep in mixed_tables.items()
    }
    rebuilt = {}  # (kind, fragment, prep, letters) -> per-wire-shaped corr

    def measure_corr(frag_idx, letters):
        key = ("m", frag_idx, letters)
        if key not in rebuilt:
            cuts, cut_qubits = measure_cuts[frag_idx]
            rebuilt[key] = _corr_from_pauli_table(
                measure_tables[frag_idx],
                measure_original[frag_idx],
                cuts,
                cut_qubits,
                letters,
            )
        return rebuilt[key], measure_original[frag_idx]

    def mixed_corr(frag_idx, letters, prep_combo):
        key = ("x", frag_idx, letters, prep_combo)
        if key not in rebuilt:
            m_cuts, cut_qubits, _p = mixed_cuts[frag_idx]
            rebuilt[key] = _corr_from_pauli_table(
                mixed_tables[frag_idx][prep_combo],
                mixed_original[frag_idx],
                m_cuts,
                cut_qubits,
                letters,
            )
        return rebuilt[key], mixed_original[frag_idx]

    configs = []
    for index in range(8**num_cuts):
        settings = [(index // (8 ** (num_cuts - 1 - k))) % 8 for k in range(num_cuts)]
        per_component = []
        for frag_idx, (cuts, _q) in enumerate(measure_cuts):
            letters = tuple(basis_of[settings[cut - 1]] for cut in cuts)
            per_component.append(measure_corr(frag_idx, letters))
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
        for frag_idx, (m_cuts, _q, p_cuts) in enumerate(mixed_cuts):
            letters = tuple(basis_of[settings[cut - 1]] for cut in m_cuts)
            expansions = [_WIRE_PREP[settings[cut - 1]] for cut in p_cuts]
            combined, original = {}, None
            for choice in itertools.product(*expansions):
                coefficient = 1.0
                base_combo = []
                for weight, base_state in choice:
                    coefficient *= weight
                    base_combo.append(base_state)
                corr, original = mixed_corr(frag_idx, letters, tuple(base_combo))
                for mask, value in corr.items():
                    combined[mask] = combined.get(mask, 0.0) + coefficient * value
            per_component.append((combined, original))
        per_component.extend(plain_flat)
        configs.append(per_component)

    return _reconstruct_wire_cut_distribution(
        configs, num_cuts, num_qubits, cut_source_index, prepared_to_orig
    )


def _wire_sampling_optimized(
    dag, wire_cuts, num_qubits, shots, software, qpu, qpu_dict, seed, joint=False
):
    """Fragment-separated sampling reconstruction.

    The sampling counterpart of :func:`_wire_cutting_optimized`: sample each measure
    fragment once per measurement basis, each prepare fragment once per base
    preparation, and each mixed fragment (owns both an ``O_i`` and a ``PS_i``
    marker, e.g. a chain of cuts in series) once per ``(basis_combo, prep_combo)``
    pair -- applying both marker substitutions to the same fragment copy before
    sampling, since :func:`_sample_gate_subcircuit` is agnostic to what the graph
    represents -- then reconstruct the distribution over the ``8**n``
    configurations from those cached correlators. Returns the reconstruction
    future, or ``None`` only when a cut's marker went missing entirely (a
    bookkeeping inconsistency).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _sample_gate_subcircuit

    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, mixed_frags, plain_frags = fragments

    cut_source_index, prepared_to_orig = _wire_cut_line_maps(
        _place_wire_cuts(dag.copy(), wire_cuts, num_qubits), num_qubits
    )

    counter = itertools.count()

    def _sample(component):
        derived = next_seed(seed, next(counter))
        return _sample_gate_subcircuit(
            component, software, shots, qpu, qpu_dict, seed=derived
        )

    measure_flat, measure_keys = [], []
    for frag_idx, (sub, cuts, cut_qubits) in enumerate(measure_frags):
        for setting in _measurement_settings(len(cuts), joint):
            graph = sub.copy()
            _apply_measurement_setting(graph, cuts, cut_qubits, setting)
            measure_flat.append(_sample(graph))
            measure_keys.append((frag_idx, setting))
    prepare_flat, prepare_keys = [], []
    for frag_idx, (sub, cuts) in enumerate(prepare_frags):
        for prep_combo in itertools.product(_PREP_STATES, repeat=len(cuts)):
            graph = sub.copy()
            for cut, state in zip(cuts, prep_combo):
                _apply_wire_prep_state(graph, cut, state)
            prepare_flat.append(_sample(graph))
            prepare_keys.append((frag_idx, prep_combo))
    mixed_flat, mixed_keys = [], []
    for frag_idx, (sub, m_cuts, cut_qubits, p_cuts) in enumerate(mixed_frags):
        for setting in _measurement_settings(len(m_cuts), joint):
            for prep_combo in itertools.product(_PREP_STATES, repeat=len(p_cuts)):
                graph = sub.copy()
                for cut, state in zip(p_cuts, prep_combo):
                    _apply_wire_prep_state(graph, cut, state)
                _apply_measurement_setting(graph, m_cuts, cut_qubits, setting)
                mixed_flat.append(_sample(graph))
                mixed_keys.append((frag_idx, setting, prep_combo))
    plain_flat = [_sample(sub.copy()) for sub in plain_frags]

    return _reconstruct_optimized_wire_distribution_task(
        measure_flat,
        measure_keys,
        [(cuts, q) for _s, cuts, q in measure_frags],
        prepare_flat,
        prepare_keys,
        [cuts for _s, cuts in prepare_frags],
        mixed_flat,
        mixed_keys,
        [(m_cuts, q, p_cuts) for _s, m_cuts, q, p_cuts in mixed_frags],
        plain_flat,
        len(wire_cuts),
        num_qubits,
        cut_source_index,
        prepared_to_orig,
    )


def _wire_cutting_sampling(
    circuit,
    cut,
    shots=1024,
    software="qiskit",
    sync=True,
    qpu=False,
    qpu_dict=None,
    seed=None,
    optimized: bool = False,
    joint: bool = True,
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
        distinct derived seed; ``None`` keeps the previous
        non-deterministic sampling. Ignored on a QPU.
    :param optimized: Sample each fragment once per measurement setting /
        preparation (``sum_f 3**i_f * 4**j_f``, or fewer with ``joint``) instead of
        every one of the ``8**n`` configurations.
    :param joint: With ``optimized=True``, measure the cut wires each fragment
        sends together: ``2**k + 1`` settings for ``k`` wires instead of ``3**k``
        (see :func:`wire_cutting`). ``False`` measures each wire on its own, with no
        entangling gates.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _sample_gate_subcircuit,
        _circuit_width,
    )

    dag = convert_to_dag(circuit)
    num_qubits = _circuit_width(dag, circuit)
    wire_cuts = _expand_wire_cuts(dag, cut)
    _validate_wire_cuts_partition(dag, wire_cuts, num_qubits)

    if not wire_cuts:
        from Qdislib.core.backend.registry import build_execution_backend

        backend = build_execution_backend(
            software, num_qubits=num_qubits, qpu=bool(qpu), qpu_dict=qpu_dict
        )
        counts = backend.evaluate_samples(dag, num_qubits, shots, seed=seed)
        total = sum(counts.values()) or 1
        return {k: v / total for k, v in counts.items()}

    if optimized:
        # Fragment-separated sampling: sum_f 3**i_f * 4**j_f fragment samplings
        # instead of 8**n (measure, prepare, and mixed -- chain-topology --
        # fragments all handled). Falls back to the standard method only on
        # error, so optimized=True never changes the (shot-noise) result.
        try:
            optimized_distribution = _wire_sampling_optimized(
                dag.copy(),
                wire_cuts,
                num_qubits,
                shots,
                software,
                qpu,
                qpu_dict,
                seed,
                joint=joint,
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
        _warn_not_optimized("sampling", reason)

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
            seed=next_seed(seed, i),
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

        subs = qd.wire_cutting_sampling_subcircuits(circuit, cut)
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(sampled)

    .. deprecated:: 2.0.0
        What a 1.x list was: no function returns it any more (their lists are
        :class:`~Qdislib.CutCircuits`), and it will be removed.
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


def _validate_wire_cuts_partition(dag, wire_cuts, num_qubits):
    """Raise unless every wire cut separates the circuit (see
    :func:`~Qdislib.core.cutting_algorithms._reconstruction.check_cuts_separate`).

    Wire cutting trades circuit count for reduced per-circuit qubit width; that
    trade only pays off if the cut(s) actually split the circuit into more
    independently-executable pieces than it already had. A cut whose wire isn't
    the sole connectivity between the two sides it's meant to separate (some
    other uncut path still links them) leaves the graph exactly as connected as
    before, so reconstructing through it pays wire cutting's shot-count cost for
    none of its benefit -- running the uncut circuit directly is strictly
    better, and cheaper, and less noisy.

    Compares the number of connected components before and after inserting the
    cut markers (on a throwaway copy), so this is correct even when the circuit
    already had multiple disconnected pieces before any cut -- it only cares
    whether THIS cut set added separation beyond that baseline. An empty
    ``wire_cuts`` is exempt: that means "no cut," which runs the circuit as-is
    and needs no partition.

    :param dag: The (uncut) circuit DAG.
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples, as
        returned by :func:`_expand_wire_cuts`.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :raises InvalidCutError: If a cut's measured and prepared sides stay in one
        fragment.
    """
    if not wire_cuts:
        return
    from Qdislib.core.cutting_algorithms._reconstruction import check_cuts_separate

    marked = dag.copy()
    _insert_wire_cut_markers(marked, wire_cuts, num_qubits)
    check_cuts_separate(
        _remove_red_edges(marked),
        {cut: (f"O_{i}", f"PS_{i}") for i, cut in enumerate(wire_cuts, start=1)},
    )


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
def _execute_subcircuit(prepared, software, shots, backend, qpu, seed=None):
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
    return exec_backend.evaluate_expectation(s_new, highest_qubit, (), shots, seed=seed)


def _build_and_run_subcircuit(component, software, shots, backend, qpu, seed=None):
    """Prepare then execute one wire-cut subcircuit as two COMPSs tasks.

    Orchestrates :func:`_create_subcircuit` (prepare) then :func:`_execute_subcircuit`
    (build + run), so subcircuit creation and execution are separate, independently
    routable tasks (the prepared-spec future feeds create -> execute directly,
    task -> task). Called from the master, so it submits both without nesting.
    Returns the execution future (the subcircuit's expectation value).
    """
    prepared = _create_subcircuit(component)
    return _execute_subcircuit(prepared, software, shots, backend, qpu, seed)


def _relabel_qubit_downstream(dag, target, qubit, fresh_qubit):
    """Relabel ``qubit`` onto ``fresh_qubit`` on ``target`` and every node
    downstream of it that still carries ``qubit`` forward.

    Used by wire-cut marker insertion to disconnect the "prepared" side of a
    cut from the "measured" side: everything from ``target`` onward that still
    references the cut qubit moves to a fresh line. This has to be transitive
    -- not just ``target``'s immediate successors -- since a node more than one
    hop downstream (e.g. an Observable-I terminal marker
    :func:`~Qdislib.utils.backends.change_basis_dag` inserted after further
    gates on that qubit) still carries the qubit forward and needs relabeling
    too, or it's silently left referencing the stale original index -- putting
    it in the wrong fragment, or at the wrong qubit position within the right
    one, once the graph is actually cut apart.

    :param dag: The circuit DAG (mutated in place).
    :param target: The node where the cut's downstream side starts.
    :param qubit: The original qubit index being relabeled.
    :param fresh_qubit: The fresh line index to relabel onto.
    """
    frontier = [target]
    visited = set()
    while frontier:
        node = frontier.pop()
        if node in visited:
            continue
        visited.add(node)
        node_qubits = dag.nodes[node].get("qubits")
        if not (node_qubits and qubit in node_qubits):
            continue
        dag.nodes[node]["qubits"] = tuple(
            fresh_qubit if q == qubit else q for q in node_qubits
        )
        frontier.extend(dag.successors(node))


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
    _place_wire_cuts(dag, wire_cuts, num_qubits)
    return dag


def _place_wire_cuts(dag, wire_cuts, num_qubits):
    """Do the cut surgery for ``wire_cuts``, in order; return the line each severed.

    For cut ``i`` = ``(source, target, qubit)``: remove the ``source -> target``
    edge, move ``target`` and everything downstream on that wire to the fresh line
    ``num_qubits + i``, and add ``O_i`` (measure, after ``source``) and ``PS_i``
    (prepare, before ``target``). Mutates ``dag``.

    ``qubit`` names the wire as in the uncut circuit. Once an earlier cut on the
    same wire has moved ``source`` onto a fresh line, the wire runs there instead,
    so each cut looks up its line at ``source`` by position (relabelling keeps a
    node's qubit order). Using the stale index, the second of two cuts in series
    severed nothing and the reconstruction was wrong.

    :return: ``[(source, target, line), ...]``, the line each cut actually severed.
    """
    original = {node: tuple(dag.nodes[node].get("qubits") or ()) for node in dag}
    placed = []
    for index, (source, target, qubit) in enumerate(wire_cuts, start=1):
        line = dag.nodes[source]["qubits"][original[source].index(qubit)]
        # The shared edge is removed once; further wires between the same pair
        # are already disconnected.
        if dag.has_edge(source, target):
            dag.remove_edge(source, target)

        # Move this wire on the target and everything downstream of it onto a
        # fresh line so the prepared (PS) side is disconnected from the measured
        # (O) side.
        _relabel_qubit_downstream(dag, target, line, num_qubits + index)

        dag.add_node(f"O_{index}", gate="s", qubits=(line,), parameters=())
        dag.add_node(
            f"PS_{index}", gate="t", qubits=(num_qubits + index,), parameters=()
        )
        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")
        placed.append((source, target, line))
    return placed


def _wire_cut_line_maps(placed, num_qubits):
    """Sampling bookkeeping for placed cuts (see :func:`_place_wire_cuts`).

    :return: ``(cut_source_index, prepared_to_orig)``: the line each cut measures
        -> its 1-based cut index, and each cut's fresh line -> the circuit qubit it
        carries on (followed through earlier cuts on the same wire).
    """
    cut_source_index, prepared_to_orig = {}, {}
    for index, (_, _, line) in enumerate(placed, start=1):
        cut_source_index[line] = index
        prepared_to_orig[num_qubits + index] = prepared_to_orig.get(line, line)
    return cut_source_index, prepared_to_orig


def _generate_wire_cutting(
    dag: typing.Any,
    wire_cuts: typing.List[typing.Tuple[str, str, int]],
    num_qubits: int,
    shots: int = 10000,
    backend: str = "numpy",
    qpu=None,
    software: str = "qibo",
    seed=None,
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
        copy_graph = _remove_red_edges(dag)

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
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
            _build_and_run_subcircuit(
                s, software, shots, backend, qpu, next_seed(seed, index, i)
            )
            for i, s in enumerate(graph_components)
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
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "x"

        # X +
        elif index == 2:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "h"
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "h"

        # X -
        elif index == 3:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "h"
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "h"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="x",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")

        # Y +i
        elif index == 4:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "sdg"
            updated_dag.add_node(
                f"O2_{idx2}",
                gate="h",
                qubits=updated_dag.nodes[f"O_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"O_{idx2}", f"O2_{idx2}", color="blue")
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "s"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="h",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")

        # Y -i
        elif index == 5:
            updated_dag.nodes[f"O_{idx2}"]["gate"] = "sdg"
            updated_dag.add_node(
                f"O2_{idx2}",
                gate="h",
                qubits=updated_dag.nodes[f"O_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"O_{idx2}", f"O2_{idx2}", color="blue")
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "s"
            updated_dag.add_node(
                f"PS2_{idx2}",
                gate="h",
                qubits=updated_dag.nodes[f"PS_{idx2}"].get("qubits"),
                parameters=(),
            )
            updated_dag.add_edge(f"PS2_{idx2}", f"PS_{idx2}", color="blue")
            updated_dag.add_node(
                f"PS3_{idx2}",
                gate="x",
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
            updated_dag.nodes[f"PS_{idx2}"]["gate"] = "x"

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
    for i, c in enumerate(
        networkx.connected_components(updated_dag.to_undirected(as_view=True))
    ):
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

#: The Pauli each base-8 cut setting needs from the measure side (0/1: identity).
_SETTING_PAULI = {0: "I", 1: "I", 2: "X", 3: "X", 4: "Y", 5: "Y", 6: "Z", 7: "Z"}

# ---------------------------------------------------------------------------
# Joint measurement of the cut wires a fragment sends. The measure side must supply
# <O (x) P> for every Pauli P on those k wires. Per wire, that takes 3**k circuits
# (one basis per wire; a qubit's I is read by ignoring its bit). Measured jointly,
# the 4**k - 1 non-identity Paulis split into 2**k + 1 commuting groups of 2**k - 1
# (mutually unbiased bases): one circuit per group reads all of it, every Pauli
# exactly once. The groups come from the field GF(2**k): group t holds X^a Z^(A_t a)
# with (A_t)_ij = Tr(t e_i e_j), plus the all-Z group. A_t is symmetric, so each
# group is measured by SDG (where A_t has a 1 on the diagonal), CZ (off-diagonal 1s)
# and H on every wire.
# ---------------------------------------------------------------------------

#: Irreducible polynomials over GF(2) (bit i = coefficient of x**i), per degree k.
_IRREDUCIBLE = {
    1: 0b11,
    2: 0b111,
    3: 0b1011,
    4: 0b10011,
    5: 0b100101,
    6: 0b1000011,
    7: 0b10000011,
    8: 0b100011011,
}


def _gf_mul(a, b, k):
    """Product in GF(2**k), elements as bit masks over the polynomial basis."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if (a >> k) & 1:
            a ^= _IRREDUCIBLE[k]
    return out


def _gf_trace(y, k):
    """Field trace GF(2**k) -> GF(2): ``y + y**2 + ... + y**(2**(k-1))``."""
    total, power = 0, y
    for _ in range(k):
        total ^= power
        power = _gf_mul(power, power, k)
    return total


def _pauli_matrix(label):
    matrices = {
        "I": np.eye(2),
        "X": np.array([[0, 1], [1, 0]]),
        "Y": np.array([[0, -1j], [1j, 0]]),
        "Z": np.diag([1, -1]),
    }
    out = np.array([[1.0 + 0j]])
    for letter in label:
        out = np.kron(out, matrices[letter])
    return out


def _gate_matrix(name, qubits, k):
    """``name`` on local ``qubits`` of a k-qubit register (qubit 0 most significant)."""
    single = {
        "H": np.array([[1, 1], [1, -1]]) / np.sqrt(2),
        "SDG": np.diag([1, -1j]),
    }
    name = name.upper()
    if name in single:
        out = np.array([[1.0 + 0j]])
        for q in range(k):
            out = np.kron(out, single[name] if q == qubits[0] else np.eye(2))
        return out
    i, j = qubits  # CZ
    return np.diag(
        [
            -1 if ((x >> (k - 1 - i)) & 1 and (x >> (k - 1 - j)) & 1) else 1
            for x in range(2**k)
        ]
    ).astype(complex)


@functools.lru_cache(maxsize=None)
def _joint_measurement_settings(k):
    """The ``2**k + 1`` joint measurement settings for ``k`` cut wires.

    :return: A tuple of ``(gates, readout)``: ``gates`` is ``[(name, local_qubits),
        ...]`` over the cut wires ``0..k-1``; ``readout`` is ``[(pauli, sign,
        bits), ...]`` with ``<O (x) pauli> = sign * <O (x) Z_bits>`` measured after
        the gates, for every Pauli the setting gives (the all-identity one too).
    """
    if k not in _IRREDUCIBLE:
        raise ValueError(
            f"joint measurement supports 1 to 8 wires per fragment, not {k}"
        )
    groups = []
    for t in range(2**k):
        a = [
            [_gf_trace(_gf_mul(t, _gf_mul(1 << i, 1 << j, k), k), k) for j in range(k)]
            for i in range(k)
        ]
        x_parts = [
            [int(i == j) for i in range(k)] for j in range(k)
        ]  # generator j: X_j
        z_parts = [[a[i][j] for i in range(k)] for j in range(k)]  # times Z^(column j)
        gates = [("sdg", (i,)) for i in range(k) if a[i][i]]
        gates += [("cz", (i, j)) for i in range(k) for j in range(i + 1, k) if a[i][j]]
        gates += [("h", (i,)) for i in range(k)]
        groups.append((gates, x_parts, z_parts))
    identity = [[int(i == j) for i in range(k)] for j in range(k)]
    groups.append(([], [[0] * k for _ in range(k)], identity))  # the all-Z group

    settings, seen = [], set()
    for gates, x_parts, z_parts in groups:
        circuit = np.eye(2**k, dtype=complex)
        for name, qubits in gates:
            circuit = _gate_matrix(name, qubits, k) @ circuit
        readout = [(("I",) * k, 1, ())]
        for coeffs in itertools.product((0, 1), repeat=k):
            if not any(coeffs):
                continue
            x = [sum(c * g[q] for c, g in zip(coeffs, x_parts)) % 2 for q in range(k)]
            z = [sum(c * g[q] for c, g in zip(coeffs, z_parts)) % 2 for q in range(k)]
            pauli = tuple("IXZY"[x[q] + 2 * z[q]] for q in range(k))
            image = circuit @ _pauli_matrix(pauli) @ circuit.conj().T
            sign = int(round(image[0, 0].real))
            bits = tuple(
                q
                for q in range(k)
                if round(image[1 << (k - 1 - q)][1 << (k - 1 - q)].real) == -sign
            )
            expected = np.diag(
                [
                    sign * (-1) ** sum((x_ >> (k - 1 - q)) & 1 for q in bits)
                    for x_ in range(2**k)
                ]
            )
            if not np.allclose(image, expected):  # pragma: no cover - construction bug
                raise RuntimeError(f"joint setting does not diagonalise {pauli}")
            readout.append((pauli, sign, bits))
            seen.add(pauli)
        settings.append((tuple(gates), tuple(readout)))
    if len(seen) != 4**k - 1:  # pragma: no cover - construction bug
        raise RuntimeError("joint settings do not cover every Pauli once")
    return tuple(settings)


def _per_wire_readout(basis_combo):
    """Readout of a per-wire setting: every Pauli with letters ``b_i`` or ``I``
    (``I`` only on Z-basis wires, so each Pauli is read by exactly one setting)."""
    options = [(["P", "I"] if basis == "Z" else ["P"]) for basis in basis_combo]
    readout = []
    for extraction in itertools.product(*options):
        pauli = tuple(b if e == "P" else "I" for b, e in zip(basis_combo, extraction))
        bits = tuple(i for i, e in enumerate(extraction) if e == "P")
        readout.append((pauli, 1, bits))
    return readout


def _measurement_settings(num_cuts, joint=False):
    """The settings that measure ``num_cuts`` cut wires of one fragment: per-wire
    basis combos (``3**k``), or joint ones (``2**k + 1``) when ``joint``."""
    if joint and num_cuts > 1:
        return [("joint", t) for t in range(2**num_cuts + 1)]
    return list(itertools.product(_MEASURE_BASES, repeat=num_cuts))


def _apply_measurement_setting(graph, cuts, cut_qubits, setting):
    """Turn the fragment's ``O_i`` markers into one measurement setting, in place.

    A per-wire setting rotates each cut qubit into its own basis
    (:func:`_apply_wire_measure_option`). A joint setting removes the markers and
    appends its gates on the cut qubits after their last gate -- nothing else acts
    on a cut qubit in this fragment after the cut, so measuring later is the same.

    :return: The setting's readout (see :func:`_joint_measurement_settings`).
    """
    if setting and setting[0] == "joint":
        gates, readout = _joint_measurement_settings(len(cuts))[setting[1]]
        last = {}
        for local, cut in enumerate(cuts):
            node = f"O_{cut}"
            preds = list(graph.predecessors(node))
            last[local] = preds[0] if preds else None
            graph.remove_node(node)
        for index, (name, qubits) in enumerate(gates):
            node = f"JM_{cuts[0]}_{index}"
            graph.add_node(
                node,
                gate=name,
                qubits=tuple(cut_qubits[cuts[q]] for q in qubits),
                parameters=(),
            )
            for q in qubits:
                if last[q] is not None:
                    graph.add_edge(last[q], node, color="blue")
                last[q] = node
        return readout
    for cut, basis in zip(cuts, setting):
        _apply_wire_measure_option(graph, cut, basis + "P")
    return _per_wire_readout(setting)


def _measure_table(counts, cut_positions, readout):
    """``{pauli: <O (x) pauli>}`` from one setting's counts: ``O`` is the parity of
    the fragment's other qubits, each Pauli the setting's ``sign`` times the parity
    of its ``bits`` among the cut qubits."""
    total = sum(counts.values()) or 1
    width = len(next(iter(counts), "")) if counts else 0
    observable_positions = [p for p in range(width) if p not in cut_positions]
    table = {}
    for pauli, sign, bits in readout:
        selected = [cut_positions[b] for b in bits]
        value = 0.0
        for bitstring, count in counts.items():
            parity = 0
            for pos in observable_positions:
                parity ^= int(bitstring[pos])
            for pos in selected:
                parity ^= int(bitstring[pos])
            value += -count if parity else count
        table[pauli] = sign * value / total
    return table


def _merge_tables(tables):
    """One ``{pauli: value}`` table from several settings' tables; a Pauli read by
    more than one setting (the all-identity one, when joint) is averaged."""
    sums, counts = {}, {}
    for table in tables:
        for pauli, value in table.items():
            sums[pauli] = sums.get(pauli, 0.0) + value
            counts[pauli] = counts.get(pauli, 0) + 1
    return {pauli: sums[pauli] / counts[pauli] for pauli in sums}


def _setting_readout(setting, num_cuts):
    """The readout of a measurement setting (per-wire basis combo or joint)."""
    if setting and setting[0] == "joint":
        return _joint_measurement_settings(num_cuts)[setting[1]][1]
    return _per_wire_readout(setting)


def _corr_pauli_table(corr, original, cuts, cut_qubits, readout):
    """``{(output_mask, pauli): value}`` from one setting's Z correlators.

    ``corr`` maps a mask over ``original`` (bit i = Z on ``original[i]``) to its
    correlator. After the setting's gates, a mask whose cut bits equal a readout
    entry's ``bits`` gives ``sign * <Z^output (x) pauli>`` on the remaining bits.
    """
    positions = [original.index(cut_qubits[cut]) for cut in cuts]
    cut_mask = sum(1 << p for p in positions)
    table = {}
    for pauli, sign, bits in readout:
        bits_mask = sum(1 << positions[b] for b in bits)
        for mask, value in corr.items():
            if mask & cut_mask == bits_mask:
                table[(mask & ~cut_mask, pauli)] = sign * value
    return table


def _corr_from_pauli_table(table, original, cuts, cut_qubits, letters):
    """The per-wire-shaped correlators the distribution reconstruction reads: a
    cut bit set means that wire's Pauli ``letters[i]``, unset means identity."""
    positions = [original.index(cut_qubits[cut]) for cut in cuts]
    cut_mask = sum(1 << p for p in positions)
    corr = {}
    for mask in range(2 ** len(original)):
        pauli = tuple(
            letters[i] if (mask >> p) & 1 else "I" for i, p in enumerate(positions)
        )
        corr[mask] = table.get((mask & ~cut_mask, pauli), 0.0)
    return corr


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
        graph.nodes[node]["gate"] = "x"
    elif state == "+":
        graph.nodes[node]["gate"] = "h"
    elif state == "+i":
        graph.nodes[node]["gate"] = "s"
        rotate = f"{node}_Hi"
        graph.add_node(
            rotate, gate="h", qubits=graph.nodes[node]["qubits"], parameters=()
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
        graph.nodes[node]["gate"] = "h"
    elif option == "YP":
        graph.nodes[node]["gate"] = "sdg"
        rotate = f"{node}_HY"
        graph.add_node(
            rotate, gate="h", qubits=graph.nodes[node]["qubits"], parameters=()
        )
        graph.add_edge(node, rotate, color="blue")
    elif option == "ZP":
        graph.remove_node(node)


def _sample_wire_fragment(graph, shots, backend, software, seed=None):
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
    return exec_backend.evaluate_samples(graph, highest, shots, seed=seed)


@task(returns=1)
def _run_wire_measure_fragment(
    sub, cuts, cut_qubits, basis_combo, shots, backend, software, seed=None
):
    """Evaluate one measure ("top") fragment for a measurement-basis configuration.

    Rotates each cut qubit's ``O_i`` marker into the requested basis (X -> H,
    Y -> SDG + H, Z -> measure directly), samples the fragment, then extracts, for
    every ``P``/``I`` combination of its cut qubits, the joint expectation of those
    Pauli/identity factors with ``Z`` on the fragment's other (observable) qubits.

    :param sub: The measure fragment DAG (carries the ``O_i`` markers).
    :param cuts: The cut indices this fragment owns (sorted).
    :param cut_qubits: ``{cut_index: original qubit}`` for each owned cut.
    :param basis_combo: The measurement setting: a per-wire ``("X"|"Y"|"Z", ...)``
        aligned with ``cuts``, or ``("joint", t)`` (see :func:`_measurement_settings`).
    :param shots: Number of measurement shots.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :return: ``{pauli: expectation}``, ``pauli`` a tuple of ``"I"/"X"/"Y"/"Z"`` over
        ``cuts``, for every Pauli the setting reads.
    """
    graph = sub.copy()
    used = sorted(
        {q for node in graph.nodes for q in graph.nodes[node].get("qubits", ())}
    )
    new_position = {q: i for i, q in enumerate(used)}
    cut_positions = [new_position[cut_qubits[cut]] for cut in cuts]

    readout = _apply_measurement_setting(graph, cuts, cut_qubits, basis_combo)
    counts = _sample_wire_fragment(graph, shots, backend, software, seed)
    return _measure_table(counts, cut_positions, readout)


class _OptimizedWireFragment(Subcircuit):
    """A fragment circuit of the optimized three-step wire cut: its markers set to
    one measurement setting (the wires it sends) and one prepared state per wire it
    receives. Running it returns the ``{pauli: value}`` table of every Pauli that
    setting reads, so ``reconstruct`` needs fewer circuits than ``8**n``."""

    def __init__(
        self, sub, measure_cuts, cut_qubits, setting, prepare_cuts, prep, software
    ):
        qubits = {q for _, d in sub.nodes(data=True) for q in d.get("qubits") or ()}
        super().__init__(sub, len(qubits), software)
        self.measure_cuts = measure_cuts
        self.cut_qubits = cut_qubits
        self.basis_combo = setting
        self.prepare_cuts = prepare_cuts
        self.prep_combo = prep

    @property
    def sub(self):
        return self.dag

    @property
    def cuts(self):
        return self.measure_cuts

    def to_circuit(self):
        from Qdislib.utils.backends import _materialize

        graph = self.dag.copy()
        for cut, state in zip(self.prepare_cuts, self.prep_combo):
            _apply_wire_prep_state(graph, cut, state)
        _apply_measurement_setting(
            graph, self.measure_cuts, self.cut_qubits, self.basis_combo
        )
        graph, width = _update_qubits(graph)
        return _materialize(graph, width, self.software)

    @property
    def observable(self):
        return None

    def _expectation(self, shots, backend, seed, cache, cache_mode):
        software = self.software if backend is None else backend
        if not self.prepare_cuts:
            return _run_wire_measure_fragment(
                self.dag,
                self.measure_cuts,
                self.cut_qubits,
                self.basis_combo,
                shots,
                None,
                software,
                seed,
            )
        return _run_wire_mixed_fragment(
            self.dag,
            self.measure_cuts,
            self.cut_qubits,
            self.basis_combo,
            self.prepare_cuts,
            self.prep_combo,
            shots,
            None,
            software,
            seed,
        )

    def evaluate(self, shots):  # the 1.x method
        return self._expectation(shots, None, None, None, None)


def _OptimizedWireMeasure(sub, cuts, cut_qubits, basis_combo, software):
    """A measure fragment (sends wires only) of the optimized three-step cut."""
    return _OptimizedWireFragment(sub, cuts, cut_qubits, basis_combo, [], (), software)


def _OptimizedWireMixed(
    sub, measure_cuts, cut_qubits, basis_combo, prepare_cuts, prep_combo, software
):
    """A fragment that sends and receives wires, for one prepared-state combo."""
    return _OptimizedWireFragment(
        sub, measure_cuts, cut_qubits, basis_combo, prepare_cuts, prep_combo, software
    )


@task(returns=1)
def _run_wire_prepare_fragment(
    sub, cuts, prep_combo, shots, backend, software, seed=None
):
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

    counts = _sample_wire_fragment(graph, shots, backend, software, seed)
    total = sum(counts.values()) or 1
    value = 0.0
    for bitstring, count in counts.items():
        parity = 1 if bitstring.count("1") % 2 == 0 else -1
        value += parity * (count / total)
    return value


@task(returns=1)
def _run_wire_mixed_fragment(
    sub,
    measure_cuts,
    cut_qubits,
    basis_combo,
    prepare_cuts,
    prep_combo,
    shots,
    backend,
    software,
    seed=None,
):
    """Evaluate one *mixed* fragment -- owns both ``O_i`` and ``PS_i`` markers.

    A fragment ends up mixed when a cut's downstream continuation feeds another
    cut's upstream side without any other separation in between (e.g. two cuts in
    series on one wire, or more generally any chain topology): the middle
    fragment(s) receive one wire via a ``PS_i`` marker and send another via an
    ``O_i`` marker, so they can't be split into an independent measure fragment
    and prepare fragment the way :func:`_run_wire_measure_fragment` /
    :func:`_run_wire_prepare_fragment` assume.

    Generalizes both: apply the ``prepare_cuts`` substitutions first (exactly as
    :func:`_run_wire_prepare_fragment` does), then rotate the ``measure_cuts``
    markers into ``basis_combo`` and extract per-cut P/I correlators exactly as
    :func:`_run_wire_measure_fragment` does -- the prepared qubits simply become
    part of this fragment's own observable positions, no different from a
    fragment that was never touched by a prepare marker. One circuit run covers
    one ``(basis_combo, prep_combo)`` pair, so a fragment owning ``i`` measure-cuts
    and ``j`` prepare-cuts costs ``3**i * 4**j`` runs -- which is exactly
    :func:`_run_wire_measure_fragment`'s ``3**i`` when ``j=0`` and
    :func:`_run_wire_prepare_fragment`'s ``4**j`` when ``i=0``; this function
    generalizes both rather than replacing them, since the pure cases stay
    cheaper run standalone (a mixed fragment's cost is always >= splitting it
    would have cost, when splitting is possible at all).

    :param sub: The fragment DAG (carries both ``O_i`` and ``PS_i`` markers).
    :param measure_cuts: The measure-type cut indices this fragment owns (sorted).
    :param cut_qubits: ``{cut_index: original qubit}`` for each owned measure-type cut.
    :param basis_combo: The measurement setting for ``measure_cuts`` (per-wire or
        joint, see :func:`_measurement_settings`).
    :param prepare_cuts: The prepare-type cut indices this fragment owns (sorted).
    :param prep_combo: One of ``("0"|"1"|"+"|"+i", ...)`` aligned with ``prepare_cuts``.
    :param shots: Number of measurement shots.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :return: ``{pauli: expectation}``, same shape as
        :func:`_run_wire_measure_fragment` (``pauli`` is ``()`` when
        ``measure_cuts`` is empty, matching :func:`_run_wire_prepare_fragment`'s
        single ``<Z...Z>`` parity).
    """
    graph = sub.copy()
    for cut, prep in zip(prepare_cuts, prep_combo):
        _apply_wire_prep_state(graph, cut, prep)

    used = sorted(
        {q for node in graph.nodes for q in graph.nodes[node].get("qubits", ())}
    )
    new_position = {q: i for i, q in enumerate(used)}
    cut_positions = [new_position[cut_qubits[cut]] for cut in measure_cuts]

    readout = _apply_measurement_setting(graph, measure_cuts, cut_qubits, basis_combo)
    counts = _sample_wire_fragment(graph, shots, backend, software, seed)
    return _measure_table(counts, cut_positions, readout)


def _optimized_wire_fragments(dag, wire_cuts, num_qubits):
    """Insert cut markers and split the circuit into measure / prepare / mixed / plain fragments.

    Shared by the all-in-one and three-step optimized paths.

    A fragment carrying BOTH a measure and a prepare marker (e.g. two cuts in
    series on one wire, or any chain topology where a cut's downstream
    continuation feeds another cut's upstream side) used to make the whole
    circuit "not fragment-separable" and fall back to the standard ``8**n``
    method. It no longer does: such a fragment is classified as *mixed* and run
    via :func:`_run_wire_mixed_fragment` (``3**i * 4**j`` for the ``i`` measure- +
    ``j`` prepare-type cuts it owns) instead of bailing out for the entire
    circuit. Every fragment -- measure, prepare, mixed, or plain -- is still
    independently, minimally re-run, so a chain of ``n`` cuts costs
    ``sum_f 3**i_f * 4**j_f`` over its fragments rather than ``8**n`` globally;
    that sum is always strictly less than ``8**n`` for ``n >= 1``, though a mixed
    fragment's own ``3**i * 4**j`` is more expensive than splitting it into a
    pure measure fragment (``3**i``) and pure prepare fragment (``4**j``) would
    have been, on the (usual) rest of the circuit where that split is possible.

    :param dag: The circuit DAG (mutated: cut markers inserted).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :return: ``(measure_frags, prepare_frags, mixed_frags, plain_frags)`` where
        ``measure_frags`` is ``[(sub, owned_cuts, {cut: qubit}), ...]``,
        ``prepare_frags`` is ``[(sub, owned_cuts), ...]``, ``mixed_frags`` is
        ``[(sub, measure_cuts, {cut: qubit}, prepare_cuts), ...]`` and
        ``plain_frags`` is ``[sub, ...]``; or ``None`` only when a cut's marker is
        missing entirely (a bookkeeping inconsistency, not a topology this
        function can't handle).
    """
    n_cuts = len(wire_cuts)
    _insert_wire_cut_markers(dag, wire_cuts, num_qubits)
    graph = _remove_red_edges(dag)

    measure_frags, prepare_frags, mixed_frags, plain_frags = [], [], [], []
    measured, prepared = set(), set()
    for component in networkx.connected_components(graph.to_undirected(as_view=True)):
        sub = graph.subgraph(component).copy()
        measures = sorted(i for i in range(1, n_cuts + 1) if f"O_{i}" in sub.nodes)
        preps = sorted(i for i in range(1, n_cuts + 1) if f"PS_{i}" in sub.nodes)
        if measures and preps:
            cut_qubits = {i: sub.nodes[f"O_{i}"]["qubits"][0] for i in measures}
            mixed_frags.append((sub, measures, cut_qubits, preps))
            measured.update(measures)
            prepared.update(preps)
        elif measures:
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
    return measure_frags, prepare_frags, mixed_frags, plain_frags


@task(
    returns=1,
    measure_flat=COLLECTION_IN,
    prepare_flat=COLLECTION_IN,
    mixed_flat=COLLECTION_IN,
    plain_flat=COLLECTION_IN,
)
def _reconstruct_optimized_wire_cut(
    measure_flat,
    measure_keys,
    measure_cuts,
    prepare_flat,
    prepare_keys,
    prepare_cuts,
    mixed_flat,
    mixed_keys,
    mixed_cuts,
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
    :param mixed_flat: Mixed-fragment results, one ``{extraction: value}`` dict per
        (fragment, basis_combo, prep_combo), in ``mixed_keys`` order.
    :param mixed_keys: ``(fragment_index, basis_combo, prep_combo)`` for each
        ``mixed_flat`` entry.
    :param mixed_cuts: ``(measure_cuts, prepare_cuts)`` owned by each mixed fragment.
    :param plain_flat: ``<Z...Z>`` of each cut-free fragment (a constant factor).
    :param n_cuts: Number of wire cuts.
    :return: The reconstructed expectation value.
    """
    # One {pauli: value} table per measure fragment (and per prepared input of a
    # mixed one), merged over its measurement settings.
    measure_tables = {}
    for (frag_idx, _setting), value in zip(measure_keys, measure_flat):
        measure_tables.setdefault(frag_idx, []).append(value)
    measure_values = {f: _merge_tables(t) for f, t in measure_tables.items()}
    prepare_values = {}
    for (frag_idx, prep_combo), value in zip(prepare_keys, prepare_flat):
        prepare_values.setdefault(frag_idx, {})[prep_combo] = value
    mixed_tables = {}
    for (frag_idx, _setting, prep_combo), value in zip(mixed_keys, mixed_flat):
        mixed_tables.setdefault(frag_idx, {}).setdefault(prep_combo, []).append(value)
    mixed_values = {
        f: {prep: _merge_tables(t) for prep, t in by_prep.items()}
        for f, by_prep in mixed_tables.items()
    }
    plain_factor = 1.0
    for value in plain_flat:
        plain_factor *= value

    total = 0.0
    for index in range(8**n_cuts):
        settings = [(index // (8 ** (n_cuts - 1 - k))) % 8 for k in range(n_cuts)]
        sign = -1.0 if sum(s in _WIRE_SIGN_FLIP for s in settings) % 2 else 1.0
        term = sign

        for frag_idx, cuts in enumerate(measure_cuts):
            pauli = tuple(_SETTING_PAULI[settings[cut - 1]] for cut in cuts)
            term *= measure_values[frag_idx][pauli]

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

        # Mixed fragments own cuts of BOTH types at once: the basis_combo/extraction
        # lookup (like a measure fragment) and the |-> / |-i> linear expansion (like
        # a prepare fragment) both apply to the SAME fragment's joint table -- one
        # physical circuit was run per (basis_combo, prep_base_combo) pair.
        for frag_idx, (m_cuts, p_cuts) in enumerate(mixed_cuts):
            pauli = tuple(_SETTING_PAULI[settings[cut - 1]] for cut in m_cuts)
            expansions = [_WIRE_PREP[settings[cut - 1]] for cut in p_cuts]
            fragment_value = 0.0
            for choice in itertools.product(*expansions):
                coefficient = 1.0
                base_combo = []
                for weight, base_state in choice:
                    coefficient *= weight
                    base_combo.append(base_state)
                fragment_value += (
                    coefficient * mixed_values[frag_idx][tuple(base_combo)][pauli]
                )
            term *= fragment_value

        total += term * plain_factor

    return total / (2**n_cuts)


def _wire_cutting_optimized(
    dag, wire_cuts, num_qubits, shots, backend, software, joint=False, seed=None
):
    """Fragment-separated wire-cut reconstruction.

    Splits the cut circuit into measure fragments (own cuts via their ``O_i``
    markers), prepare fragments (via ``PS_i``), and mixed fragments (own both --
    e.g. a chain of cuts in series). Each measure fragment runs its
    ``3**(#cuts)`` basis configurations, each prepare fragment its
    ``4**(#cuts)`` preparations, and each mixed fragment the joint
    ``3**(#measure-cuts) * 4**(#prepare-cuts)`` -- always cheaper than the
    standard method's global ``8**n``, though a mixed fragment's own product is
    costlier than splitting those same cuts into pure fragments would have been,
    were that split possible. The ``8**n`` signed sum (base-8, sign flip on
    ``{3, 5, 7}``, prefactor ``1/2**n``) is then evaluated classically over the
    cached fragment values, with |-> / |-i> expanded as linear combinations.
    Backend-agnostic; each fragment run is a PyCOMPSs task.

    :param dag: The circuit DAG (mutated: cut markers inserted).
    :param wire_cuts: Single-wire cuts as ``(source, target, qubit)`` tuples.
    :param num_qubits: Highest qubit index; fresh lines are numbered above it.
    :param shots: Number of measurement shots per fragment.
    :param backend: Qibo platform name.
    :param software: Execution backend selector.
    :param joint: Measure each fragment's sent wires jointly (``2**k + 1`` circuits
        instead of ``3**k``, see :func:`_joint_measurement_settings`).
    :return: The reconstructed expectation value, or ``None`` only when a cut's
        marker went missing entirely (a bookkeeping inconsistency).
    """
    n_cuts = len(wire_cuts)
    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, mixed_frags, plain_frags = fragments

    # Fan out every fragment execution (PyCOMPSs tasks), flattened with their keys.
    measure_flat, measure_keys = [], []
    for frag_idx, (sub, cuts, cut_qubits) in enumerate(measure_frags):
        for combo in _measurement_settings(len(cuts), joint):
            measure_flat.append(
                _run_wire_measure_fragment(
                    sub,
                    cuts,
                    cut_qubits,
                    combo,
                    shots,
                    backend,
                    software,
                    next_seed(seed, "m", len(measure_flat)),
                )
            )
            measure_keys.append((frag_idx, combo))
    prepare_flat, prepare_keys = [], []
    for frag_idx, (sub, cuts) in enumerate(prepare_frags):
        for combo in itertools.product(_PREP_STATES, repeat=len(cuts)):
            prepare_flat.append(
                _run_wire_prepare_fragment(
                    sub,
                    cuts,
                    combo,
                    shots,
                    backend,
                    software,
                    next_seed(seed, "p", len(prepare_flat)),
                )
            )
            prepare_keys.append((frag_idx, combo))
    mixed_flat, mixed_keys = [], []
    for frag_idx, (sub, m_cuts, cut_qubits, p_cuts) in enumerate(mixed_frags):
        for basis_combo in _measurement_settings(len(m_cuts), joint):
            for prep_combo in itertools.product(_PREP_STATES, repeat=len(p_cuts)):
                mixed_flat.append(
                    _run_wire_mixed_fragment(
                        sub,
                        m_cuts,
                        cut_qubits,
                        basis_combo,
                        p_cuts,
                        prep_combo,
                        shots,
                        backend,
                        software,
                        next_seed(seed, "x", len(mixed_flat)),
                    )
                )
                mixed_keys.append((frag_idx, basis_combo, prep_combo))
    plain_flat = [
        _run_wire_prepare_fragment(
            sub, [], (), shots, backend, software, next_seed(seed, "plain", n)
        )
        for n, sub in enumerate(plain_frags)
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
        mixed_flat,
        mixed_keys,
        [(m_cuts, p_cuts) for _sub, m_cuts, _q, p_cuts in mixed_frags],
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

    Appends a measurement on every qubit to a copy, so the caller's circuit is left
    untouched and can be run again. One with a collapse measurement is the exception:
    a copy has new measurement gates, and the circuit's symbols would stay bound to
    the old ones, so it gets the measurement in place and is single-use.

    :param qc: A Qibo ``Circuit``.
    :param shots: Number of measurement shots.
    :param backend: Qibo backend name.
    :return: A ``{bitstring: count}`` dict (qubit-0-first keys).
    """
    qibo.set_backend(backend)
    if not any(isinstance(g, gates.M) and g.collapse for g in qc.queue):
        qc = qc.copy(deep=True)
    qc.add(gates.M(*range(qc.nqubits)))
    result = qc(nshots=shots)
    return dict(result.frequencies(binary=True))


def _expectation_from_counts_qubit0_first(freq, observables, shots=1024):
    """Compute ``<O>`` for a Z/I Pauli string from measurement counts whose keys are
    qubit-0-first bit strings (qibo, CUDA-Q, PennyLane), which align directly with
    the observable string.

    :param freq: A ``{bitstring: count}`` dict, qubit-0-first (:func:`_counts_qibo`,
        :func:`_counts_cudaq`, :func:`_counts_pennylane`).
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


_expectation_from_counts_qibo = _expectation_from_counts_qubit0_first
_expectation_from_counts_cudaq = _expectation_from_counts_qubit0_first
_expectation_from_counts_pennylane = _expectation_from_counts_qubit0_first


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

    if not isinstance(subcirc, qibo.models.Circuit):
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


def _counts_cudaq(kernel, shots=1024, target=None):
    """Sample a CUDA-Q kernel and return its measurement counts.

    On ``target``; by default the GPU (``nvidia``) target, or the CPU (``qpp-cpu``)
    with a warning.

    :param kernel: A built CUDA-Q kernel.
    :param shots: Number of measurement shots.
    :param target: The CUDA-Q target, as a :class:`CudaqBackend` is configured.
    :return: A ``{bitstring: count}`` dict.
    """
    from Qdislib.core.backend.cudaq_backend import select_cudaq_target
    from Qdislib.utils.graph_cudaq import cudaq_without_measurements

    select_cudaq_target(target)
    freq = cudaq.sample(cudaq_without_measurements(kernel), shots_count=shots)
    return {str(key): int(value) for key, value in freq.items()}


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
    positions = {}  # every cut-point outcome of each qubit (a qubit can be cut twice)
    for i, qb in enumerate(measured_order):
        positions.setdefault(qb, []).append(i)
    total = 0.0
    drawn = 0
    for record in records:
        drawn += 1
        contribution = 1
        # record is an int; bit k = outcome k of [cut-points..., finals...].
        for q, obs in enumerate(observables):
            if q in prep_qubits and q in positions:
                # Measured cut qubit: decomposition sign -(-1)**m from each
                # cut-point outcome m for *any* observable; Z then multiplies by
                # the final measurement.
                for pos in positions[q]:
                    contribution *= -((-1) ** ((record >> pos) & 1))
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
def _expec_value_cudaq(
    subcirc: typing.Any, shots=1024, prep_qubits=None, target=None
) -> float:
    """Execute the given circuit.

    :param subcirc: Circuit to execute.
    :param prep_qubits: Renumbered cut qubits measured for a gate cut; when any is
        observed with identity, its decomposition sign is recovered from the
        cut-point measurement (see :func:`_expectation_cudaq_decoupled`).
    :param target: The CUDA-Q target (see :func:`_counts_cudaq`).
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
        from Qdislib.core.backend.cudaq_backend import select_cudaq_target

        select_cudaq_target(target)
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

    freq = _counts_cudaq(subcirc, shots, target)
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


def _fold_cut_outcomes(outcomes):
    """One bit standing for all of a qubit's cut-point outcomes.

    Each cut measurement on a qubit contributes ``-(-1)**m`` to the
    decomposition sign, and they multiply: ``(-1)**n * (-1)**(m_1+...+m_n)``.
    The readers apply ``-(-1)**bit`` once per qubit, which equals that product
    for ``bit = (m_1 + ... + m_n + n - 1) % 2``. Works on ints or per-shot arrays.

    :param outcomes: The qubit's cut-point outcomes, in any order (at least one).
    """
    return (sum(outcomes) + len(outcomes) - 1) % 2


def _folded_mid_samples(samples, collapse_order):
    """``{qubit: folded per-shot cut-point bit}`` from one sample array per collapse."""
    per_qubit = {}
    for qubit, sample in zip(collapse_order, samples):
        per_qubit.setdefault(qubit, []).append(np.asarray(sample).ravel())
    return {q: _fold_cut_outcomes(arrays) for q, arrays in per_qubit.items()}


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
    # Every collapse, in circuit order: a qubit cut twice is measured twice.
    collapse_order = [
        qubits[0] for gate, qubits, _params in ops if gate.upper() == "MEASURE"
    ]

    def qfunc():
        mids = []
        for gate, qubits, params in ops:
            value = _apply_pennylane_gate(gate, qubits, params)
            if value is not None:  # a MEASURE returns its MeasurementValue
                mids.append(value)
        # Finals for every qubit, then each cut-point outcome -- one execution keeps
        # them correlated per shot (built here so the MeasurementValues are in scope).
        return [qml.sample(wires=range(num_qubits))] + [
            qml.sample(value) for value in mids
        ]

    qnode = _build_pennylane_qnode(qfunc, num_qubits, shots, device, seed)
    out = qnode()
    final = np.asarray(out[0]).reshape(shots, num_qubits)  # (shots, n), qubit-0-first
    mid_samples = _folded_mid_samples(out[1:], collapse_order)

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
    subcirc, shots=1024, prep_qubits=None, device="default.qubit", seed=None
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
    :param seed: Seed for the device's sampling.
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
            circuit, observables, shots, set(prep_qubits), device, seed=seed
        )

    freq = _counts_pennylane(circuit, shots, device, seed=seed)
    return _expectation_from_counts_pennylane(freq, observables, shots)


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
    n_cuts,
    measure_index,
    prepare_index,
    mixed_index,
    plain_positions,
    measure_cuts,
    prepare_cuts,
    mixed_cuts,
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
        measure_tables = {}
        for (frag_idx, _setting), position in measure_index.items():
            measure_tables.setdefault(frag_idx, []).append(values[position])
        measure_values = {f: _merge_tables(t) for f, t in measure_tables.items()}
        prepare_values = {}
        for (frag_idx, combo), position in prepare_index.items():
            prepare_values.setdefault(frag_idx, {})[combo] = values[position]
        mixed_tables = {}
        for (frag_idx, _setting, prep_combo), position in mixed_index.items():
            mixed_tables.setdefault(frag_idx, {}).setdefault(prep_combo, []).append(
                values[position]
            )
        mixed_values = {
            f: {prep: _merge_tables(t) for prep, t in by_prep.items()}
            for f, by_prep in mixed_tables.items()
        }
        plain_factor = 1.0
        for position in plain_positions:
            plain_factor *= values[position]

        total = 0.0
        for index in range(8**n_cuts):
            settings = [(index // (8 ** (n_cuts - 1 - k))) % 8 for k in range(n_cuts)]
            sign = -1.0 if sum(s in _WIRE_SIGN_FLIP for s in settings) % 2 else 1.0
            term = sign
            for frag_idx, cuts in enumerate(measure_cuts):
                # Each measure circuit yields a {pauli: value} table (a Z-basis one
                # gives <Z> and <I>, a joint one a whole commuting group).
                pauli = tuple(_SETTING_PAULI[settings[cut - 1]] for cut in cuts)
                term *= measure_values[frag_idx][pauli]
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
            # Mixed fragments own cuts of BOTH types: basis_combo/extraction lookup
            # (like a measure fragment) and the |-> / |-i> expansion (like a prepare
            # fragment) both apply to the SAME fragment's joint table -- see
            # _reconstruct_optimized_wire_cut for the all-in-one counterpart.
            for frag_idx, (m_cuts, p_cuts) in enumerate(mixed_cuts):
                pauli = tuple(_SETTING_PAULI[settings[cut - 1]] for cut in m_cuts)
                expansions = [_WIRE_PREP[settings[cut - 1]] for cut in p_cuts]
                fragment_value = 0.0
                for choice in itertools.product(*expansions):
                    coefficient = 1.0
                    base_combo = []
                    for weight, base_state in choice:
                        coefficient *= weight
                        base_combo.append(base_state)
                    fragment_value += (
                        coefficient * mixed_values[frag_idx][tuple(base_combo)][pauli]
                    )
                term *= fragment_value
            total += term * plain_factor
        return total / (2**n_cuts)

    return reconstruct


def _optimized_wire_cutting_subcircuits(dag, cut, num_qubits, software, joint=False):
    """Three-step counterpart of :func:`_wire_cutting_optimized`.

    Returns a :class:`~Qdislib.utils.backends.CutSubcircuits` holding the reduced
    fragment subcircuits -- ``3**n`` measure probes (each an
    :class:`_OptimizedWireMeasure` returning all its P/I observables from one
    circuit), ``4**n`` base preparations, and (for chain topologies) mixed-fragment
    probes (each an :class:`_OptimizedWireMixed`), whose ``reconstruct`` does the
    fragment-separated base-8 sum. Every subcircuit is a normal ``expectation_value``
    target, so the standard three-step loop is unchanged. Returns ``None`` only
    when a cut's marker is missing entirely (a bookkeeping inconsistency).
    """
    if software not in _sdk_names():
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(_sdk_names())}."
        )
    wire_cuts = _expand_wire_cuts(dag, cut)
    fragments = _optimized_wire_fragments(dag, wire_cuts, num_qubits)
    if fragments is None:
        return None
    measure_frags, prepare_frags, mixed_frags, plain_frags = fragments

    subcircuits = []
    measure_index, prepare_index, mixed_index, plain_positions = {}, {}, {}, []
    for frag_idx, (sub, cuts, cut_qubits) in enumerate(measure_frags):
        # One _OptimizedWireMeasure per basis combo (3**n) -- evaluating it returns
        # every P/I observable of that circuit, so <Z> and <I> share one circuit.
        for basis_combo in _measurement_settings(len(cuts), joint):
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
    for frag_idx, (sub, m_cuts, cut_qubits, p_cuts) in enumerate(mixed_frags):
        # One _OptimizedWireMixed per (basis_combo, prep_combo) pair --
        # 3**(#m_cuts) * 4**(#p_cuts) probes for this fragment.
        for basis_combo in _measurement_settings(len(m_cuts), joint):
            for prep_combo in itertools.product(_PREP_STATES, repeat=len(p_cuts)):
                mixed_index[(frag_idx, basis_combo, prep_combo)] = len(subcircuits)
                subcircuits.append(
                    _OptimizedWireMixed(
                        sub.copy(),
                        m_cuts,
                        cut_qubits,
                        basis_combo,
                        p_cuts,
                        prep_combo,
                        software,
                    )
                )
    for sub in plain_frags:
        plain_positions.append(len(subcircuits))
        subcircuits.append(_make_subcircuit_task(sub.copy(), software))

    reconstructor = _make_optimized_wire_reconstructor(
        len(wire_cuts),
        measure_index,
        prepare_index,
        mixed_index,
        plain_positions,
        [cuts for _s, cuts, _q in measure_frags],
        [cuts for _s, cuts in prepare_frags],
        [(m_cuts, p_cuts) for _s, m_cuts, _q, p_cuts in mixed_frags],
    )
    return CutSubcircuits(
        subcircuits,
        number_cuts=len(wire_cuts),
        number_components=len(subcircuits),
        reconstructor=reconstructor,
    )


def _wire_cutting_subcircuits(
    rand_qc: typing.Any,
    cut: typing.List[typing.Any],
    observables=None,
    software="qiskit",
    optimized: bool = False,
    joint: bool = True,
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

    observables : str or list of (float, str), optional
        Per-qubit observable string (e.g. ``"ZZ"``). When given, the measurement
        basis is rotated accordingly before cutting. Defaults to the all-``Z``
        observable. A weighted Pauli sum ``[(coef, "PAULI"), ...]`` gives one
        subcircuit set per term, and ``reconstruct`` returns the weighted sum.

    software : str, optional
        The backend the subcircuits are emitted for, which must match the backend
        used to evaluate them:
        - ``"qiskit"`` (default): Qiskit ``QuantumCircuit`` objects
        - ``"qibo"``: Qibo ``Circuit`` objects
        - ``"cudaq"``: CUDA-Q ``(dag, num_qubits)`` payloads

    joint : bool, default=True
        With ``optimized=True``, measure the cut wires each fragment sends together
        instead of one by one: ``2**k + 1`` measurement circuits for ``k`` wires
        leaving a fragment instead of ``3**k`` (5 instead of 9 at ``k = 2``, 9
        instead of 27 at ``k = 3``), reading every Pauli of the wires exactly once.
        Same value; usually fewer shots for the same error, not always for
        ``k = 2``. Some of those circuits add CZ gates between the fragment's cut
        qubits: pass ``False`` to measure each wire on its own (no entangling gates),
        e.g. on hardware with costly two-qubit gates. Ignored without ``optimized``.

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

    # A weighted Pauli sum: one subcircuit set per term (each term rotates its own
    # measurement basis), recombined as sum_k c_k <P_k>.
    if observables is not None and not isinstance(observables, str):
        terms = list(normalize_observables(observables))
        parts = [
            _wire_cutting_subcircuits(
                rand_qc.copy() if isinstance(rand_qc, networkx.Graph) else rand_qc,
                cut,
                observables=pauli,
                software=software,
                optimized=optimized,
                joint=joint,
            )
            for _coeff, pauli in terms
        ]
        return CombinedCutSubcircuits(parts, [coeff for coeff, _ in terms])

    dag = convert_to_dag(rand_qc, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    _validate_wire_cuts_partition(dag, _expand_wire_cuts(dag, cut), _max_qubit(dag))

    if optimized and cut:
        # Fragment-separated subcircuits (~4**n + 4**n instead of 8**n * components);
        # the returned object's reconstruct() does the classical base-8 sum. Falls
        # back to the standard generator on error.
        try:
            optimized_subcircuits = _optimized_wire_cutting_subcircuits(
                dag.copy(), cut, _max_qubit(dag), software, joint=joint
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
        _warn_not_optimized("subcircuits", reason)

    if networkx.number_connected_components(dag.to_undirected(as_view=True)) > 1:
        # Disconnected pieces are cut (or not) independently and multiply: keep each
        # piece's subcircuits as its own part so reconstruct() can't mix them up.
        S = [
            dag.subgraph(c).copy()
            for c in networkx.connected_components(dag.to_undirected(as_view=True))
        ]
        parts = []
        for s in S:
            num_qubits = _max_qubit(s)
            tmp_cuts = [c for c in cut if s.has_edge(*c)]
            piece_cuts = len(_expand_wire_cuts(s, tmp_cuts)) if tmp_cuts else 0
            if tmp_cuts:
                piece = _generate_wire_cutting_subcircuits(
                    s, tmp_cuts, num_qubits=num_qubits, software=software
                )
            else:
                s_new, highest_qubit = _update_qubits(s)
                piece = [make_subcircuit(s_new, highest_qubit, software)]
            parts.append(
                CutSubcircuits(
                    piece,
                    number_cuts=piece_cuts,
                    number_components=max(1, len(piece) // (8**piece_cuts)),
                    reconstructor=_wire_cutting_subcircuit_reconstruction,
                )
            )
        return CombinedCutSubcircuits(parts)
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
            subcircuits = [make_subcircuit(s_new, highest_qubit, software)]

    number_components = max(1, len(subcircuits) // (8**number_cuts))
    return CutSubcircuits(
        subcircuits,
        number_cuts=number_cuts,
        number_components=number_components,
        reconstructor=_wire_cutting_subcircuit_reconstruction,
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

    _place_wire_cuts(dag, wire_cuts, num_qubits)

    if software not in _sdk_names():
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(_sdk_names())}."
        )

    # Each subcircuit is *created* in its own task (_make_subcircuit_task), separate
    # from the per-subcircuit evaluation task -- so under PyCOMPSs these are
    # Subcircuit futures fed straight into expectation_value's eval task (task -> task).
    subcircuits = []
    for index in range(8 ** len(wire_cuts)):
        copy_graph = _remove_red_edges(dag)

        num_components = networkx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
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


def _wire_cutting_subcircuit_reconstruction(
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
        >>> _wire_cutting_subcircuit_reconstruction(values, number_cuts=1, number_components=2)
    """
    from Qdislib.core.cutting_algorithms._reconstruction import signed_terms_sum

    return signed_terms_sum(
        expected_values, number_cuts, number_components, _change_sign
    )


# --------------------------------------------------------------------------- #
# Public API (the calling convention shared by every cut type, see _api.py)
# --------------------------------------------------------------------------- #


def _exact_or_legacy(legacy, custom_strategy=False):
    from Qdislib.core.cutting_algorithms.gate_cutting import _route_to_circuit_cutting

    _route_to_circuit_cutting(legacy, custom_strategy)


@standard_arguments
def wire_cutting(
    circuit,
    cuts,
    observable=None,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    optimized=False,
    joint=True,
    num_samples=None,
    cache=None,
    cache_key=None,
    cache_refresh=False,
    error=False,
    total_shots=None,
    circuits_per_task=None,
    parameter_values=None,
    cache_mode=None,
    _legacy=None,
):
    """Reconstruct ``<observable>`` from a circuit cut at wires (8 terms per wire,
    ``gamma = 4`` per wire).

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Wires to cut, as ``(source, target)`` gate pairs; each is cut on
        every qubit the two gates share.
    :param observable: A Pauli string or a weighted Pauli sum; ``None`` is all ``Z``.
    :param backend: Where the subcircuits run (see :func:`gate_cutting`).
    :param shots: Shots per subcircuit; ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling.
    :param sync: ``False`` returns the value as a future (PyCOMPSs).
    :param optimized: Run each fragment once per measurement basis / prepared state
        of its own wires (``3**i * 4**j`` per fragment instead of ``8**n`` overall).
    :param joint: With ``optimized``, measure the wires a fragment sends together
        (``2**k + 1`` circuits instead of ``3**k``).
    :param cache: Reuse fragment circuits already run (see :func:`circuit_cutting`
        and :func:`~Qdislib.use_cache`).
    :param cache_key: ``"exact"`` (default) or ``"semantic"``.
    :param cache_refresh: Run every circuit again and replace what is stored.
    :param error: ``True`` (or a number of resamples): return an
        :class:`~Qdislib.Estimate` with its standard error (see
        :func:`circuit_cutting`).
    :param total_shots: One shot budget instead of ``shots`` per circuit (see
        :func:`circuit_cutting`).
    :param circuits_per_task: Fragment circuits per task instead of one each (see
        :func:`circuit_cutting`).
    :param parameter_values: A list of value sets for the circuit's free
        parameters, one result each, the cut planned once (see
        :func:`circuit_cutting`).
    :param cache_mode: Deprecated and ignored (the cache keeps the whole result).
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :raises InvalidCutError: For a cut that does not join two gates of the circuit on
        a shared qubit, is given twice or does not separate its two sides.
    :raises ValueError: For an observable without one letter per qubit, ``shots`` below
        1, or ``shots=None`` with a backend that is no simulator.
    :return: The reconstructed expectation value (a list of them with
        ``parameter_values``).

    Example::

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)  # gate 2, "CZ_2"
        qc.ry(0.3, 1)  # gate 3, "RY_3"
        qc.cz(1, 2)  # gate 4, "CZ_4"
        value = qd.wire_cutting(qc, [("RY_3", "CZ_4")], "ZZZ", shots=4096, seed=1)
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _cache_mode_option,
        _legacy_run,
        _legacy_store,
    )

    custom = _cache_mode_option(cache_mode)
    if not _legacy_run(_legacy, custom) or shots is None or num_samples is not None:
        from Qdislib.core.cutting_algorithms.mixed_cutting import circuit_cutting

        _exact_or_legacy(_legacy, custom)
        return circuit_cutting(
            circuit,
            {"wire": cuts},
            observable,
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            optimized=optimized,
            joint=joint,
            num_samples=num_samples,
            cache=cache,
            cache_key=cache_key,
            cache_refresh=cache_refresh,
            error=error,
            total_shots=total_shots,
            circuits_per_task=circuits_per_task,
            parameter_values=parameter_values,
        )
    if (
        error
        or total_shots is not None
        or circuits_per_task is not None
        or parameter_values is not None
    ):
        raise ValueError(
            "error=, total_shots=, circuits_per_task= and parameter_values= run through circuit_cutting: they are not available with the "
            "deprecated qpu=/gpu=/method= options."
        )
    legacy = _legacy or {}
    return _wire_cutting(
        circuit,
        cuts,
        sync=sync,
        observables=checked_observable(observable, circuit),
        shots=shots,
        backend=legacy.get("platform", "numpy"),
        qpu=legacy.get("qpu"),
        software=resolve_backend(backend, circuit),
        cache=_legacy_store(cache),
        cache_mode=cache_mode or "expectation",
        optimized=optimized,
        joint=joint,
        seed=seed,
    )


@standard_arguments
def wire_cutting_sampling(
    circuit,
    cuts,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    optimized=False,
    joint=True,
    num_samples=None,
    cache=None,
    cache_key=None,
    cache_refresh=False,
    error=False,
    circuits_per_task=None,
    parameter_values=None,
    qubits=None,
    _legacy=None,
):
    """Reconstruct the output **probability distribution** of a circuit cut at
    wires.

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Wires to cut, as in :func:`wire_cutting`.
    :param backend: Where the subcircuits run (see :func:`gate_cutting`).
    :param shots: Shots per subcircuit; ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling.
    :param sync: ``False`` returns the distribution as a future (PyCOMPSs).
    :param optimized: Run each fragment once per measurement basis / prepared state
        of its own wires (``3**i * 4**j`` per fragment instead of ``8**n`` overall).
    :param joint: With ``optimized``, measure the wires a fragment sends together
        (``2**k + 1`` circuits instead of ``3**k``).
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :param cache: Reuse fragment circuits already run (see :func:`circuit_cutting`
        and :func:`~Qdislib.use_cache`).
    :param cache_key: ``"exact"`` (default) or ``"semantic"``.
    :param cache_refresh: Run every circuit again and replace what is stored.
    :param error: ``True`` (or a number of resamples): return a
        :class:`~Qdislib.Distribution` with each probability's standard error.
    :param circuits_per_task: Fragment circuits per task instead of one each (see
        :func:`circuit_cutting`).
    :param parameter_values: A list of value sets for the circuit's free
        parameters, one distribution each, the cut planned once.
    :param qubits: The distribution over only these qubits, as in
        :func:`circuit_cutting_sampling`.
    :raises InvalidCutError: For a cut that does not join two gates of the circuit on
        a shared qubit, is given twice or does not separate its two sides.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings), a
        :class:`~Qdislib.FactoredDistribution` over more than 24 qubits; a list of them
        with ``parameter_values``.
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _legacy_run

    if not _legacy_run(_legacy) or shots is None or num_samples is not None:
        from Qdislib.core.cutting_algorithms.mixed_cutting import (
            circuit_cutting_sampling,
        )

        _exact_or_legacy(_legacy)
        return circuit_cutting_sampling(
            circuit,
            {"wire": cuts},
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            optimized=optimized,
            joint=joint,
            num_samples=num_samples,
            cache=cache,
            cache_key=cache_key,
            cache_refresh=cache_refresh,
            error=error,
            circuits_per_task=circuits_per_task,
            parameter_values=parameter_values,
            qubits=qubits,
        )
    if (
        error
        or circuits_per_task is not None
        or parameter_values is not None
        or qubits is not None
    ):
        raise ValueError(
            "error=, circuits_per_task=, parameter_values= and qubits= run through circuit_cutting: they are not available with the "
            "deprecated qpu= / qpu_dict= options."
        )
    legacy = _legacy or {}
    return _wire_cutting_sampling(
        circuit,
        cuts,
        shots=shots,
        software=resolve_backend(backend, circuit),
        sync=sync,
        qpu=legacy.get("qpu", False),
        qpu_dict=legacy.get("qpu_dict"),
        seed=seed,
        optimized=optimized,
        joint=joint,
    )


@standard_arguments
def wire_cutting_subcircuits(
    circuit,
    cuts,
    observable=None,
    *,
    _software=None,
    optimized=False,
    joint=True,
    num_samples=None,
    seed=None,
):
    """The subcircuits of :func:`wire_cutting`, to run yourself (three-step)::

        subs = qd.wire_cutting_subcircuits(circuit, [("RY_7", "CZ_8")], "ZZZZ")
        values = [qd.expectation_value(s, shots=4096) for s in subs]
        value = subs.reconstruct(values)

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Wires to cut, as in :func:`wire_cutting`.
    :param observable: A Pauli string or a weighted Pauli sum; ``None`` is all ``Z``.
    :param optimized: Run each fragment once per measurement basis / prepared state
        of its own wires.
    :param joint: With ``optimized``, measure the wires a fragment sends together.
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :param seed: Seed for the ``num_samples`` draws.
    :raises InvalidCutError: For a cut that does not join two gates of the circuit on
        a shared qubit, is given twice or does not separate its two sides.
    :return: A :class:`~Qdislib.CutCircuits` (a list of subcircuits with
        ``reconstruct(values)``); for a Pauli sum a
        :class:`~Qdislib.CombinedCutSubcircuits`.
    """
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        circuit_cutting_subcircuits,
    )

    return circuit_cutting_subcircuits(
        circuit,
        {"wire": cuts},
        observable,
        _software=_software,
        optimized=optimized,
        joint=joint,
        num_samples=num_samples,
        seed=seed,
    )


@standard_arguments
def wire_cutting_sampling_subcircuits(
    circuit,
    cuts,
    *,
    _software=None,
    optimized=False,
    joint=True,
    num_samples=None,
    seed=None,
):
    """The subcircuits of :func:`wire_cutting_sampling`, to sample yourself::

        subs = qd.wire_cutting_sampling_subcircuits(circuit, [("RY_7", "CZ_8")])
        results = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(results)

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Wires to cut, as in :func:`wire_cutting`.
    :param optimized: Run each fragment once per measurement basis / prepared state
        of its own wires.
    :param joint: With ``optimized``, measure the wires a fragment sends together.
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :param seed: Seed for the ``num_samples`` draws.
    :raises InvalidCutError: For a cut that does not join two gates of the circuit on
        a shared qubit, is given twice or does not separate its two sides.
    :return: A :class:`~Qdislib.CutCircuits` (a list of subcircuits with
        ``reconstruct(results)``).
    """
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        circuit_cutting_sampling_subcircuits,
    )

    return circuit_cutting_sampling_subcircuits(
        circuit,
        {"wire": cuts},
        _software=_software,
        optimized=optimized,
        num_samples=num_samples,
        seed=seed,
        joint=joint,
    )


# Deprecated in 2.0.0: ``subcircuits.reconstruct(results)`` is the one way.


def wire_cutting_subcircuit_reconstruction(
    expected_values, number_cuts, number_components=2
):
    """Deprecated: use ``subcircuits.reconstruct(results)``.

    :param expected_values: One value per subcircuit of a 1.x wire cut, in order.
    :param number_cuts: The wires cut.
    :param number_components: The subcircuits of each of its ``8**number_cuts``
        configurations.
    :raises TypeError: For the results of a 2.0 list, which only its ``reconstruct``
        can combine.
    :return: The reconstructed expectation value.
    """
    import warnings

    warnings.warn(
        "wire_cutting_subcircuit_reconstruction() is deprecated: use subcircuits.reconstruct(results); it is "
        "removed in the next release.",
        DeprecationWarning,
        stacklevel=2,
    )
    from Qdislib.core.cutting_algorithms._reconstruction import reject_engine_results

    reject_engine_results(expected_values, "wire_cutting_subcircuit_reconstruction")
    return _wire_cutting_subcircuit_reconstruction(
        expected_values, number_cuts, number_components
    )
