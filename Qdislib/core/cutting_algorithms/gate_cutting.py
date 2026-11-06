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

"""Gate cutting algorithms."""

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "gate_cutting",
    "gate_cutting_sampling",
    "gate_cutting_sampling_subcircuits",
    "gate_cutting_sampling_reconstruction",
    "sample_subcircuit",
    "SampledCutSubcircuits",
    "gate_cutting_subcircuits",
    "gate_cutting_subcircuit_reconstruction",
]

import networkx as nx
import cmath
import math
import numpy as np

try:
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit_aer import AerSimulator
except ImportError:  # qiskit / qiskit-aer are optional
    generate_preset_pass_manager = AerSimulator = None

import re

from Qdislib.utils.sdks import sdk_names as _sdk_names
from Qdislib.utils.gates import canonical
from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN


from Qdislib.utils.graph_qibo import (
    _update_qubits,
    _remove_red_edges,
    _max_qubits_graph,
)
from Qdislib.core.cutting_algorithms.wire_cutting import (
    _fold_cut_outcomes,
    _folded_mid_samples,
    _sum_results,
)
from Qdislib.core.find_cut.find_cut import _find_nodes_with_qubit
from Qdislib.utils.graph_pennylane import _apply_pennylane_gate, _build_pennylane_qnode

try:
    from qiskit_ibm_runtime import SamplerV2 as Sampler
except ImportError:  # qiskit-ibm-runtime optional (only used for IBM QPU)
    Sampler = None


try:
    import qiskit
    from qiskit import transpile

    qiskit_available = True
except ImportError:  # qiskit is optional; functions that need it fail when called
    qiskit = transpile = None
    qiskit_available = False

try:
    import qibo

    qibo_available = True
except ImportError:  # qibo is optional; functions that need it fail when called
    qibo = None
    qibo_available = False

try:
    import pennylane as qml

    pennylane_available = True
except ImportError:  # pennylane is optional; functions that need it fail when called
    qml = None
    pennylane_available = False


import typing


from Qdislib.utils.backends import (
    convert_to_dag,
    change_basis_dag,
    normalize_observables,
    make_subcircuit,
    SampledSubcircuit,
    _make_subcircuit_task,
    CutSubcircuits,
    CombinedCutSubcircuits,
)
from Qdislib.utils.exceptions import UnsupportedBackendError, InvalidCutError
from Qdislib.core.cutting_algorithms._reconstruction import next_seed
from Qdislib.core.cutting_algorithms._api import (
    standard_arguments,
    resolve_backend,
    default_observable,
)


@task(returns=1, values=COLLECTION_IN)
def _weighted_sum(coeffs, values):
    """Reduce a Hamiltonian's terms ``<sum_i c_i P_i> = sum_i c_i <P_i>`` (COMPSs task).

    Fed the per-term expectation futures directly via ``COLLECTION_IN`` (task ->
    task), so every ``(term x cut-configuration)`` subcircuit schedules together
    under PyCOMPSs and this reduction fires only once they finish -- no master-side
    barrier between Hamiltonian terms.
    """
    return float(sum(c * v for c, v in zip(coeffs, values)))


def _gate_cutting(
    dag: typing.Any,
    gates_cut: typing.List[typing.Any],
    observables: typing.Optional[typing.Any] = None,
    shots: int = 1024,
    method: str = "automatic",
    sync: bool = True,
    gpu: bool = False,
    gpu_min_qubits: typing.Optional[int] = None,
    qpu: bool = False,
    qpu_dict: typing.Optional[dict] = None,
    qpu_config: typing.Optional[dict] = None,
    software: str = "qiskit",
    cache=None,
    cache_mode="expectation",
    seed=None,
):
    """
    Apply gate cutting to a quantum circuit to enable distributed smaller executions of subcircuits.

    This function partitions a quantum circuit by cutting a list of specified two-qubit gates,
    evaluates the resulting subcircuits individually (optionally on CPU/GPU or QPU resources),
    and reconstructs the final expectation value from the partial results.
    It supports circuits represented in Qiskit or Qibo formats, as well as their internal DAG representations.

    If no gates are cut, the function directly computes the expectation value of the original circuit.

    Parameters
    ----------
    dag : qiskit.QuantumCircuit or qibo.models.Circuit or networkx.Graph
        The quantum circuit to cut and evaluate. Can be a Qiskit circuit, a Qibo circuit,
        or an already constructed DAG (as a NetworkX graph).
    gates_cut : list
        List of two-qubit gates (edges in the DAG) to cut. Each gate is typically identified by its node ID.
    observables : list or dict, optional
        Observables to measure. If specified, the circuit is transformed into the measurement basis.
        Observables may include the identity ("I"), handled for any qubit on all
        three CPU backends (``"qiskit"``, ``"qibo"``, ``"cudaq"``; see Notes). The
        GPU/QPU evaluation paths and the 3-step ``gate_cutting_subcircuits``
        workflow do not yet decouple this case and raise
        :class:`~Qdislib.utils.exceptions.InvalidCutError` for an identity on a cut
        qubit that feeds a multi-qubit gate after the cut.
    shots : int, default=1024
        Number of measurement shots to use for circuit evaluations.
    method : str, default='automatic'
        Simulation method to use. Possible values are backend-dependent (e.g., "statevector", "qasm_simulator").
    sync : bool, default=True
        Whether to synchronize and wait for all computations to finish (useful in distributed environments like PyCOMPSs).
    gpu : bool, default=False
        If True, attempt to run circuit evaluations on GPU simulators when possible.
    gpu_min_qubits : int, optional
        Minimum number of qubits needed to offload execution to GPU simulators. Default is 0.
    qpu : bool, default=False
        If True, attempt to run circuit evaluations on real Quantum Processing Units (QPUs).
    qpu_dict : dict, optional
        Dictionary specifying QPU backends and their maximum supported qubits.
        Example: `{"IBM_Quantum": 153, "MN_Ona": 5}`.
    qpu_config : dict, optional
        Explicit IBM Quantum configuration (preferred over `IBM_QUANTUM_*`
        environment variables, which remain the fallback). Keys: ``token``,
        ``channel`` (``"ibm_cloud"``/``"ibm_quantum"``), ``instance``,
        ``backend_name`` (omit for the least-busy backend), ``mode``
        (``"job"``, the default: plain jobs; ``"batch"`` / ``"session"``),
        ``max_time`` (batch / session only), ``proxies``, or a pre-built
        ``service`` (``QiskitRuntimeService``). The
        connection and backend selection run once (see
        :class:`~Qdislib.core.backend.qpu_backends.IBMBackend`).
    software : str, default="qiskit"
        Backend used to build and execute the subcircuits when no QPU/GPU
        applies: "qiskit", "qibo" or "cudaq".
    cache : CircuitCache, optional
        Enable the semantic subcircuit cache (see
        :func:`~Qdislib.utils.circuit_cache.default_circuit_cache`). When given,
        evaluation runs the CPU-simulation subcircuit workflow (GPU/QPU dispatch
        is bypassed). ``None`` (default) disables caching.
    cache_mode : str or CacheStrategy, default="expectation"
        Cache granularity; ``"expectation"``, ``"counts"``, or a custom strategy.

    Diagnostic details are logged at DEBUG level; enable them with
    ``qd.set_log_level("DEBUG")``.

    Returns
    -------
    final_recons : float
        The reconstructed expectation value after gate cutting and execution of subcrcuits.

    Notes
    -----
    - If the circuit consists of multiple disconnected components, each component is processed separately.
    - Depending on the device availability and size of the subcircuits, they are evaluated on CPU, GPU, or QPU.
    - Results are automatically scaled by the number of gates cut (following standard gate cutting reconstruction rules).
    - GPU and QPU support is optional and requires appropriate setup (e.g., Qiskit Aer with GPU support, IBM Q account credentials).
    - Identity observable on a cut qubit: a gate cut measures its cut qubits as
      part of the measure-and-prepare decomposition, and that measurement carries
      the reconstruction sign. If the cut qubit feeds a multi-qubit gate *after*
      the cut, an identity there would suppress it if read off the final
      measurement. All three CPU backends instead recover the sign from the
      cut-point measurement (``qiskit``: the cut-measure classical register;
      ``qibo``: the collapse gate's per-shot samples; ``cudaq``:
      ``sample(explicit_measurements=True)``), decoupling it from the final
      observable, so identity on any cut qubit reconstructs correctly. The
      GPU/QPU paths and the 3-step subcircuit workflow do not do this yet and
      raise :class:`~Qdislib.utils.exceptions.InvalidCutError` rather than
      returning a silently-wrong value.

    Examples
    --------
    Cutting a simple Qiskit circuit:

    .. code-block:: python

        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)     # Gate "H_1"
        qc.cz(0, 1) # Gate "CZ_2"
        qc.h(1)     # Gate "H_3"

        # Define the gate to cut (e.g., the CZ gate)
        gates_to_cut = ["CZ_2"]  # You must identify the gate in order of definition

        reconstruction = _gate_cutting(qc, gates_to_cut, shots=2048)

    """

    # A weighted sum of Pauli observables (a Hamiltonian) is evaluated term by
    # term: <sum_i c_i P_i> = sum_i c_i <P_i>. Each term reuses the single
    # Pauli-string path below. A plain string (or None) skips this and runs once.
    if observables is not None and not isinstance(observables, str):
        # <sum_i c_i P_i> = sum_i c_i <P_i>. Submit every term's cut WITHOUT waiting
        # (sync=False) so all (term x cut-configuration) subcircuits schedule
        # together under PyCOMPSs, then reduce with one task fed by the term futures
        # (task -> task, no master barrier between terms).
        terms = list(normalize_observables(observables))
        coeffs = [coeff for coeff, _ in terms]
        values = [
            _gate_cutting(
                dag.copy() if isinstance(dag, nx.Graph) else dag,
                gates_cut,
                observables=pauli,
                shots=shots,
                method=method,
                sync=False,
                gpu=gpu,
                gpu_min_qubits=gpu_min_qubits,
                qpu=qpu,
                qpu_dict=qpu_dict,
                qpu_config=qpu_config,
                software=software,
                cache=cache,
                cache_mode=cache_mode,
                seed=next_seed(seed, "term", k),
            )
            for k, (_, pauli) in enumerate(terms)
        ]
        total = _weighted_sum(coeffs, values)
        return compss_wait_on(total) if sync else total

    # Semantic cache: reuse the recommended subcircuit workflow, which routes
    # each subcircuit through expectation_value(cache=...). This runs the CPU
    # simulation path, so GPU/QPU dispatch is bypassed when caching is enabled.
    if cache is not None:
        from Qdislib.utils.backends import expectation_value

        subcircuits = _gate_cutting_subcircuits(
            dag, gates_cut, observables=observables, software=software
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

    # Normalise optional dict/int args before they are indexed below.
    if qpu_dict is None:
        qpu_dict = {}
    if gpu_min_qubits is None:
        gpu_min_qubits = 0

    if qpu and "IBM_Quantum" in qpu_dict:
        batch, backend = _check_ibm_qc(qpu_config)
    else:
        batch, backend = None, None

    dag = convert_to_dag(dag, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    # A non-CZ/CX cut gate is decomposed into CZ + single-qubit gates and its CZs are
    # cut instead (transparently grows gates_cut, so 6**n / 1/2**n scale accordingly).
    dag, gates_cut = _expand_general_gate_cuts(dag, gates_cut)
    _validate_gate_cuts_partition(dag, gates_cut)

    if nx.number_connected_components(dag.to_undirected(as_view=True)) > 1:
        S = [
            dag.subgraph(c).copy()
            for c in nx.connected_components(dag.to_undirected(as_view=True))
        ]
        results = []
        for piece, s in enumerate(S):
            tmp_cuts = [c for c in gates_cut if s.has_node(c)]
            if tmp_cuts:
                graphs = _execute_gate_cutting(
                    s,
                    tmp_cuts,
                    shots=shots,
                    method=method,
                    gpu=gpu,
                    gpu_min_qubits=gpu_min_qubits,
                    qpu=qpu,
                    qpu_dict=qpu_dict,
                    batch=batch,
                    backend=backend,
                    software=software,
                    seed=next_seed(seed, "piece", piece),
                )
                graphs = _sum_results(graphs)
                results.append(graphs)
            else:
                results.append(
                    _build_and_run_gate_subcircuit(
                        s,
                        _max_qubits_graph(s),
                        software,
                        shots,
                        method,
                        gpu,
                        gpu_min_qubits,
                        qpu,
                        qpu_dict,
                        batch,
                        backend,
                        next_seed(seed, "piece", piece),
                    )
                )
        from Qdislib.core.cutting_algorithms._reconstruction import scaled_product

        final_recons = scaled_product(results, 1 / (2 ** len(gates_cut)))
        return compss_wait_on(final_recons) if sync else final_recons
    else:
        if gates_cut:
            results = _execute_gate_cutting(
                dag,
                gates_cut,
                shots=shots,
                method=method,
                gpu=gpu,
                gpu_min_qubits=gpu_min_qubits,
                qpu=qpu,
                qpu_dict=qpu_dict,
                batch=batch,
                backend=backend,
                software=software,
                seed=seed,
            )
            from Qdislib.core.cutting_algorithms._reconstruction import scaled_sum

            final_recons = scaled_sum(results, 1 / (2 ** len(gates_cut)))
            return compss_wait_on(final_recons) if sync else final_recons
        else:
            final_recons = _build_and_run_gate_subcircuit(
                dag,
                _max_qubits_graph(dag),
                software,
                shots,
                method,
                gpu,
                gpu_min_qubits,
                qpu,
                qpu_dict,
                batch,
                backend,
                seed,
            )
            if sync:
                final_recons = compss_wait_on(final_recons)
            logger.debug("No cut performed; expected value: %s", final_recons)
        return final_recons


# The measure-and-prepare substitution below implements the CZ decomposition; a
# CX/CNOT is cut via the identity CX(c, t) = H(t) CZ(c, t) H(t) (see ``_generate_cut``).
# Any OTHER two-qubit gate is decomposed into CZ + single-qubit gates up front
# (``_expand_general_gate_cuts``), so only these reach the direct substitution.
_CUTTABLE_GATES = {"cz", "cx", "cnot"}


def _assert_cuttable_gate(dag, gate_name, rotations=False):
    """Guard the direct-substitution path: the SUBS insertion below only handles a
    two-qubit ``CZ`` or ``CX``/``CNOT``. General two-qubit gates ARE supported -- the
    public entry points call :func:`_expand_general_gate_cuts` first, which decomposes
    them into ``CZ`` + single-qubit gates -- so by the time a cut reaches here it must
    already be a ``CZ``/``CX``.

    :param dag: DAG the cut is applied to (after ``_remove_red_edges``).
    :param gate_name: Node id from ``gates_cut``.
    :raises InvalidCutError: if the node is not a two-qubit CZ or CX/CNOT.
    """
    gate = dag.nodes[gate_name].get("gate") or ""
    qubits = dag.nodes[gate_name].get("qubits") or ()
    allowed = _CUTTABLE_GATES | ({"rzz"} if rotations else set())
    if gate.lower() not in allowed or len(qubits) != 2:
        raise InvalidCutError(
            f"Gate cutting substitutes a two-qubit CZ or CX/CNOT (a CX is cut via "
            f"CX = H*CZ*H), but node {gate_name!r} is a {gate!r} on qubits "
            f"{tuple(qubits)}. General gates are decomposed to CZ automatically by the "
            f"public entry points; call _expand_general_gate_cuts first, or use "
            f"gate_cutting_hadamard (which cuts any gate directly)."
        )


def _generate_cut(dag, gates_cut, rotations=False):
    """Replace each cut gate by its two side markers ``SUBS1_k`` / ``SUBS2_k``.
    ``rotations`` also accepts ``rzz`` (the angle-aware cut of ``circuit_cutting``,
    whose terms the caller supplies)."""
    dag_copy = dag.copy()
    dag_copy = _remove_red_edges(dag_copy)
    for index, gate_name in enumerate(gates_cut, start=1):
        _assert_cuttable_gate(dag_copy, gate_name, rotations)
        # A CX/CNOT is cut as H(t) CZ(c, t) H(t): the CZ measure-and-prepare below
        # is reused unchanged, and the target qubit (qubits[1]) is conjugated by a
        # Hadamard on each side of its cut node -- so it is effectively measured
        # and prepared in the X basis, which is exactly what the identity says.
        is_cx = (dag_copy.nodes[gate_name].get("gate") or "").lower() in ("cx", "cnot")
        target_qubits = dag_copy.nodes[gate_name]["qubits"]
        pred_0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[0], direction="predecessor"
        )

        pred_1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[1], direction="predecessor"
        )

        succ_0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[0], direction="successor"
        )

        succ_1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=target_qubits[1], direction="successor"
        )

        dag_copy.remove_node(gate_name)

        dag_copy.add_node(
            f"SUBS1_{index}", gate="s", qubits=(target_qubits[0],), parameters=()
        )
        dag_copy.add_node(
            f"SUBS2_{index}", gate="s", qubits=(target_qubits[1],), parameters=()
        )

        if pred_0:
            dag_copy.add_edge(pred_0[0], f"SUBS1_{index}", color="blue")

        if succ_0:
            dag_copy.add_edge(f"SUBS1_{index}", succ_0[0], color="blue")

        if is_cx:
            # pred_1 -> H -> SUBS2 -> H -> succ_1 (Hadamard conjugation on target)
            h_pre = f"CXHPRE_{index}"
            h_post = f"CXHPOST_{index}"
            dag_copy.add_node(
                h_pre, gate="h", qubits=(target_qubits[1],), parameters=()
            )
            dag_copy.add_node(
                h_post, gate="h", qubits=(target_qubits[1],), parameters=()
            )
            if pred_1:
                dag_copy.add_edge(pred_1[0], h_pre, color="blue")
            dag_copy.add_edge(h_pre, f"SUBS2_{index}", color="blue")
            dag_copy.add_edge(f"SUBS2_{index}", h_post, color="blue")
            if succ_1:
                dag_copy.add_edge(h_post, succ_1[0], color="blue")
        else:
            if pred_1:
                dag_copy.add_edge(pred_1[0], f"SUBS2_{index}", color="blue")

            if succ_1:
                dag_copy.add_edge(f"SUBS2_{index}", succ_1[0], color="blue")

    return dag_copy


def _validate_gate_cuts_partition(dag, gates_cut):
    """Raise unless every gate cut separates the circuit (see
    :func:`~Qdislib.core.cutting_algorithms._reconstruction.check_cuts_separate`).

    Gate cutting (like wire cutting -- see
    :func:`~Qdislib.core.cutting_algorithms.wire_cutting._validate_wire_cuts_partition`)
    only pays off if it increases the number of independently executable
    pieces beyond what the circuit already had. A cut gate whose two qubits
    stay linked through some OTHER, uncut gate (e.g. a third gate bridging the
    same two qubits through a different path) leaves the graph exactly as
    connected as before, so reconstructing through it pays the cut's shot-count
    cost (``6**n`` configurations) for none of its qubit-width benefit -- running
    the uncut circuit directly is strictly better, cheaper, and less noisy.

    :param dag: The circuit DAG the cut will actually be applied to (i.e. after
        any general-gate expansion the caller has already done).
    :param gates_cut: Gate node ids to cut, matching what will be passed to
        :func:`_generate_cut`.
    :raises InvalidCutError: If a cut gate's two sides stay in one fragment.
    """
    if not gates_cut:
        return
    from Qdislib.core.cutting_algorithms._reconstruction import check_cuts_separate

    check_cuts_separate(
        _generate_cut(dag, gates_cut),
        {
            gate: (f"SUBS1_{k}", f"SUBS2_{k}")
            for k, gate in enumerate(gates_cut, start=1)
        },
    )


def _decompose_gate_to_cz(matrix):
    """Decompose a general two-qubit unitary into ``CZ`` + single-qubit ``RZ``/``RY``.

    Returns qibo-convention ops ``[(gate_name, (local_qubits...), [params]), ...]`` that
    reconstruct ``matrix`` (qibo big-endian, ``q0 (x) q1``) up to a global phase -- so a
    general gate can be cut in the measure-and-prepare method by cutting the ``CZ``\\ s of
    its decomposition (the same way the Hadamard cut reuses the KAK decomposition).
    Uses qiskit's ``TwoQubitBasisDecomposer`` with the ``CZ`` basis; qiskit is
    little-endian, so both qubit orderings are tried and the reconstructing one is kept.

    :raises RuntimeError: if no ordering reconstructs ``matrix`` (should not happen).
    """
    import qibo
    from qibo import gates as qibo_gates, models
    from qiskit.synthesis import TwoQubitBasisDecomposer
    from qiskit.circuit.library import CZGate

    matrix = np.asarray(matrix, dtype=complex)
    decomposer = TwoQubitBasisDecomposer(CZGate(), euler_basis="ZYZ")
    qc = decomposer(matrix)
    backend = qibo.get_backend()

    def _ops(swap):
        ops = []
        for instr in qc.data:
            qubits = [qc.find_bit(b).index for b in instr.qubits]
            if swap:  # qiskit little-endian -> qibo big-endian
                qubits = [1 - q for q in qubits]
            ops.append(
                (
                    instr.operation.name.upper(),
                    tuple(qubits),
                    [float(p) for p in instr.operation.params],
                )
            )
        return ops

    for swap in (True, False):  # the little->big-endian flip first (the usual case)
        ops = _ops(swap)
        circuit = models.Circuit(2)
        for name, qubits, params in ops:
            cls = getattr(qibo_gates, name)
            circuit.add(cls(*qubits, *params) if params else cls(*qubits))
        rebuilt = np.asarray(circuit.unitary(backend))
        if abs(np.trace(rebuilt.conj().T @ matrix)) / 4 > 1 - 1e-6:
            return ops
    raise RuntimeError("could not decompose the two-qubit gate into CZ + single-qubit")


def _decompose_gate_in_dag(dag, gate_name):
    """Replace a general two-qubit gate node with its ``CZ`` + ``RZ``/``RY``
    decomposition (mutates ``dag`` in place); return the new ``CZ`` node ids (the cut
    targets). The single-qubit gates ride along in whichever fragment their qubit lands
    in, exactly like any other gate on the cut qubits."""
    from Qdislib.utils.graph_qibo import gate_matrix

    node = dag.nodes[gate_name]
    params = node.get("parameters") or ()
    if len(node["qubits"]) != 2:
        raise InvalidCutError(
            f"Gate cutting needs a two-qubit gate; {gate_name!r} is not one"
        )
    matrix = gate_matrix(str(node["gate"]), params, 2)
    if matrix.shape != (4, 4):
        raise InvalidCutError(
            f"Gate cutting needs a two-qubit gate; {gate_name!r} is {matrix.shape}"
        )
    return _replace_with_ops(dag, gate_name, _decompose_gate_to_cz(matrix), "CZ")


def _expand_general_gate_cuts(dag, gates_cut):
    """Expand any non-``CZ``/``CX`` gate in ``gates_cut`` into the ``CZ``\\ s of its
    decomposition (mutating ``dag``). Returns ``(dag, expanded_gates_cut)`` where every
    entry is now a cuttable ``CZ``/``CX`` node -- so the rest of the measure-and-prepare
    pipeline (and its ``6 ** n`` enumeration) handles a general gate transparently.

    A generic two-qubit gate needs up to 3 ``CZ``\\ s (so up to 3 cuts on the same qubit
    pair). Reconstructing that correctly for an *identity* observable on a cut qubit
    needs the per-cut sign recovery, which the qibo evaluator does; the Qiskit/CUDA-Q
    evaluators record one cut-point bit per qubit, so for general gates observed with
    identity on a cut qubit use ``software="qibo"`` (or the ancilla-based
    :func:`gate_cutting_hadamard`, which handles any gate on any backend)."""
    expanded = []
    for name in gates_cut:
        gate = (dag.nodes[name].get("gate") or "").lower()
        if gate in ("cz", "cx", "cnot"):
            expanded.append(name)
        else:
            expanded.extend(_decompose_gate_in_dag(dag, name))
    return dag, expanded


# CPU-simulation backends whose expectation reconstruction decouples the internal
# cut-measurement sign from the final observable (see
# :func:`_expectation_from_counts_qiskit`), and so handle an identity observable
# on any cut qubit correctly. All three CPU backends do; the guard below still
# applies to GPU/QPU evaluation paths (not decoupled) and to the 3-step
# subcircuit workflow (which evaluates without the decoupled reconstruction).
class _Simulators:
    """Every SDK of ``Qdislib.utils.sdks`` (its simulator decouples the sign)."""

    def __contains__(self, software):
        return software in _sdk_names()


_IDENTITY_CUT_DECOUPLED = _Simulators()


def _assert_identity_observable_supported(cut_dag, software=None, gpu=False, qpu=False):
    """Reject identity-observable cut cases a path cannot yet reconstruct.

    A gate cut measures its cut qubits as part of the measure-and-prepare
    quasiprobability decomposition, and that measurement carries the
    reconstruction sign. If the cut qubit ends on single-qubit gates it stays
    isolated after the cut, so an identity observable on it is a trivial
    single-qubit subcircuit and always works. But if the cut qubit feeds a
    **multi-qubit** gate *after* the cut, an identity there would suppress the
    load-bearing measurement if it were read off the final measurement.

    The three CPU backends (:data:`_IDENTITY_CUT_DECOUPLED`) and the **GPU** path
    recover that sign from the cut-point measurement (decoupling it from the final
    observable) and handle this correctly -- both the all-in-one and three-step
    workflows. CPU reads the joint final/cut-point counts directly; GPU routes to
    Aer's GPU simulator and runs the *same* joint-counts reconstruction
    (:func:`_expec_value_qiskit_gpu`). So the check is skipped for CPU and GPU.

    It still fires for **QPU** paths: a real device returns only aggregated final
    counts (the per-shot mid-circuit cut-point outcome is lost), so a cut qubit
    that keeps evolving after the cut cannot be decoupled -- raising here beats
    returning a silently-wrong value.

    :param cut_dag: DAG after :func:`_generate_cut` (carries the ``SUBS*`` nodes).
    :param software: CPU-simulation backend that will evaluate the subcircuits.
    :param gpu: Whether GPU evaluation is in use (Aer GPU -- supported).
    :param qpu: Whether QPU evaluation is in use (unsupported -- aggregated counts).
    :raises InvalidCutError: if an identity observable sits on a cut qubit that
        feeds a multi-qubit gate after the cut, on the (QPU) path that can't
        decouple it.
    """
    # CPU decoupled backends handle it; so does GPU (it routes to Aer's GPU
    # simulator and runs the same joint-counts reconstruction). Only real QPUs
    # can't -- they return aggregated final counts with no per-shot cut-point.
    if not qpu and (gpu or software in _IDENTITY_CUT_DECOUPLED):
        return

    for subs in [n for n in cut_dag if isinstance(n, str) and n.startswith("SUBS")]:
        qubits = cut_dag.nodes[subs].get("qubits") or ()
        if not qubits:
            continue
        cut_qubit = qubits[0]
        descendants = nx.descendants(cut_dag, subs)
        on_this_qubit = any(
            cut_dag.nodes[d].get("gate") == "Observable I"
            and (cut_dag.nodes[d].get("qubits") or (None,))[0] == cut_qubit
            for d in descendants
        )
        if not on_this_qubit:
            continue
        feeds_multiqubit = any(
            len(cut_dag.nodes[d].get("qubits") or ()) > 1 for d in descendants
        )
        if feeds_multiqubit:
            raise InvalidCutError(
                f"This QPU gate-cut path cannot honour an identity ('I') observable "
                f"on cut qubit {cut_qubit}: it feeds a multi-qubit gate after the "
                f"cut, so the sign must come from the cut-point measurement -- but a "
                f"real QPU returns only aggregated final counts, losing the per-shot "
                f"cut-point outcome needed to decouple it. Supported on CPU (qiskit, "
                f"qibo, cudaq) and GPU, in both the all-in-one and three-step "
                f"workflows. Use a CPU/GPU run, wire cutting, or a cut whose qubits "
                f"are not measured with identity."
            )


def _decimal_to_base6(num):
    if num == 0:
        return "0"

    base6 = ""
    while num > 0:
        base6 = str(num % 6) + base6
        num //= 6
    return base6


#: Measure-and-prepare gate cut: what each side (1: first qubit, 2: second) of a cut
#: CZ becomes in each of the 6 terms. "rz" is Rz(-pi/2); "z" Z then Rz(-pi/2);
#: "measure" a mid-circuit measurement then Rz(-pi/2); "rzpi" Rz(-pi); "remove"
#: nothing. Terms 3 and 5 carry a minus sign. Side 1 is the same for terms 2 and 3,
#: side 2 for terms 4 and 5.
_GATE_SIDE_OPS = {
    1: ("rz", "z", "measure", "measure", "rzpi", "remove"),
    2: ("rz", "z", "rzpi", "remove", "measure", "measure"),
}

# The measure-and-prepare cut of RZZ(phi) = exp(-i phi/2 Z(x)Z) (Mitarai & Fujii):
# with theta = -phi/2, c = cos(theta), s = sin(theta), the gate's channel is
#   c**2 [I (x) I] + s**2 [Z (x) Z] + c s ([M (x) Rz+] - [M (x) Rz-]
#                                         + [Rz+ (x) M] - [Rz- (x) M]),
# M the signed Z measurement of a cut point (outcome 0 -> -1, as the engine reads
# it), Rz+- = RZ(+-pi/2). gamma = 1 + 2|sin(phi)|: 3 for a CZ-like angle, -> 1 for
# a weak one. Per side, the operation of each of the six terms:
_RZZ_SIDE_OPS = {
    1: ("remove", "pz", "pmeasure", "pmeasure", "rz+", "rz-"),
    2: ("remove", "pz", "rz+", "rz-", "pmeasure", "pmeasure"),
}


def _rzz_cut_weights(phi):
    """The six term weights of the measure-and-prepare cut of ``RZZ(phi)``."""
    theta = -float(phi) / 2
    c, s = math.cos(theta), math.sin(theta)
    return np.array([c * c, s * s, c * s, -c * s, c * s, -c * s])


def _u_angles(matrix):
    """``(theta, phi, lam)`` of OpenQASM's ``u`` equal to a 2x2 unitary up to a
    global phase."""
    m = np.asarray(matrix, dtype=complex)
    theta = 2 * math.atan2(abs(m[1, 0]), abs(m[0, 0]))
    if abs(m[0, 0]) > 1e-9:
        phase = cmath.phase(m[0, 0])
        lam = cmath.phase(-m[0, 1]) - phase if abs(m[0, 1]) > 1e-9 else 0.0
        if abs(m[1, 0]) > 1e-9:
            phi = cmath.phase(m[1, 0]) - phase
        else:  # theta = 0: only phi + lam is defined
            phi = cmath.phase(m[1, 1]) - phase - lam
    else:  # theta = pi: only phi - lam is defined
        lam = 0.0
        phase = cmath.phase(-m[0, 1])
        phi = cmath.phase(m[1, 0]) - phase
    return theta, phi, lam


_H2 = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
_S2 = np.diag([1, 1j])
#: V with V Z V^dag = P: exp(i t P(x)P) = (V(x)V) exp(i t Z(x)Z) (V(x)V)^dag.
_TO_PAULI = {"X": _H2, "Y": _S2 @ _H2, "Z": np.eye(2, dtype=complex)}


def _sequence_unitary(ops):
    """The 4x4 unitary of ``[(name, local qubits, params)]`` on two qubits, the
    first most significant."""
    from Qdislib.utils.statevector import gate_unitary

    total = np.eye(4, dtype=complex)
    for name, qubits, params in ops:
        m = gate_unitary(name, params)
        if len(qubits) == 1:
            m = np.kron(m, np.eye(2)) if qubits[0] == 0 else np.kron(np.eye(2), m)
        elif tuple(qubits) == (1, 0):
            swap = np.eye(4)[[0, 2, 1, 3]]
            m = swap @ m @ swap
        total = m @ total
    return total


def _same_up_to_phase(a, b):
    return abs(np.trace(np.asarray(a).conj().T @ np.asarray(b))) / 4 > 1 - 1e-9


def _rotation_decomposition(matrix):
    """A two-qubit unitary as single-qubit gates and at most three ``rzz``
    interactions, ``[(name, local qubits, params)]`` in time order, so a
    measure-and-prepare cut can cut each ``rzz`` at its own angle instead of
    paying ``gamma = 3`` per CZ of a CZ decomposition.

    A diagonal gate (``cp``, ``crz``, ``rzz``, ...) is ``e^{i a}(RZ (x) RZ) RZZ``
    in closed form; any other goes through the KAK (Weyl) decomposition, its
    ``XX`` / ``YY`` / ``ZZ`` parts each an ``rzz`` between basis changes.
    Checked against ``matrix`` (up to a global phase).

    :raises RuntimeError: if no decomposition rebuilds the gate.
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _diag_angles

    matrix = np.asarray(matrix, dtype=complex)
    if np.allclose(matrix, np.diag(np.diag(matrix)), atol=1e-12):
        _, beta, gamma, delta = _diag_angles(np.diag(matrix))
        ops = [
            ("rz", (0,), [-2 * gamma]),
            ("rz", (1,), [-2 * beta]),
            ("rzz", (0, 1), [-2 * delta]),
        ]
        ops = [op for op in ops if op[0] == "rzz" or abs(op[2][0]) > 1e-12]
        if _same_up_to_phase(_sequence_unitary(ops), matrix):
            return ops
    from qiskit.synthesis import TwoQubitWeylDecomposition

    weyl = TwoQubitWeylDecomposition(matrix)
    for k1a, k1b, k2a, k2b in (
        (weyl.K1l, weyl.K1r, weyl.K2l, weyl.K2r),
        (weyl.K1r, weyl.K1l, weyl.K2r, weyl.K2l),
    ):
        # time order: K2, then exp(i a XX), exp(i b YY), exp(i c ZZ), then K1
        pending = {0: np.asarray(k2a), 1: np.asarray(k2b)}
        ops = []
        for coefficient, pauli in ((weyl.a, "X"), (weyl.b, "Y"), (weyl.c, "Z")):
            if abs(coefficient) < 1e-12:
                continue
            v = _TO_PAULI[pauli]
            for q in (0, 1):
                ops.append(("u", (q,), list(_u_angles(v.conj().T @ pending[q]))))
                pending[q] = v
            ops.append(("rzz", (0, 1), [-2 * coefficient]))
        for q, k1 in ((0, k1a), (1, k1b)):
            ops.append(("u", (q,), list(_u_angles(np.asarray(k1) @ pending[q]))))
        if _same_up_to_phase(_sequence_unitary(ops), matrix):
            return ops
    raise RuntimeError("could not decompose the two-qubit gate into rzz rotations")


def _gate_matrix_4x4(name, params=()):
    """A two-qubit gate's matrix (first qubit most significant): the OpenQASM
    table, else qibo's (``fsim``, ``generalizedfsim``, ...); ``None`` if neither
    knows it."""
    from Qdislib.utils.statevector import gate_unitary

    matrix = gate_unitary(name, params)
    if matrix is None:
        try:
            from Qdislib.core.cutting_algorithms.hadamard_cut import _two_qubit_matrix

            matrix = _two_qubit_matrix(name, params)
        except Exception:  # noqa: BLE001
            return None
    return matrix if np.shape(matrix) == (4, 4) else None


def _rotation_gamma(name, params=()):
    """``gamma`` of the angle-aware measure-and-prepare cut of a two-qubit gate:
    ``prod (1 + 2|sin phi|)`` over the ``rzz`` rotations it is cut as."""
    kind = canonical(name, params)[0]
    if kind in ("cz", "cx", "cnot"):
        return 3.0
    matrix = _gate_matrix_4x4(name, params)
    if matrix is None:
        raise ValueError(f"no matrix for the gate {name!r}")
    ops = _rotation_decomposition(matrix)
    return math.prod(1 + 2 * abs(math.sin(op[2][0])) for op in ops if op[0] == "rzz")


def _replace_with_ops(dag, gate_name, ops, mark):
    """Replace a two-qubit gate node with ``ops`` (``[(name, local qubits,
    params)]``, in time order) in place; return the ids of the ops named ``mark``
    (the cut targets)."""
    node = dag.nodes[gate_name]
    qubits = node["qubits"]
    qmap = {0: qubits[0], 1: qubits[1]}
    preds = {
        q: _find_nodes_with_qubit(dag, gate_name, qubit=q, direction="predecessor")
        for q in qubits
    }
    succs = {
        q: _find_nodes_with_qubit(dag, gate_name, qubit=q, direction="successor")
        for q in qubits
    }
    dag.remove_node(gate_name)

    last = {q: (preds[q][0] if preds[q] else None) for q in qubits}
    marked = []
    for i, (name, local_qubits, op_params) in enumerate(ops):
        node_qubits = tuple(qmap[lq] for lq in local_qubits)
        # A symmetric cut gate still has a SUBS1 / SUBS2 side by qubit *position*
        # (and the identity-observable sign recovery reads it), so pin the first
        # side on the lower qubit -- a (1, 0) orientation otherwise mis-reconstructs.
        if name.lower() == mark.lower():
            node_qubits = tuple(sorted(node_qubits))
        node_id = f"{name.upper()}_{gate_name}_{i}"
        dag.add_node(
            node_id,
            gate=canonical(name)[0],
            qubits=node_qubits,
            parameters=list(op_params),
        )
        for q in node_qubits:
            if last[q] is not None:
                dag.add_edge(last[q], node_id, color="blue")
            last[q] = node_id
        if name.lower() == mark.lower():
            marked.append(node_id)
    for q in qubits:
        if succs[q] and last[q] is not None:
            dag.add_edge(last[q], succs[q][0], color="blue")
    return marked


def _expand_rotation_cut(dag, gate_name):
    """Replace a gate other than CZ/CX by single-qubit gates and the ``rzz``
    rotations of :func:`_rotation_decomposition` (mutating ``dag``); return the
    ``rzz`` node ids, each cut on its own at its angle. An ``rzz`` is left as is."""
    node = dag.nodes[gate_name]
    if len(node.get("qubits") or ()) != 2:
        raise InvalidCutError(
            f"Gate cutting needs a two-qubit gate; {gate_name!r} is not one"
        )
    kind = canonical(str(node["gate"]), node.get("parameters") or ())[0]
    if kind == "rzz":
        return [gate_name]
    matrix = _gate_matrix_4x4(str(node["gate"]), node.get("parameters") or ())
    if matrix is None:
        raise InvalidCutError(
            f"Gate cutting needs a two-qubit gate; cannot build {gate_name!r}"
        )
    return _replace_with_ops(dag, gate_name, _rotation_decomposition(matrix), "rzz")


def _apply_gate_side(dag, node, op, after):
    """Turn one side marker of a measure-and-prepare gate cut into its part of a term
    (see :data:`_GATE_SIDE_OPS`), in place. ``after`` names the Rz(-pi/2) node that
    "z" and "measure" add behind the marker."""
    preds, succs = list(dag.pred[node]), list(dag.succ[node])
    data = dag.nodes[node]
    if op == "remove":
        dag.remove_node(node)
        if preds and succs:
            dag.add_edge(preds[0], succs[0], color="blue")
    elif op in ("rz", "rzpi"):
        data["gate"] = "rz"
        data["parameters"] = [-np.pi / 2 if op == "rz" else -np.pi]
    elif op in ("rz+", "rz-"):  # the angle-aware cut's RZ(+-pi/2)
        data["gate"] = "rz"
        data["parameters"] = [np.pi / 2 if op == "rz+" else -np.pi / 2]
    elif op in ("pz", "pmeasure"):  # ... and its plain Z / cut-point measurement
        data["gate"] = op[1:]
        data["parameters"] = ()
    elif op in ("z", "measure"):
        data["gate"] = op
        dag.add_node(
            after, gate="rz", qubits=data.get("qubits"), parameters=[-np.pi / 2]
        )
        if succs:
            dag.remove_edge(node, succs[0])
            dag.add_edge(after, succs[0], color="blue")
        dag.add_edge(node, after, color="blue")
    else:
        raise ValueError(f"Unexpected gate-cut side operation: {op}")


# Plain function (not a @task): the gate-cut substitution is cheap graph surgery,
# run on the master, populating graph_components in place with concrete DiGraphs.
# Keeping it task-free (no COLLECTION_OUT futures) lets the generated subcircuits
# be wrapped as serializable Subcircuit objects under PyCOMPSs.
def _apply_gate_substitution(updated_dag, gates_cut, index):
    """Apply the base-6 measure-and-prepare substitution for gate-cut configuration
    ``index`` to the SUBS1/SUBS2 marker nodes (mutates ``updated_dag`` in place, no
    component split). Factored out of :func:`_generate_gate_cutting` so the mixed
    wire+gate driver can apply the gate substitution and the wire substitution to the
    same DAG before a single split."""
    digits = list(map(int, _decimal_to_base6(index).zfill(len(gates_cut))))
    for k, digit in enumerate(digits, start=1):
        if digit > 5:
            raise ValueError(f"Unexpected gate-cut substitution index: {digit}")
        for side in (1, 2):
            _apply_gate_side(
                updated_dag,
                f"SUBS{side}_{k}",
                _GATE_SIDE_OPS[side][digit],
                f"SUBS{side}{side}_{k}",
            )


def _generate_gate_cutting(updated_dag, gates_cut, index, graph_components):
    _apply_gate_substitution(updated_dag, gates_cut, index)
    for i, c in enumerate(
        nx.connected_components(updated_dag.to_undirected(as_view=True))
    ):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")

    return updated_dag


# @constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign_gate_cutting(expectation_value, index):
    """Signed product of subcircuit values for one gate-cut term (base 6, {3, 5}).

    See :func:`Qdislib.core.cutting_algorithms._reconstruction.signed_product`.
    """
    # Resolved and reduced on the worker (not the master), so the subcircuit
    # futures feed straight into this task -- no barrier before the sign step.
    from Qdislib.core.cutting_algorithms._reconstruction import signed_product

    return signed_product(expectation_value, index, 6, {3, 5})


def _measured_cut_qubits(dag):
    """Renumbered qubits that a gate-cut term *measures* (``SUBS`` -> ``measure``).

    Their measure-and-prepare sign must be recovered from the cut-point
    measurement when the user observes them with identity -- see
    :func:`_expectation_from_counts_qiskit`. Returns an empty set for a ``None``
    component (an empty gate-cut term, tolerated downstream).
    """
    return frozenset(_cut_measurements(dag))


def _cut_measurements(dag):
    """Like :func:`_measured_cut_qubits`, but one entry per cut measurement, sorted.

    A qubit cut twice appears twice, so each cut's sign is applied (see
    :func:`_cut_point_rounds`).
    """
    if dag is None:
        return []
    return sorted(
        q
        for q in (
            (dag.nodes[n].get("qubits") or (None,))[0]
            for n in dag
            if dag.nodes[n].get("gate") == "measure"
        )
        if q is not None
    )


def _dag_num_qubits(dag):
    """Number of qubits spanned by a DAG (max qubit index + 1)."""
    qubits = set()
    for node in dag:
        qubits.update(dag.nodes[node].get("qubits") or ())
    return (max(qubits) + 1) if qubits else 0


def _component_original_qubits(component):
    """Sorted original qubit indices a component acts on (pre-renumbering)."""
    qubits = set()
    for node in component:
        qubits.update(component.nodes[node].get("qubits") or ())
    return sorted(qubits)


def _joint_counts_qibo(
    dag, num_qubits, prep_sorted, shots, seed=None, platform="numpy"
):
    """Per-component joint counts ``{(final_bits, cut-point_bits): count}`` for Qibo.

    ``final_bits`` are qubit-0-first over the component's qubits; ``cut-point_bits``
    are the mid-circuit collapse-measurement outcomes for ``prep_sorted`` (the
    measured cut qubits), captured per shot so they stay correlated with the final
    measurement (see :func:`_expectation_from_samples_qibo`). Pass ``seed`` for a
    reproducible sample (Qibo's numpy backend draws from the global numpy RNG).
    ``platform`` is the Qibo backend that runs it.
    """
    import qibo
    from qibo import gates as qgates
    from Qdislib.utils.graph_qibo import dag_to_circuit_qibo

    built = dag_to_circuit_qibo(dag, num_qubits)
    if built is None:
        return None
    qc = built[0] if isinstance(built, (list, tuple)) else built
    qibo.set_backend(platform)
    if seed is not None:
        np.random.seed(seed)
    collapse = {}  # every collapse gate of each qubit (a qubit can be cut twice)
    for gate in qc.queue:
        if isinstance(gate, qgates.M) and getattr(gate, "collapse", False):
            collapse.setdefault(gate.qubits[0], []).append(gate)
    qc.add(qgates.M(*range(qc.nqubits)))
    final = np.asarray(qc(nshots=shots).samples())  # (shots, nqubits), qubit-0-first
    cut = {
        q: (
            _fold_cut_outcomes(
                [np.asarray(c.result.samples()).ravel() for c in collapse[q]]
            )
            if q in collapse
            else np.zeros(len(final), dtype=int)
        )
        for q in prep_sorted
    }
    joint = {}
    for shot in range(len(final)):
        ft = tuple(int(final[shot, q]) for q in range(qc.nqubits))
        ct = tuple(int(cut[q][shot]) for q in prep_sorted)
        joint[(ft, ct)] = joint.get((ft, ct), 0) + 1
    return joint


def _joint_counts_qiskit(dag, num_qubits, prep_sorted, shots, seed=None):
    """Per-component joint counts for Qiskit, reusing the cut-measure register ``c``.

    :func:`_counts_qiskit` returns keys ``"<measure-all> <cut-measure c>"`` (each
    little-endian), so the cut-point outcome of qubit ``q`` lives in the second
    register at position ``num_qubits-1-q`` -- exactly the register the expectation
    path reads for the decoupled sign. ``seed`` seeds the Aer sampler.
    """
    from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

    built = dag_to_circuit_qiskit(dag, num_qubits)
    if built is None:
        return None
    qc = built[0] if isinstance(built, (list, tuple)) else built
    return _joint_from_qiskit_counts(
        _counts_qiskit(qc, shots=shots, seed=seed), qc, num_qubits, prep_sorted
    )


def _joint_from_qiskit_counts(counts, qc, num_qubits, prep_sorted):
    """``{(final_bits, cut_point_bits): count}`` from the Qiskit-style counts of
    ``qc`` run with ``measure_all`` appended (keys ``"<meas> <c{K}> ... <c>"``,
    each little-endian), as Aer and the IBM / IQM samplers return them. ``qc`` is
    the circuit *without* that final measurement: its measurements are the cut
    points."""
    rounds = _cut_point_rounds(
        qc.find_bit(inst.qubits[0]).index
        for inst in qc.data
        if inst.operation.name == "measure"
    )
    joint = {}
    for key, cnt in counts.items():
        measure_all = key.split()[0]
        ft = tuple(int(measure_all[num_qubits - 1 - q]) for q in range(num_qubits))
        # One bit per qubit standing for all its cut-point outcomes: the sign
        # -(-1)**bit the correlators apply is -1 for bit 0, +1 for bit 1.
        ct = tuple(
            (
                0
                if _cut_point_sign_qiskit(key, q, rounds.get(q, 1), num_qubits) < 0
                else 1
            )
            for q in prep_sorted
        )
        joint[(ft, ct)] = joint.get((ft, ct), 0) + cnt
    return joint


def _joint_counts_cudaq(dag, num_qubits, prep_sorted, shots, seed=None, target=None):
    """Per-component joint counts for CUDA-Q via ``explicit_measurements``.

    With a measured cut qubit the kernel is built with ``explicit_final=True`` so
    each per-shot key is ``[cut-point outcomes in mz order] + [final q0..q(n-1)]``
    (see :func:`_expectation_cudaq_decoupled`); otherwise plain counts suffice.
    Pass ``seed`` for a reproducible sample; ``target`` is the CUDA-Q target (see
    :func:`~Qdislib.core.backend.cudaq_backend.select_cudaq_target`).
    """
    import cudaq
    from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel
    from Qdislib.core.cutting_algorithms.wire_cutting import _counts_cudaq

    from Qdislib.core.backend.cudaq_backend import select_cudaq_target

    select_cudaq_target(target)
    if seed is not None:
        cudaq.set_random_seed(seed)

    if prep_sorted:
        # cudaq.run (not sample) so the mid-circuit cut measurements collapse; each
        # per-shot record is [cut-point outcomes in mz order] + [final q0..q(n-1)].
        kernel, _obs, measured_order = dag_to_circuit_cudaq_kernel(
            dag, num_qubits, explicit_final=True
        )
        positions = {}  # every cut-point outcome of each qubit
        for i, qb in enumerate(measured_order):
            positions.setdefault(qb, []).append(i)
        n_mid = len(measured_order)
        records = cudaq.run(kernel, shots_count=shots)
        joint = {}
        for record in records:
            # record is an int; bit k = outcome k of [cut-points..., finals...].
            ft = tuple((record >> (n_mid + q)) & 1 for q in range(num_qubits))
            ct = tuple(
                (
                    _fold_cut_outcomes([(record >> p) & 1 for p in positions[q]])
                    if q in positions
                    else 0
                )
                for q in prep_sorted
            )
            joint[(ft, ct)] = joint.get((ft, ct), 0) + 1
        return joint

    built = dag_to_circuit_cudaq_kernel(dag, num_qubits)
    kernel = built[0] if isinstance(built, (list, tuple)) else built
    freq = _counts_cudaq(kernel, shots, target)
    joint = {}
    for key, cnt in freq.items():
        ft = tuple(int(key[q]) for q in range(num_qubits))
        joint[(ft, ())] = joint.get((ft, ()), 0) + cnt
    return joint


def _joint_counts_pennylane(
    dag, num_qubits, prep_sorted, shots, seed=None, device="default.qubit"
):
    """Per-component joint counts for PennyLane via mid-circuit measurements.

    Runs one QNode that samples the final measurements and each measured cut
    qubit's mid-circuit ``MeasurementValue``, kept correlated per shot (see
    :func:`~Qdislib.core.cutting_algorithms.wire_cutting._expectation_from_samples_pennylane`);
    reduces to plain final counts when there is no measured cut qubit. Both
    ``final_bits`` and ``cut-point_bits`` are qubit-0-first. Pass ``seed`` for a
    reproducible sample.
    """
    from Qdislib.utils.graph_pennylane import dag_to_circuit_pennylane

    built = dag_to_circuit_pennylane(dag, num_qubits)
    if built is None:
        return None
    circuit = built[0] if isinstance(built, (list, tuple)) else built
    ops, width = circuit.ops, circuit.num_qubits

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
        return [qml.sample(wires=range(width))] + [qml.sample(value) for value in mids]

    qnode = _build_pennylane_qnode(qfunc, width, shots, device, seed)
    out = qnode()
    final = np.asarray(out[0]).reshape(shots, width)
    mid_samples = _folded_mid_samples(out[1:], collapse_order)

    joint = {}
    for shot in range(shots):
        ft = tuple(int(final[shot, q]) for q in range(width))
        ct = tuple(
            int(mid_samples[q][shot]) if q in mid_samples else 0 for q in prep_sorted
        )
        joint[(ft, ct)] = joint.get((ft, ct), 0) + 1
    return joint


def _z_correlators(joint, prep_sorted, k, masks=None):
    """The Z-string correlators ``<Z^b>`` of one component's samples: all ``2**k``
    (a dict), or only ``masks``.

    ``corr[mask]`` (bit ``i`` of ``mask`` <=> local qubit ``i`` in the Z-string) is
    computed with the gate-cut prep-qubit sign rule: a measured cut qubit observed
    with **Z** uses its final bit ``(-1)**f``, but observed with **identity** uses
    the *cut-point* outcome ``-(-1)**m``. This is what lets a cut qubit that keeps
    evolving after the cut reconstruct correctly (its final bit no longer equals its
    cut-point outcome). Reduces to the plain Z/I correlators when there is no prep.

    That sign, ``prod -(-1)**m``, is the same for every mask, so the correlators are
    the Walsh-Hadamard transform of the signed histogram of final outcomes: vectorised,
    and for ``masks`` a sum over the outcomes only (no ``2**k`` vector).
    """
    if joint is None or k == 0:
        return {0: 1.0}
    total = sum(joint.values()) or 1
    prep_index = [i for i, _ in enumerate(prep_sorted)]
    items = list(joint.items())
    finals = np.array([[int(b) for b in final[:k]] for (final, _), _ in items])
    outcome = (finals << np.arange(k)).sum(axis=1)
    signed = np.array([cnt / total for _, cnt in items])
    if prep_index:
        cuts = np.array([[int(cut[i]) for i in prep_index] for (_, cut), _ in items])
        signed = signed * np.prod(-((-1.0) ** cuts), axis=1)
    if masks is not None:
        masks = [int(m) for m in masks]
        union = 0
        for mask in masks:
            union |= mask
        bits = [q for q in range(k) if (union >> q) & 1]
        if len(masks) > 2 * len(bits) and len(bits) <= 22:
            # many masks (a marginal distribution): one histogram over the bits
            # they use, one transform, instead of a pass per mask
            packed = np.zeros(len(outcome), dtype=np.int64)
            for j, q in enumerate(bits):
                packed |= ((outcome >> q) & 1) << j
            histogram = np.zeros(2 ** len(bits))
            np.add.at(histogram, packed, signed)
            full = _fast_walsh_hadamard(histogram)

            def pack(mask):
                return sum(((mask >> q) & 1) << j for j, q in enumerate(bits))

            return {mask: float(full[pack(mask)]) for mask in masks}
        corr = {}
        for mask in masks:
            parity = np.zeros(len(outcome), dtype=np.int64)
            for q in range(k):
                if (mask >> q) & 1:
                    parity ^= (outcome >> q) & 1
            corr[int(mask)] = float(np.sum(signed * (1 - 2 * parity)))
        return corr
    histogram = np.zeros(2**k)
    np.add.at(histogram, outcome, signed)
    return dict(enumerate(_fast_walsh_hadamard(histogram).tolist()))


class _JointCounts:
    """A fragment run's joint final / cut-point counts: the result the cache
    keeps, which serves any observable. A plain object, not a container, so that
    PyCOMPSs ships it whole (a dict inside a ``COLLECTION_IN`` is rebuilt as a
    plain dict); the reconstruction reads it through :class:`_CountCorrelators`.

    It also holds its outcomes as arrays, built once: a bootstrap resample
    (:meth:`resampled`) is a new count vector over the same outcomes.
    """

    def __init__(self, joint, prep_sorted, k, _arrays=None):
        self._joint = joint
        self.prep_sorted = list(prep_sorted)
        self.k = int(k)
        self._arrays = _arrays

    @property
    def joint(self):
        """``{(final_bits, cut_point_bits): count}``."""
        if self._joint is None:
            outcomes, _, _, counts = self._arrays
            self._joint = {o: int(c) for o, c in zip(outcomes, counts) if c}
        return self._joint

    def arrays(self):
        """``(outcomes, bits, sign, counts)``: the outcomes, their final bits (one
        row each), their cut-point sign ``prod -(-1)**m`` and their counts."""
        if self._arrays is None:
            items = list(self._joint.items())
            outcomes = [outcome for outcome, _ in items]
            bits = np.array(
                [[int(b) for b in final[: self.k]] for final, _ in outcomes],
                dtype=np.uint8,
            ).reshape(len(items), self.k)
            sign = np.ones(len(items))
            if self.prep_sorted:
                cuts = np.array(
                    [[int(b) for b in cut] for _, cut in outcomes], dtype=float
                ).reshape(len(items), -1)
                sign = np.prod(-((-1.0) ** cuts), axis=1)
            counts = np.array([count for _, count in items], dtype=float)
            self._arrays = (outcomes, bits, sign, counts)
        return self._arrays

    def resampled(self, rng):
        """The same number of shots drawn again from these counts (bootstrap)."""
        outcomes, bits, sign, counts = self.arrays()
        total = counts.sum()
        if total <= 0:
            return self
        drawn = rng.multinomial(int(round(total)), counts / total).astype(float)
        return _JointCounts(
            None, self.prep_sorted, self.k, (outcomes, bits, sign, drawn)
        )

    def correlators(self):
        return _CountCorrelators(counts=self)

    def __getstate__(self):  # ship the counts, not the arrays
        return {
            "_joint": self.joint,
            "prep_sorted": self.prep_sorted,
            "k": self.k,
            "_arrays": None,
        }


class _CountCorrelators(dict):
    """The correlators of :class:`_JointCounts`, computed when first read (any
    mask)."""

    #: Up to this many qubits every correlator is computed at once (a
    #: Walsh-Hadamard transform); above, only the ones read.
    FULL_TRANSFORM_QUBITS = 10

    def __init__(self, joint=None, prep_sorted=(), k=0, counts=None):
        super().__init__()
        self.counts = (
            counts if counts is not None else _JointCounts(joint, prep_sorted, k)
        )

    @property
    def joint(self):
        return self.counts.joint

    @property
    def prep_sorted(self):
        return self.counts.prep_sorted

    @property
    def k(self):
        return self.counts.k

    def __missing__(self, mask):
        mask = int(mask)
        if self.k <= self.FULL_TRANSFORM_QUBITS and 0 <= mask < 2**self.k:
            value = float(self._transform()[mask])
            dict.__setitem__(self, mask, value)
            return value
        _, bits, sign, counts = self.counts.arrays()
        total = counts.sum() or 1.0
        picked = [q for q in range(self.k) if (mask >> q) & 1]
        parity = (
            bits[:, picked].sum(axis=1, dtype=np.int64) & 1 if picked else 0
        )  # signed: 1 - 2 * parity must not wrap
        value = float(np.sum(counts * sign * (1 - 2 * parity)) / total)
        dict.__setitem__(self, mask, value)
        return value

    def values_at(self, masks):
        """The correlators at ``masks`` (an array); when they are most of them,
        all ``2**k`` at once by a Walsh-Hadamard transform (a distribution)."""
        masks = np.asarray(masks, dtype=np.int64)
        if self.k <= 22 and len(masks) * 16 >= 2**self.k:
            return self._transform()[masks]
        return np.array([self[int(m)] for m in masks])

    def _transform(self):
        """Every correlator, by a Walsh-Hadamard transform of the signed histogram
        (once): from the count arrays, so a bootstrap resample builds no dict."""
        if getattr(self, "_full", None) is None:
            _, bits, sign, counts = self.counts.arrays()
            total = counts.sum() or 1.0
            outcome = bits.astype(np.int64) @ (1 << np.arange(self.k, dtype=np.int64))
            histogram = np.zeros(2**self.k)
            np.add.at(histogram, outcome, counts * sign / total)
            self._full = _fast_walsh_hadamard(histogram)
        return self._full

    def __reduce__(self):
        return (_CountCorrelators, (self.joint, self.prep_sorted, self.k))


def _prepare_sampled_component(component):
    """Renumber a fragment DAG to a standalone ``0..n`` subcircuit and resolve its
    measure-and-prepare (cut-point) qubits: ``(dag, num_qubits, prep_sorted,
    original)``, ``dag`` ``None`` for an empty component."""
    original = _component_original_qubits(component)
    # _update_qubits renumbers in place: work on a copy, so a DAG the caller keeps
    # (a subcircuit's, read again by the next result_from_counts) is left as it was
    s_new, highest_qubit = _update_qubits(component.copy())
    if s_new is None:  # empty component (a qubit touched only by a removed cut node)
        return (None, 0, [], original)
    prep_sorted = sorted(_measured_cut_qubits(s_new))
    return (s_new, highest_qubit, prep_sorted, original)


def _sample_prepared_component(
    prepared,
    software,
    shots,
    qpu=False,
    qpu_dict=None,
    seed=None,
    masks=None,
    keep_counts=False,
):
    """Build and sample a prepared subcircuit: ``(corr, original)``, ``corr`` its
    ``<Z^b>`` correlators with the cut-point sign baked in (only ``masks``, if
    given; with ``keep_counts`` a :class:`_CountCorrelators` over the counts)."""
    s_new, highest_qubit, prep_sorted, original = prepared
    if s_new is None:  # empty component
        return ({0: 1.0}, [])
    from Qdislib.core.backend.registry import build_execution_backend

    backend = build_execution_backend(
        software, num_qubits=highest_qubit, qpu=bool(qpu), qpu_dict=qpu_dict
    )
    from Qdislib.core.cutting_algorithms._api import is_simulator

    deferred = None
    if prep_sorted and not qpu and is_simulator(software):
        deferred = _defer_cut_measurements(s_new, highest_qubit, prep_sorted)
    if deferred is not None:
        joint = _deferred_joint_counts(
            backend, deferred, highest_qubit, len(prep_sorted), shots, seed
        )
    else:
        joint = backend.evaluate_joint_counts(
            s_new, highest_qubit, prep_sorted, shots, seed=seed
        )
    if keep_counts:
        return (_JointCounts(joint, prep_sorted, highest_qubit), original)
    corr = _z_correlators(joint, prep_sorted, highest_qubit, masks)
    return (corr, original)


def _defer_cut_measurements(dag, num_qubits, prep_sorted):
    """The DAG with each cut-point measurement deferred: a CNOT onto a fresh qubit
    (``num_qubits + j`` for the ``j``-th measured qubit), read at the end.

    A mid-circuit measurement whose outcome is only recorded is the same as that
    CNOT (deferred measurement), and without it a simulator samples the final state
    once instead of simulating every shot. ``None`` if a qubit is measured twice.
    """
    measures = [
        node
        for node in dag.nodes
        if (dag.nodes[node].get("gate") or "").lower() == "measure"
    ]
    qubits = [dag.nodes[node]["qubits"][0] for node in measures]
    if sorted(qubits) != list(prep_sorted):
        return None
    deferred = dag.copy()
    for node, qubit in zip(measures, qubits):
        fresh = num_qubits + prep_sorted.index(qubit)
        deferred.nodes[node].update(gate="cx", qubits=(qubit, fresh), parameters=())
    return deferred


def _defer_all_cut_measurements(dag, num_qubits):
    """The DAG with every cut-point measurement deferred onto a fresh qubit of its
    own, a qubit cut twice included: ``(dag, cut_qubits)``, fresh qubit
    ``num_qubits + j`` holding a measurement of ``cut_qubits[j]`` (by qubit, then in
    time order: the order of :func:`_cut_measurements`)."""
    import networkx

    measures = [
        node
        for node in networkx.topological_sort(dag)
        if (dag.nodes[node].get("gate") or "").lower() == "measure"
    ]
    if not measures:
        return dag, []
    qubits = [dag.nodes[node]["qubits"][0] for node in measures]
    ranked = sorted(range(len(measures)), key=lambda i: (qubits[i], i))
    deferred = dag.copy()
    for j, i in enumerate(ranked):
        deferred.nodes[measures[i]].update(
            gate="cx", qubits=(qubits[i], num_qubits + j), parameters=()
        )
    return deferred, [qubits[i] for i in ranked]


def _fold_cut_bits(bits, cut_qubits, prep_sorted):
    """One bit per cut qubit from one bit per cut measurement: the bit whose sign
    ``-(-1)**bit`` is the product of its measurements' (as the engine keeps them)."""
    folded = []
    for q in prep_sorted:
        sign = 1
        for b, owner in zip(bits, cut_qubits):
            if owner == q:
                sign *= -((-1) ** int(b))
        folded.append(0 if sign < 0 else 1)
    return tuple(folded)


def _deferred_joint_counts(backend, dag, num_qubits, num_cuts, shots, seed):
    """``{(final_bits, cut_point_bits): count}`` from one end-of-circuit sample of
    the deferred DAG (the cut-point bits are the fresh qubits)."""
    counts = backend.evaluate_samples(dag, num_qubits + num_cuts, shots, seed=seed)
    joint = {}
    for bits, count in counts.items():
        key = (
            tuple(int(b) for b in bits[:num_qubits]),
            tuple(int(b) for b in bits[num_qubits:]),
        )
        joint[key] = joint.get(key, 0) + count
    return joint


@task(returns=1)
def _create_sample_subcircuit(component):
    """Prepare one gate-cut sampling subcircuit for execution (COMPSs task).

    The *creation* half of sampling a subcircuit, split from execution so the two
    are independently schedulable/routable (see :func:`_prepare_sampled_component`).
    The native circuit is built and sampled later, in
    :func:`_execute_sample_subcircuit` -- a Qibo collapse circuit / CUDA-Q kernel
    cannot cross a task boundary, so it never leaves the execution task.
    """
    return _prepare_sampled_component(component)


@task(returns=1)
def _execute_sample_subcircuit(
    prepared, software, shots, qpu=False, qpu_dict=None, seed=None
):
    """Build + sample a prepared gate-cut subcircuit; return its Z-correlators.

    The *execution* half (COMPSs task) -- this is where the native circuit is built
    and sampled, so it is the task to route to control where subcircuits actually
    run. ``build_execution_backend`` turns ``software``/``qpu``/``qpu_dict`` into a
    concrete backend (local SDK, GPU Aer, or a QPU) with no per-backend branching,
    and its :meth:`evaluate_joint_counts` captures the final and mid-circuit
    cut-point measurements. Returns ``(corr, original)``: ``corr`` maps each Z-string
    mask to its ``<Z^b>`` correlator (``2**|component|`` floats) and ``original``
    lists the component's global qubit indices. Returning correlators (not raw
    counts) is what keeps the reconstruction correct for cut qubits that keep
    evolving after the cut (the cut-point sign is baked in here).
    """
    return _sample_prepared_component(prepared, software, shots, qpu, qpu_dict, seed)


def _sample_gate_subcircuit(
    component, software, shots, qpu=False, qpu_dict=None, seed=None
):
    """Sample one connected-component subcircuit into its Z-correlators.

    Orchestrates two COMPSs tasks -- :func:`_create_sample_subcircuit` (prepare the
    subcircuit) then :func:`_execute_sample_subcircuit` (build + sample it) -- so
    subcircuit creation and execution are separate, independently routable tasks
    (the prepared-spec future feeds create -> execute directly, task -> task). Called
    from the master, so it submits both tasks without nesting. ``qpu``/``qpu_dict``
    route execution to a QPU (via ``build_execution_backend``). Returns the execution
    future ``(corr, original)``.
    """
    prepared = _create_sample_subcircuit(component)
    return _execute_sample_subcircuit(prepared, software, shots, qpu, qpu_dict, seed)


def _gate_cut_sampling_components(dag, gates_cut):
    """Build the connected-component subcircuits for every 6**n_cuts gate-cut config.

    Returns ``(components, number_components)`` where ``components`` is a *flat*,
    configuration-major list of component DAGs (config 0's components, then config
    1's, ...) and ``number_components`` is how many components each configuration
    splits into (uniform -- the cut topology is the same across configurations).
    Sampling is kept separate (:func:`_sample_gate_subcircuit`) so the all-in-one
    :func:`gate_cutting_sampling` and the three-step
    :func:`gate_cutting_sampling_subcircuits` share one generation path.
    """
    new_dag = _generate_cut(dag, gates_cut)
    components = []
    number_components = 1
    for index in range(6 ** len(gates_cut)):
        copy_graph = _remove_red_edges(new_dag)
        number_components = nx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
        )
        graph_components = [nx.DiGraph() for _ in range(number_components)]
        _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)
        components.extend(graph_components)
    return components, number_components


def _fast_walsh_hadamard(vector):
    """In-place fast Walsh-Hadamard transform: ``out[x] = sum_A (-1)**popcount(A&x) v[A]``."""
    a = np.asarray(vector, dtype=float).copy()
    n = len(a)
    h = 1
    while h < n:
        a = a.reshape(-1, 2 * h)
        left = a[:, :h].copy()
        right = a[:, h:].copy()
        a[:, :h] = left + right
        a[:, h:] = left - right
        a = a.reshape(-1)
        h *= 2
    return a


def _reconstruct_gate_cut_distribution(configs, num_cuts, num_qubits):
    """Reconstruct the output distribution from sampled gate-cut configurations.

    The correct distribution is recovered from Pauli-Z correlators, not from
    tensoring the components' *final* measurement distributions: a cut qubit that
    keeps evolving after the cut has a forced final measurement, so its raw final
    distribution is wrong -- but its ``<Z^b>`` correlator (built with the cut-point
    sign in :func:`_z_correlators`) is right. Each component supplies ``<Z^b>`` over
    its qubits; per configuration these multiply (signed by ``sign/2**n_cuts``) into
    the global correlator ``<Z^A>``, and ``P(x) = 2**-n * sum_A (-1)**(A.x) <Z^A>``
    (an inverse Walsh-Hadamard transform). Qubits covered by no component are
    identity there, i.e. ``|0>``. The signed result is abs-normalised into a valid
    distribution. Validated against the exact statevector distribution (TVD ~ shot
    noise) for terminal *and* non-terminal cuts.
    """
    dim = 2**num_qubits
    z_correlators = np.zeros(dim)
    for index, per_component in enumerate(configs):
        # Sign from the base-6 digits (3 or 5 flip), mirroring _change_sign_gate_cutting.
        sign, digits = 1, index
        while digits:
            if digits % 6 in (3, 5):
                sign = -sign
            digits //= 6

        contribution = np.full(dim, sign / (2**num_cuts), dtype=float)
        for corr, original in per_component:
            # Map each global Z-string A to this component's local mask, then read
            # off its correlator. Uncovered qubits fall through as identity (mask 0).
            component_values = np.empty(dim)
            for global_mask in range(dim):
                local_mask = 0
                for local_index, qubit in enumerate(original):
                    if (global_mask >> qubit) & 1:
                        local_mask |= 1 << local_index
                component_values[global_mask] = corr[local_mask]
            contribution *= component_values
        z_correlators += contribution

    probabilities = _fast_walsh_hadamard(z_correlators) / dim
    keys = ["".join(str((x >> q) & 1) for q in range(num_qubits)) for x in range(dim)]
    norm = float(np.abs(probabilities).sum()) or 1.0
    return {key: abs(float(p)) / norm for key, p in zip(keys, probabilities)}


@task(returns=1, sampled=COLLECTION_IN)
def _reconstruct_gate_cut_distribution_task(
    sampled, num_cuts, num_qubits, number_components
):
    """Task wrapper around :func:`_reconstruct_gate_cut_distribution`.

    Takes the flat, configuration-major list of per-component ``(corr, original)``
    samples via ``COLLECTION_IN`` -- so the :func:`_sample_gate_subcircuit` futures
    feed straight into reconstruction (task -> task, no master barrier) -- regroups
    them by configuration, and returns the reconstructed distribution.
    """
    configs = [
        sampled[i : i + number_components]
        for i in range(0, len(sampled), number_components)
    ]
    return _reconstruct_gate_cut_distribution(configs, num_cuts, num_qubits)


def _gate_cutting_sampling(
    circuit,
    gates_cut,
    shots=1024,
    software="qiskit",
    sync=True,
    qpu=False,
    qpu_dict=None,
    seed=None,
):
    """Reconstruct the output **probability distribution** of a gate-cut circuit.

    Unlike :func:`gate_cutting` (which returns an expectation value), this samples
    each subcircuit and reconstructs the full distribution over measurement
    bitstrings. Every subcircuit is sampled independently -- so under PyCOMPSs the
    subcircuits run in parallel -- and combined via the gate-cut quasiprobability
    decomposition, then abs-normalised into a valid distribution.

    :param circuit: A Qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: The gates to cut (as in :func:`gate_cutting`).
    :param shots: Measurement shots per subcircuit. Sampling reconstruction has a
        quasiprobability overhead, so it needs more shots than an expectation.
    :param software: Execution backend, or the name of any registered backend.
    :param sync: Under PyCOMPSs, wait for the sampling tasks and return the
        distribution; set ``False`` to keep the futures unresolved (advanced).
    :param qpu: Route each subcircuit's sampling to a QPU (via
        ``build_execution_backend``) instead of a local simulator.
    :param qpu_dict: QPU-name -> max-qubit map selecting the QPU (e.g.
        ``{"MN_Ona": 4}``), as in :func:`gate_cutting`.
    :param seed: Optional seed for a reproducible distribution. Each subcircuit
        gets a distinct derived seed (``seed + i``) so the pieces stay independent;
        ``None`` keeps the previous non-deterministic sampling. Ignored on a QPU.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    dag = convert_to_dag(circuit)
    num_qubits = _dag_num_qubits(dag)
    _validate_gate_cuts_partition(dag, gates_cut)

    if not gates_cut:
        from Qdislib.core.backend.registry import build_execution_backend

        backend = build_execution_backend(
            software, num_qubits=num_qubits, qpu=bool(qpu), qpu_dict=qpu_dict
        )
        counts = backend.evaluate_samples(dag, num_qubits, shots, seed=seed)
        total = sum(counts.values()) or 1
        return {k: v / total for k, v in counts.items()}

    components, number_components = _gate_cut_sampling_components(dag, gates_cut)
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
    distribution = _reconstruct_gate_cut_distribution_task(
        sampled, len(gates_cut), num_qubits, number_components
    )
    return compss_wait_on(distribution) if sync else distribution


def sample_subcircuit(subcircuit, shots=1024, *, backend=None, seed=None):
    """Run one subcircuit of a ``*_sampling_subcircuits`` list (step two of three)::

        subs = qd.gate_cutting_sampling_subcircuits(circuit, cuts)
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(sampled)

    :param subcircuit: A subcircuit from any ``*_sampling_subcircuits`` list (or a
        ``circuit_cutting`` fragment).
    :param shots: Shots; ``None`` is exact (statevector) where supported.
    :param backend: Where to run it, as in :func:`~Qdislib.expectation_value`.
    :param seed: Seed for a simulator's sampling.
    :raises TypeError: For what is not a subcircuit, or one of a list made for
        expectation values with no distribution to sample (1.x lists).
    :raises ValueError: For ``shots=None`` with a backend that is no simulator.
    :return: ``(corr, original)`` -- its ``<Z^b>`` correlators and its qubits in the
        uncut circuit; ``reconstruct`` takes the list of these.
    """
    from Qdislib.utils.backends import Subcircuit

    if isinstance(subcircuit, Subcircuit):
        # an engine subcircuit keeps its counts: any observable reads back from them
        run = getattr(subcircuit, "_sample_counts", subcircuit._sample)
        return run(shots, backend, seed)
    if hasattr(subcircuit, "sample"):
        return subcircuit.sample(shots, seed)
    raise TypeError(
        f"sample_subcircuit takes a subcircuit of a *_sampling_subcircuits list, "
        f"not {type(subcircuit).__name__}."
    )


def _gate_cutting_sampling_reconstruction(
    sampled_results, number_cuts, number_components, num_qubits
):
    """Recombine per-subcircuit samples into the output probability distribution.

    The sampling counterpart of :func:`gate_cutting_subcircuit_reconstruction`.
    ``sampled_results`` is the flat, configuration-major list of ``(corr, original)``
    from :func:`sample_subcircuit` (futures under PyCOMPSs); it feeds straight into
    the reconstruction task (task -> task, no master barrier).

    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    distribution = _reconstruct_gate_cut_distribution_task(
        sampled_results, number_cuts, num_qubits, number_components
    )
    return compss_wait_on(distribution)


class SampledCutSubcircuits(list):
    """Sampling subcircuits of a gate cut, that also know how to be recombined.

    The sampling analogue of :class:`~Qdislib.utils.backends.CutSubcircuits`:
    behaves like the plain ``list`` of sampling subcircuits (iterate/index/``len``)
    and additionally carries the metadata to reconstruct the output distribution::

        subs = qd.gate_cutting_sampling_subcircuits(circuit, ["CZ_2"])
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(sampled)

    What a 1.x list was; the functions build it, and 2.0 lists are
    :class:`~Qdislib.CutCircuits`.
    """

    def __init__(
        self, subcircuits, *, number_cuts, number_components, num_qubits, reconstructor
    ):
        super().__init__(subcircuits)
        #: Number of gate cuts applied.
        self.number_cuts = number_cuts
        #: Components each cut configuration splits into.
        self.number_components = number_components
        #: Qubit count of the original circuit (distribution width).
        self.num_qubits = num_qubits
        self._reconstructor = reconstructor

    def reconstruct(self, sampled_results):
        """Recombine per-subcircuit samples into the output distribution.

        :param sampled_results: One :func:`sample_subcircuit` result per subcircuit,
            in the same order as the subcircuits in this object.
        :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
        """
        return self._reconstructor(
            sampled_results, self.number_cuts, self.number_components, self.num_qubits
        )


def _gate_cutting_sampling_subcircuits(
    circuit, gates_cut, software="qiskit", qpu=False, qpu_dict=None
):
    """Partition a circuit into **sampling** subcircuits (three-step workflow).

    The sampling counterpart of :func:`gate_cutting_subcircuits`: generate the
    subcircuits, sample each one yourself (any backend/hardware), then reconstruct
    the output *distribution* (rather than an expectation value)::

        import Qdislib as qd

        # 1. Cut into independently samplable subcircuits.
        subs = qd.gate_cutting_sampling_subcircuits(circuit, ["CZ_2"], software="qiskit")

        # 2. Sample each subcircuit (your loop -- real QPUs, batching, caching, ...).
        sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]

        # 3. Recombine into the original circuit's output distribution.
        distribution = subs.reconstruct(sampled)

    :param circuit: A Qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: The gates to cut (as in :func:`gate_cutting_sampling`).
    :param software: Backend the subcircuits are sampled on (must match the backend
        used to sample them).
    :param qpu: Route sampling to a QPU when each subcircuit is later sampled.
    :param qpu_dict: QPU-name -> max-qubit map selecting the QPU (e.g. ``{"MN_Ona": 4}``).
    :return: A :class:`SampledCutSubcircuits` (a ``list`` of sampling subcircuits
        plus the cut metadata and a :meth:`~SampledCutSubcircuits.reconstruct`
        method).
    """
    dag = convert_to_dag(circuit)
    num_qubits = _dag_num_qubits(dag)
    _validate_gate_cuts_partition(dag, gates_cut)
    components, number_components = _gate_cut_sampling_components(dag, gates_cut)
    units = [SampledSubcircuit(component, software) for component in components]
    return SampledCutSubcircuits(
        units,
        number_cuts=len(gates_cut),
        number_components=number_components,
        num_qubits=num_qubits,
        reconstructor=_gate_cutting_sampling_reconstruction,
    )


@task(returns=1)
def _create_gate_subcircuit(component):
    """Prepare one gate-cut subcircuit for execution (COMPSs task).

    The *creation* half, split from execution so the two are independently
    schedulable/routable. Renumbers the component to a standalone ``0..n``
    subcircuit DAG and resolves its measure-and-prepare (prep) qubits. Returns a
    serializable spec ``(dag, num_qubits, prep_qubits)``; the native circuit is
    built later, in :func:`_execute_gate_subcircuit` (a Qibo collapse circuit /
    CUDA-Q kernel cannot cross a task boundary).
    """
    s_new, highest_qubit = _update_qubits(component)
    prep_qubits = _cut_measurements(s_new)
    return (s_new, highest_qubit, prep_qubits)


@task(returns=1)
def _execute_gate_subcircuit(
    prepared,
    max_qubit,
    software,
    shots,
    method,
    gpu,
    gpu_min_qubits,
    qpu,
    qpu_dict,
    batch,
    backend,
    seed=None,
):
    """Build + run a prepared gate-cut subcircuit; return its expectation value.

    The *execution* half (COMPSs task) -- this is where the native circuit is built
    and run, so it is the task to route (e.g. ``@constraint`` or a QPU backend) to
    control where subcircuits execute. ``build_execution_backend`` selects the QPU /
    GPU / software backend (Qiskit Aer / Qibo / CUDA-Q) and evaluates the subcircuit;
    the native circuit never leaves this task.
    """
    s_new, highest_qubit, prep_qubits = prepared

    # No per-backend or per-device branching here: build_execution_backend selects
    # the QPU / GPU / software backend, and it evaluates the subcircuit.
    from Qdislib.core.backend.registry import build_execution_backend

    exec_backend = build_execution_backend(
        software,
        num_qubits=highest_qubit,
        gpu=gpu,
        gpu_min_qubits=gpu_min_qubits,
        qpu=bool(qpu),
        qpu_dict=qpu_dict,
        max_qubit=max_qubit,
        batch=batch,
        ibm_backend=backend,
    )
    return exec_backend.evaluate_expectation(
        s_new, highest_qubit, prep_qubits, shots, method, seed=seed
    )


def _build_and_run_gate_subcircuit(
    component,
    max_qubit,
    software,
    shots,
    method,
    gpu,
    gpu_min_qubits,
    qpu,
    qpu_dict,
    batch,
    backend,
    seed=None,
):
    """Prepare then execute one gate-cut subcircuit as two COMPSs tasks.

    Orchestrates :func:`_create_gate_subcircuit` (prepare) then
    :func:`_execute_gate_subcircuit` (build + run), so subcircuit creation and
    execution are separate, independently routable tasks (the prepared-spec future
    feeds create -> execute directly, task -> task). Called from the master, so it
    submits both tasks without nesting. Returns the execution future (the
    subcircuit's expectation value).
    """
    prepared = _create_gate_subcircuit(component)
    return _execute_gate_subcircuit(
        prepared,
        max_qubit,
        software,
        shots,
        method,
        gpu,
        gpu_min_qubits,
        qpu,
        qpu_dict,
        batch,
        backend,
        seed,
    )


def _execute_gate_cutting(
    dag,
    gates_cut,
    shots=10000,
    method="automatic",
    gpu=False,
    gpu_min_qubits=None,
    qpu=None,
    qpu_dict=None,
    batch=None,
    backend=None,
    software="qiskit",
    seed=None,
):
    new_dag = _generate_cut(dag, gates_cut)
    _assert_identity_observable_supported(
        new_dag, software=software, gpu=bool(gpu), qpu=bool(qpu)
    )

    S = [
        new_dag.subgraph(c).copy()
        for c in nx.connected_components(new_dag.to_undirected(as_view=True))
    ]

    max_qubit_list = [_max_qubits_graph(s) for s in S]

    reconstruction = []
    for index in range(6 ** len(gates_cut)):
        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
        )

        graph_components = [nx.DiGraph() for _ in range(num_components)]

        # Populates graph_components in place (one DiGraph per component).
        _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        exp_values = []
        for i, s in enumerate(graph_components):
            logger.debug("Subcircuit %s", index * 2 + i)
            exp_values.append(
                _build_and_run_gate_subcircuit(
                    s,
                    max_qubit_list[i],
                    software,
                    shots,
                    method,
                    gpu,
                    gpu_min_qubits,
                    qpu,
                    qpu_dict,
                    batch,
                    backend,
                    next_seed(seed, index, i),
                )
            )

        exp = _change_sign_gate_cutting(exp_values, index)
        reconstruction.append(exp)

    return reconstruction


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

    if not isinstance(subcirc, qibo.models.Circuit):
        observables = subcirc[1]
        subcirc = subcirc[0]
    else:
        observables = "Z" * subcirc.nqubits

    from Qdislib.utils.qpu_bridge import submit_and_wait
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expectation_from_counts_qibo,
    )

    result = submit_and_wait(subcirc)
    freq = result.counts()
    # Normalise by the actual shot count returned by the QPU (its num_bins may
    # differ from the classical-side ``shots`` argument).
    total_shots = sum(freq.values()) or 1
    return _expectation_from_counts_qibo(freq, observables, total_shots)


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
def _counts_qiskit(qc, shots=1024, method="automatic", seed=None):
    """Simulate a Qiskit circuit with Aer and return its measurement counts.

    Appends ``measure_all`` to a *copy*, so the caller's circuit is left
    untouched and a subcircuit can be evaluated more than once (e.g. simulate
    first, then submit the same cut to hardware). Then runs it on the CPU Aer
    simulator. Backend/method details are logged at DEBUG level.

    :param qc: A Qiskit ``QuantumCircuit``.
    :param shots: Number of measurement shots.
    :param method: Aer simulation method.
    :param seed: Optional Aer ``seed_simulator`` for a reproducible sample.
    :return: The Aer ``get_counts()`` dict (keys are ``"<measure-all reg> ..."``).
    """
    qc = qc.copy()
    qc.measure_all()

    logger.debug("Simulation using CPUs; method: %s", method)
    simulator = AerSimulator(
        device="CPU",
        method=method,
        max_parallel_threads=1,
        mps_omp_threads=1,
        mps_parallel_threshold=1,
    )

    job = simulator.run(qc, shots=shots, seed_simulator=seed)
    result = job.result()
    logger.debug("Backend: %s", result.backend_name)

    counts = result.get_counts()
    logger.debug("Counts: %s", counts)
    return counts


def _cut_point_rounds(prep_qubits):
    """``{qubit: number of cut measurements}`` from a list of measured cut qubits.

    A qubit listed twice was cut-measured twice (two gate cuts on one wire, e.g. a
    general gate decomposed into several CZs). A set counts each qubit once.
    """
    rounds = {}
    for q in prep_qubits or ():
        rounds[q] = rounds.get(q, 0) + 1
    return rounds


def _cut_point_sign_qiskit(key, q, rounds, num_qubits):
    """Decomposition sign of qubit ``q`` from a Qiskit counts key.

    Each cut-point outcome ``s`` of ``q`` contributes ``-(-1)**s`` and they
    multiply, as on the Qibo path. The key is ``"<meas> <c{K}> ... <c2> <c>"``:
    ``graph_qiskit`` puts a qubit's k-th cut measurement in register ``c{k}``
    (``c`` for the first), each little-endian.
    """
    cut_registers = key.split()[:0:-1]  # c, c2, c3, ...
    sign = 1
    for k in range(rounds):
        bits = cut_registers[k] if k < len(cut_registers) else ""
        s = int(bits[num_qubits - 1 - q]) if len(bits) == num_qubits else 0
        sign *= -((-1) ** s)
    return sign


def _expectation_from_counts_qiskit(counts, observables, shots=1024, prep_qubits=None):
    """Compute ``<O>`` for a Z/I Pauli string from Qiskit measurement counts.

    ``prep_qubits`` names the (renumbered) cut qubits that were *measured* for the
    gate-cut decomposition in this subcircuit. Such a qubit carries a
    decomposition sign that is normally read off its (forced) final measurement.
    When the user observes that qubit with identity, the final measurement is
    traced out and the sign would be lost -- so instead we recover it from the
    cut-point measurement register ``c`` as ``-(-1)**s`` (``s`` = the cut-point
    outcome). This decouples the internal cut-measurement from the final
    observable, letting identity on a cut qubit reconstruct correctly. When
    ``prep_qubits`` is empty this reduces to the plain Z/I expectation, so the
    all-Z / identity-on-non-cut-qubit paths are unchanged.

    :param counts: The ``get_counts()`` dict from :func:`_counts_qiskit`; keys are
        ``"<measure-all reg> <cut-measure reg c>"`` (space-separated, little-endian).
    :param observables: Per-qubit Pauli string over ``{I, Z}`` (qubit-0-first).
    :param shots: Number of shots the counts were drawn from.
    :param prep_qubits: Renumbered qubit indices measured for the cut, one entry
        per cut measurement (a qubit cut twice is listed twice).
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    prep = _cut_point_rounds(prep_qubits)
    num_qubits = len(observables)
    expectation_value = 0
    for key, value in counts.items():
        # Qiskit get_counts() keys are "<measure-all reg> <cut-measure regs>", each
        # little-endian (qubit 0 is the last bit); index bit q as position
        # num_qubits-1-q to align with the qubit-0-first observable string.
        measure_all = key.split()[0]
        contribution = 1
        for q, obs in enumerate(observables):
            if q in prep:
                # Measured cut qubit: each (I +/- Z)/2 projector contributes the
                # decomposition sign -(-1)**s from its cut-point outcome s, for
                # *any* observable. The observable then multiplies it: Z reads the
                # final measurement, I is trivial. Reading the cut-point (not the
                # final bit) is what makes this correct when the cut qubit keeps
                # evolving after the cut.
                contribution *= _cut_point_sign_qiskit(key, q, prep[q], num_qubits)
                if obs == "Z":
                    contribution *= (-1) ** int(measure_all[num_qubits - 1 - q])
                elif obs != "I":
                    raise ValueError(f"Unsupported observable {obs}")
            elif obs == "Z":
                contribution *= (-1) ** int(measure_all[num_qubits - 1 - q])
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    return expectation_value


def _expec_value_qiskit(
    qc, shots=1024, method="automatic", prep_qubits=None, seed=None
):

    if qc is None:
        return None

    if not isinstance(qc, qiskit.QuantumCircuit):
        observables = qc[1]
        qc = qc[0]
    else:
        observables = "Z" * qc.num_qubits

    if prep_qubits is None:
        # 3-step subcircuit workflow: prep qubits aren't threaded in, so recover
        # them from the subcircuit's mid-circuit (gate-cut) measurements. Lets the
        # workflow use the decoupled reconstruction like the all-in-one path.
        prep_qubits = [
            qc.find_bit(instr.qubits[0]).index
            for instr in qc.data
            if instr.operation.name == "measure"
        ]

    counts = _counts_qiskit(qc, shots=shots, method=method, seed=seed)
    expectation_value = _expectation_from_counts_qiskit(
        counts, observables, shots, prep_qubits=prep_qubits
    )
    logger.debug("Expected value: %s", expectation_value)
    return expectation_value


def _counts_from_sampler(pub_result, circuit):
    """``get_counts()``-style counts from a SamplerV2 result, over all registers.

    SamplerV2 returns one bit array per classical register. Join them per shot
    into Qiskit's usual ``"<last register> ... <first register>"`` keys (each
    little-endian), so :func:`_expectation_from_counts_qiskit` reads the final
    ``meas`` bits and every gate-cut register (``c``, ``c2``, ...) as on the
    simulators.
    """
    names = [register.name for register in reversed(circuit.cregs)]
    columns = [getattr(pub_result.data, name).get_bitstrings() for name in names]
    counts = {}
    for bits in zip(*columns):
        key = " ".join(bits)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _joint_counts_qiskit_qpu(
    dag, num_qubits, prep_sorted, shots, batch, backend, options=None
):
    """Joint final/cut-point counts of a fragment on IBM hardware (SamplerV2):
    the cut-point measurements stay mid-circuit, each in its own register."""
    from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

    built = dag_to_circuit_qiskit(dag, num_qubits)
    if built is None:
        return None
    base = built[0] if isinstance(built, (list, tuple)) else built
    qc = base.copy()
    qc.measure_all()
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    job = Sampler(mode=batch, options=options).run([pm.run(qc)], shots=shots)
    logger.info(f"job id: {job.job_id()}")
    counts = _counts_from_sampler(job.result()[0], qc)
    return _joint_from_qiskit_counts(counts, base, num_qubits, prep_sorted)


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
def _expec_value_qiskit_qpu(
    qc,
    shots=1024,
    method="automatic",
    batch=None,
    backend=None,
    prep_qubits=None,
    options=None,
):
    """Expectation value of a Qiskit subcircuit on IBM hardware (SamplerV2).

    Same reconstruction as :func:`_expec_value_qiskit`: a gate cut's mid-circuit
    measurements carry its decomposition sign, read from their registers. Reading
    only the final ``meas`` register, as before, dropped those signs.
    """
    if qc is None:
        return None

    if not isinstance(qc, qiskit.QuantumCircuit):
        observables = qc[1]
        qc = qc[0]
    else:
        observables = "Z" * qc.num_qubits

    if prep_qubits is None:
        # As on the simulators: the mid-circuit measurements are the cut points.
        prep_qubits = [
            qc.find_bit(instr.qubits[0]).index
            for instr in qc.data
            if instr.operation.name == "measure"
        ]

    qc = qc.copy()  # don't append measure_all to the caller's subcircuit
    qc.measure_all()
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)
    sampler = Sampler(mode=batch, options=options)
    # Pass shots: without it the sampler used its own default while the value
    # was still divided by ``shots``.
    job = sampler.run([isa_circuit], shots=shots)
    logger.info(f"job id: {job.job_id()}")
    counts = _counts_from_sampler(job.result()[0], qc)
    logger.info(counts)
    return _expectation_from_counts_qiskit(
        counts, observables, sum(counts.values()), prep_qubits=prep_qubits
    )


# Runs inside a per-subcircuit @task — keep plain (PyCOMPSs has no nested tasks).
def _expec_value_qiskit_gpu(qc, shots=1024, method="automatic"):
    if qc is None:
        return None

    if not isinstance(qc, qiskit.QuantumCircuit):
        observables = qc[1]
        subcirc = qc[0]
        qc = subcirc
    else:
        observables = ["Z"] * qc.num_qubits
        observables = "".join(observables)

    # Recover the gate-cut measured (prep) qubits from the mid-circuit measurements
    # before appending the final measure_all, so the same decoupled reconstruction
    # as the CPU path applies (the cut measure is collapse-and-continue).
    prep_qubits = [
        qc.find_bit(instr.qubits[0]).index
        for instr in qc.data
        if instr.operation.name == "measure"
    ]
    qc = qc.copy()  # don't append measure_all to the caller's subcircuit
    qc.measure_all()

    simulator = AerSimulator(device="GPU", method=method, cuStateVec_enable=True)

    circ = transpile(qc, simulator)

    job = simulator.run(circ, shots=shots)
    result = job.result()

    counts = result.get_counts()
    logger.info(counts)

    return _expectation_from_counts_qiskit(
        counts, observables, shots, prep_qubits=prep_qubits
    )


def _check_ibm_qc(qpu_config=None):
    """Connect to IBM Quantum and select a backend.

    Thin wrapper over :class:`~Qdislib.core.backend.qpu_backends.IBMBackend`,
    which is the single place that authenticates and selects the backend. Pass
    ``qpu_config`` (a dict of ``token`` / ``channel`` / ``instance`` /
    ``backend_name`` / ``mode`` / ``max_time`` / ``proxies`` / ``service``) for
    explicit configuration; each key falls back to its ``IBM_QUANTUM_*`` env var,
    so env-var setups keep working.

    :return: ``(mode, backend)`` -- what the jobs run in (the backend itself for
        plain jobs, the default; or an opened Batch / Session with ``mode=``) and
        the selected backend.
    """
    from Qdislib.core.backend.registry import build_backend
    from Qdislib.core.backend.base_backend import BackendConfig

    ibm = build_backend(
        "IBM_BACKEND", BackendConfig(num_qubits=0, extra=qpu_config or {})
    )
    return ibm.session, ibm.ibm_backend


def _gate_cutting_subcircuits(
    dag: typing.Any,
    gates_cut: typing.List[typing.Any],
    observables=None,
    software="qiskit",
):
    """
    Partition a quantum circuit into subcircuits via gate cutting.

    This is the entry point of the **recommended three-step workflow**: generate
    the subcircuits, evaluate each one (in your own loop, on whatever backend or
    hardware you like), then recombine the results. It gives you full control over
    *how* each subcircuit is executed, unlike :func:`gate_cutting`, which runs the
    whole pipeline in one call::

        import Qdislib as qd

        # 1. Cut the circuit into independently executable subcircuits.
        subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")

        # 2. Evaluate each subcircuit (here with the built-in helper, but this
        #    loop is yours: dispatch to real QPUs, batch, cache, etc.).
        values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

        # 3. Recombine into the original circuit's expectation value.
        result = subcircuits.reconstruct(values)

    The cut metadata (number of cuts/components) is carried by the returned object,
    so :meth:`~Qdislib.utils.backends.CutSubcircuits.reconstruct` needs nothing else.

    Parameters
    ----------
    dag : Union[qiskit.QuantumCircuit, qibo.models.Circuit, cudaq.PyKernelDecorator, nx.Graph]
        The input quantum circuit. It can be a Qiskit ``QuantumCircuit``, a Qibo
        ``Circuit``, a CUDA-Q kernel, or an already-built DAG/networkx graph.

    gates_cut : List[str]
        Gate identifiers (``"<GATE>_<INDEX>"``, e.g. ``"CZ_2"``) marking where the
        circuit is split into subcircuits. Use :func:`~Qdislib.find_cut` to discover
        good candidates.

    observables : str or list of (float, str), optional
        Per-qubit observable string (e.g. ``"ZZ"``). When given, the measurement
        basis is rotated accordingly before cutting; ``"I"`` entries become identity
        observables. Defaults to the all-``Z`` observable. A weighted Pauli sum
        ``[(coef, "PAULI"), ...]`` gives one subcircuit set per term, and
        ``reconstruct`` returns the weighted sum.

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
        If ``gates_cut`` is malformed or an unsupported ``software`` value is given.

    Notes
    -----
    - If the circuit has multiple connected components, each is processed separately.
    - If ``gates_cut`` is empty the whole circuit is returned as a single subcircuit
      after qubit reindexing.
    - To pass the values to :func:`gate_cutting_subcircuit_reconstruction` directly
      instead of via :meth:`reconstruct`, use ``subcircuits.number_cuts`` and
      ``subcircuits.number_components``.

    Examples
    --------
    >>> import Qdislib as qd
    >>> from qibo import models, gates
    >>> circuit = models.Circuit(2)
    >>> circuit.add(gates.H(0)); circuit.add(gates.CZ(0, 1)); circuit.add(gates.RY(0, 0.8))
    >>> subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")
    >>> values = [qd.expectation_value(s, shots=8192) for s in subcircuits]
    >>> subcircuits.reconstruct(values)  # doctest: +SKIP
    -0.717...
    """
    pattern = re.compile(r"^[A-Za-z]+_\d+$")

    if not isinstance(gates_cut, list) or not all(
        isinstance(g, str) for g in gates_cut
    ):
        raise InvalidCutError(
            "Invalid format of gates to cut: expected a list of strings."
        )
    if not all(pattern.match(g) for g in gates_cut):
        raise InvalidCutError(
            "Invalid gate name format: each entry must match '<GATE>_<NUMBER>', e.g. 'CX_1'."
        )

    # A weighted Pauli sum: one subcircuit set per term (each term rotates its own
    # measurement basis), recombined as sum_k c_k <P_k>.
    if observables is not None and not isinstance(observables, str):
        terms = list(normalize_observables(observables))
        parts = [
            _gate_cutting_subcircuits(
                dag.copy() if isinstance(dag, nx.Graph) else dag,
                gates_cut,
                observables=pauli,
                software=software,
            )
            for _coeff, pauli in terms
        ]
        return CombinedCutSubcircuits(parts, [coeff for coeff, _ in terms])

    dag = convert_to_dag(dag, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    # A non-CZ/CX cut gate is decomposed into CZ + single-qubit gates and its CZs cut.
    dag, gates_cut = _expand_general_gate_cuts(dag, gates_cut)
    _validate_gate_cuts_partition(dag, gates_cut)

    number_cuts = len(gates_cut)

    if nx.number_connected_components(dag.to_undirected(as_view=True)) > 1:
        # Disconnected pieces are cut (or not) independently and multiply: keep each
        # piece's subcircuits as its own part so reconstruct() can't mix them up.
        S = [
            dag.subgraph(c).copy()
            for c in nx.connected_components(dag.to_undirected(as_view=True))
        ]
        parts = []
        for s in S:
            tmp_cuts = [c for c in gates_cut if s.has_node(c)]
            if tmp_cuts:
                piece = _execute_gate_cutting_subcircuits(
                    s, tmp_cuts, software=software
                )
            else:
                s_new, highest_qubit = _update_qubits(s)
                piece = [make_subcircuit(s_new, highest_qubit, software)]
            parts.append(
                CutSubcircuits(
                    piece,
                    number_cuts=len(tmp_cuts),
                    number_components=max(1, len(piece) // (6 ** len(tmp_cuts))),
                    reconstructor=_gate_cutting_subcircuit_reconstruction,
                )
            )
        return CombinedCutSubcircuits(parts)
    elif gates_cut:
        subcircuits = _execute_gate_cutting_subcircuits(
            dag, gates_cut, software=software
        )
    else:
        s_new, highest_qubit = _update_qubits(dag)
        subcircuits = [make_subcircuit(s_new, highest_qubit, software)]

    number_components = max(1, len(subcircuits) // (6**number_cuts))
    return CutSubcircuits(
        subcircuits,
        number_cuts=number_cuts,
        number_components=number_components,
        reconstructor=_gate_cutting_subcircuit_reconstruction,
    )


def _execute_gate_cutting_subcircuits(dag, gates_cut, software):
    new_dag = _generate_cut(dag, gates_cut)
    # The CPU backends recover the cut-point sign in their evaluator (prep qubits
    # auto-detected from the subcircuit's mid-circuit measurements), so the 3-step
    # workflow handles an identity observable on a non-terminal cut qubit too --
    # pass ``software`` so the guard only fires for the paths that genuinely can't.
    _assert_identity_observable_supported(new_dag, software=software)

    if software not in _sdk_names():
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(_sdk_names())}."
        )

    # Each subcircuit is *created* in its own task (_make_subcircuit_task), separate
    # from the per-subcircuit evaluation task -- so under PyCOMPSs these are
    # Subcircuit futures fed straight into expectation_value's eval task (task -> task).
    subcircuits = []
    for index in range(6 ** len(gates_cut)):
        copy_graph = _remove_red_edges(new_dag)
        num_components = nx.number_connected_components(
            copy_graph.to_undirected(as_view=True)
        )
        graph_components = [nx.DiGraph() for _ in range(num_components)]

        # Populates graph_components in place (one DiGraph per component).
        _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        for s in graph_components:
            subcircuits.append(_make_subcircuit_task(s, software))
    return subcircuits


def _gate_cutting_subcircuit_reconstruction(
    expected_values, number_cuts, number_components=2
):
    """
    Reconstruct the final expectation value from subcircuit results after gate cutting.

    This function performs post-processing of expectation values obtained from
    individually simulated or measured subcircuits. It applies the gate-cutting
    sign correction and computes the final observable estimate according to the
    number of cuts made in the original quantum circuit.

    In the recommended workflow you don't call this directly: the object returned
    by :func:`gate_cutting_subcircuits` already carries ``number_cuts`` and
    ``number_components``, so ``subcircuits.reconstruct(values)`` forwards here for
    you. Call this function explicitly only when reconstructing from values whose
    cut metadata you are tracking yourself.

    Parameters
    ----------
    expected_values : List[float]
        A flat list of raw expectation values resulting from the execution
        of subcircuits, in the same order as the subcircuits were generated.
        Values are grouped into chunks of ``number_components``.

    number_cuts : int
        The total number of gate cuts performed in the original circuit.
        This determines the normalization factor used in the final result.

    number_components : int
        The number of independent parts the initial subcircuits has.

    Returns
    -------
    float
        The final reconstructed expectation value of the original circuit,
        computed as the sum of sign-adjusted subcircuit values divided by
        2 raised to the power of the number of cuts.

    Notes
    -----
    - The `expected_values` list is divided into pairs, and each pair is processed
      using the `_change_sign_gate_cutting` function, which applies sign corrections
      based on the cut structure.
    - The computation uses `compss_wait_on()` to ensure all results are synchronized
      before aggregation. This is relevant in distributed or parallel execution contexts (using PyCOMPSs).
    - The final result scales by a factor of `1 / (2 ** number_cuts)`, consistent with
      the probabilistic interpretation of cut circuits.

    Examples
    --------
    >>> raw_values = [0.5, -0.5, 0.3, -0.3, 0.5, -0.5, 0.3, -0.3, 0.5, -0.5, 0.3, -0.3]
    >>> number_of_cuts = 1
    >>> reconstructed_value = _gate_cutting_subcircuit_reconstruction(raw_values, number_of_cuts, number_components=2)
    >>> print(reconstructed_value)
    0.0
    """
    from Qdislib.core.cutting_algorithms._reconstruction import signed_terms_sum

    return signed_terms_sum(
        expected_values, number_cuts, number_components, _change_sign_gate_cutting
    )


# --------------------------------------------------------------------------- #
# Public API (the calling convention shared by every cut type, see _api.py)
# --------------------------------------------------------------------------- #


def _legacy_run(legacy, custom_strategy=False):
    """True when a call needs the pre-2.0 implementation: the deprecated execution
    options (``qpu=``, ``gpu=``, ...) or a custom ``CacheStrategy``. Everything else
    runs through :func:`circuit_cutting`, with its cache (the qibo platform of a
    ``backend="numpy"``-style call is not honoured there: it is the ``"qibo"``
    simulator)."""
    return bool(set(legacy or ()) - {"platform"}) or bool(custom_strategy)


def _route_to_circuit_cutting(legacy, custom_strategy=False):
    if _legacy_run(legacy, custom_strategy):
        raise ValueError(
            "shots=None (exact), optimized=True and num_samples run through "
            "circuit_cutting, which takes backend= only (no qpu/gpu options, no "
            "custom cache strategy)."
        )


def _cache_mode_option(cache_mode):
    """True for a custom ``CacheStrategy`` (kept on the 1.x path for one release);
    ``cache_mode`` is otherwise ignored: the engine caches each circuit's whole
    result (counts or exact correlators), which serves every observable."""
    import warnings

    if cache_mode is None:
        return False
    from Qdislib.utils.circuit_cache import CacheStrategy

    custom = isinstance(cache_mode, CacheStrategy)
    warnings.warn(
        "cache_mode is deprecated and removed in the next release: the cache keeps "
        "each fragment circuit's counts (or exact correlators), which serve every "
        "observable"
        + (
            "; a custom CacheStrategy still runs the 1.x implementation."
            if custom
            else ", so it is ignored."
        ),
        DeprecationWarning,
        stacklevel=3,
    )
    return custom


def _legacy_store(cache):
    """A store object for the 1.x path from what ``cache=`` takes."""
    if cache is None or cache is False:
        return None
    from Qdislib.core.cutting_algorithms.engine_cache import _store

    return _store(cache)


@standard_arguments
def gate_cutting(
    circuit,
    cuts,
    observable=None,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    optimized=False,
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
    """Reconstruct ``<observable>`` from a circuit cut at two-qubit gates
    (measure-and-prepare gate cut, 6 terms per gate: ``gamma = 3`` for a CZ/CX,
    ``1 + 2|sin(phi)|`` for a rotation such as ``rzz(phi)``).

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Gate names to cut (``["CZ_3", ...]``); a gate other than CZ/CX is
        cut as the ``rzz`` rotations it decomposes into, each at its own angle (one
        for ``cp``, ``crz``, ``rxx``, ...; up to three for a general gate).
    :param observable: A Pauli string (``X``/``Y`` are rotated into ``Z``) or a
        weighted Pauli sum ``[(coefficient, string), ...]``; ``None`` is all ``Z``.
    :param backend: Where the subcircuits run: an SDK (``"qibo"``, ``"qiskit"``,
        ``"cudaq"``, ``"pennylane"``, its local simulator), a registered backend
        name (``"IBM_BACKEND"``, ``"IQM_BACKEND"``, ...) or a
        :func:`~Qdislib.build_backend` object; ``None`` is the circuit's own SDK.
    :param shots: Shots per subcircuit; ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling (each subcircuit gets its own
        stream); nothing is random with ``shots=None``.
    :param sync: ``False`` returns the value as a future (PyCOMPSs).
    :param optimized: Run each fragment once per distinct setting of its cuts (5
        circuits per side of a cut gate instead of 6), see :func:`circuit_cutting`.
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
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
    :raises InvalidCutError: For a cut that is not a gate of the circuit, is given
        twice or does not separate its two sides.
    :raises ValueError: For an observable without one letter per qubit, ``shots`` below
        1, or ``shots=None`` with a backend that is no simulator.
    :return: The reconstructed expectation value (a list of them with
        ``parameter_values``).

    Example::

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)  # gate 2, "CZ_2"
        qc.cz(1, 2)  # gate 3, "CZ_3"
        value = qd.gate_cutting(qc, ["CZ_2"], "ZZZ", shots=4096, seed=1)
    """
    custom = _cache_mode_option(cache_mode)
    if (
        not _legacy_run(_legacy, custom)
        or shots is None
        or optimized
        or num_samples is not None
    ):
        from Qdislib.core.cutting_algorithms.mixed_cutting import circuit_cutting

        _route_to_circuit_cutting(_legacy, custom)
        return circuit_cutting(
            circuit,
            {"gate": cuts},
            observable,
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            optimized=optimized,
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
    return _gate_cutting(
        circuit,
        cuts,
        observables=default_observable(observable, circuit),
        shots=shots,
        method=legacy.get("method", "automatic"),
        sync=sync,
        gpu=legacy.get("gpu", False),
        gpu_min_qubits=legacy.get("gpu_min_qubits"),
        qpu=legacy.get("qpu", False),
        qpu_dict=legacy.get("qpu_dict"),
        qpu_config=legacy.get("qpu_config"),
        software=resolve_backend(backend, circuit),
        cache=_legacy_store(cache),
        cache_mode=cache_mode or "expectation",
        seed=seed,
    )


@standard_arguments
def gate_cutting_sampling(
    circuit,
    cuts,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    optimized=False,
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
    two-qubit gates.

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Gate names to cut, as in :func:`gate_cutting`.
    :param backend: Where the subcircuits run (see :func:`gate_cutting`).
    :param shots: Shots per subcircuit; ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling (each subcircuit gets its own
        stream).
    :param sync: ``False`` returns the distribution as a future (PyCOMPSs).
    :param optimized: Run each fragment once per distinct setting of its cuts, see
        :func:`circuit_cutting`.
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
    :raises InvalidCutError: For a cut that is not a gate of the circuit, is given
        twice or does not separate its two sides.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings), a
        :class:`~Qdislib.FactoredDistribution` over more than 24 qubits; a list of them
        with ``parameter_values``.
    """
    if (
        not _legacy_run(_legacy)
        or shots is None
        or optimized
        or num_samples is not None
    ):
        from Qdislib.core.cutting_algorithms.mixed_cutting import (
            circuit_cutting_sampling,
        )

        _route_to_circuit_cutting(_legacy)
        return circuit_cutting_sampling(
            circuit,
            {"gate": cuts},
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            optimized=optimized,
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
    return _gate_cutting_sampling(
        circuit,
        cuts,
        shots=shots,
        software=resolve_backend(backend, circuit),
        sync=sync,
        qpu=legacy.get("qpu", False),
        qpu_dict=legacy.get("qpu_dict"),
        seed=seed,
    )


@standard_arguments
def gate_cutting_subcircuits(
    circuit,
    cuts,
    observable=None,
    *,
    _software=None,
    optimized=False,
    num_samples=None,
    seed=None,
):
    """The subcircuits of :func:`gate_cutting`, to run yourself (three-step)::

        subs = qd.gate_cutting_subcircuits(circuit, ["CZ_3"], "ZZZZ")
        values = [qd.expectation_value(s, shots=4096) for s in subs]
        value = subs.reconstruct(values)

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Gate names to cut, as in :func:`gate_cutting`.
    :param observable: A Pauli string or a weighted Pauli sum; ``None`` is all ``Z``.
    :param optimized: Run each fragment once per distinct setting of its cuts, see
        :func:`circuit_cutting`.
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :param seed: Seed for the ``num_samples`` draws.
    :raises InvalidCutError: For a cut that is not a gate of the circuit, is given
        twice or does not separate its two sides.
    :return: A :class:`~Qdislib.CutCircuits` (a list of subcircuits with
        ``reconstruct(values)``); for a Pauli sum a
        :class:`~Qdislib.CombinedCutSubcircuits`.
    """
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        circuit_cutting_subcircuits,
    )

    return circuit_cutting_subcircuits(
        circuit,
        {"gate": cuts},
        observable,
        _software=_software,
        optimized=optimized,
        num_samples=num_samples,
        seed=seed,
    )


@standard_arguments
def gate_cutting_sampling_subcircuits(
    circuit, cuts, *, _software=None, optimized=False, num_samples=None, seed=None
):
    """The subcircuits of :func:`gate_cutting_sampling`, to sample yourself::

        subs = qd.gate_cutting_sampling_subcircuits(circuit, ["CZ_3"])
        results = [qd.sample_subcircuit(s, shots=8192) for s in subs]
        distribution = subs.reconstruct(results)

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Gate names to cut, as in :func:`gate_cutting`.
    :param optimized: Run each fragment once per distinct setting of its cuts, see
        :func:`circuit_cutting`.
    :param num_samples: Draw this many configurations by weight instead of running
        them all (Monte-Carlo sampling, see :func:`circuit_cutting`).
    :param seed: Seed for the ``num_samples`` draws.
    :raises InvalidCutError: For a cut that is not a gate of the circuit, is given
        twice or does not separate its two sides.
    :return: A :class:`~Qdislib.CutCircuits` (a list of subcircuits with
        ``reconstruct(results)``).
    """
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        circuit_cutting_sampling_subcircuits,
    )

    return circuit_cutting_sampling_subcircuits(
        circuit,
        {"gate": cuts},
        _software=_software,
        optimized=optimized,
        num_samples=num_samples,
        seed=seed,
    )


# Deprecated in 2.0.0: ``subcircuits.reconstruct(results)`` is the one way.


def gate_cutting_subcircuit_reconstruction(
    expected_values, number_cuts, number_components=2
):
    """Deprecated: use ``subcircuits.reconstruct(results)``.

    :param expected_values: One value per subcircuit of a 1.x gate cut, in order.
    :param number_cuts: The gates cut.
    :param number_components: The subcircuits of each of its ``6**number_cuts``
        configurations.
    :raises TypeError: For the results of a 2.0 list, which only its ``reconstruct``
        can combine.
    :return: The reconstructed expectation value.
    """
    import warnings

    warnings.warn(
        "gate_cutting_subcircuit_reconstruction() is deprecated: use subcircuits.reconstruct(results); it is "
        "removed in the next release.",
        DeprecationWarning,
        stacklevel=2,
    )
    from Qdislib.core.cutting_algorithms._reconstruction import reject_engine_results

    reject_engine_results(expected_values, "gate_cutting_subcircuit_reconstruction")
    return _gate_cutting_subcircuit_reconstruction(
        expected_values, number_cuts, number_components
    )


def gate_cutting_sampling_reconstruction(
    sampled_results, number_cuts, number_components, num_qubits
):
    """Deprecated: use ``subcircuits.reconstruct(results)``.

    :param sampled_results: One :func:`sample_subcircuit` result per subcircuit of a
        1.x gate cut, in order.
    :param number_cuts: The gates cut.
    :param number_components: The subcircuits of each configuration.
    :param num_qubits: The circuit's qubits.
    :raises TypeError: For the results of a 2.0 list, which only its ``reconstruct``
        can combine.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    import warnings

    warnings.warn(
        "gate_cutting_sampling_reconstruction() is deprecated: use subcircuits.reconstruct(results); it is "
        "removed in the next release.",
        DeprecationWarning,
        stacklevel=2,
    )
    from Qdislib.core.cutting_algorithms._reconstruction import reject_engine_results

    reject_engine_results(sampled_results, "gate_cutting_sampling_reconstruction")
    return _gate_cutting_sampling_reconstruction(
        sampled_results, number_cuts, number_components, num_qubits
    )
