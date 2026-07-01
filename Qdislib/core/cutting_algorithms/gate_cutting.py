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
    "gate_cutting_subcircuits",
    "gate_cutting_subcircuit_reconstruction",
]

import networkx as nx
import math
import numpy as np

try:
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit_aer import AerSimulator
except ImportError:  # qiskit / qiskit-aer are optional
    generate_preset_pass_manager = AerSimulator = None

import re

from Qdislib.utils.runtime import task, compss_wait_on, constraint, COLLECTION_IN


from Qdislib.utils.graph_qibo import (
    _update_qubits,
    _remove_red_edges,
    _max_qubits_graph,
    dag_to_circuit_qibo,
)
from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
from Qdislib.core.cutting_algorithms.wire_cutting import _sum_results
from Qdislib.core.find_cut.find_cut import _find_nodes_with_qubit

try:
    from qiskit_ibm_runtime import (
        QiskitRuntimeService,
        Batch,
        SamplerV2 as Sampler,
    )
except ImportError:  # qiskit-ibm-runtime optional (only used for IBM QPU)
    QiskitRuntimeService = Batch = Sampler = None


try:
    import qiskit
    from qiskit import transpile

    qiskit_available = True
except ImportError:  # qiskit is optional; functions that need it fail when called
    qiskit = transpile = None
    qiskit_available = False

try:
    import qibo
    from qibo import gates

    qibo_available = True
except ImportError:  # qibo is optional; functions that need it fail when called
    qibo = gates = None
    qibo_available = False


import os
import typing


from Qdislib.utils.backends import (
    convert_to_dag,
    change_basis_dag,
    normalize_observables,
    make_subcircuit,
    CutSubcircuits,
)
from Qdislib.utils.exceptions import UnsupportedBackendError, InvalidCutError


def gate_cutting(
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
    software: str = "qiskit",
    cache=None,
    cache_mode="expectation",
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
    return_subcircuits : bool, default=False
        If True, also return the list of subcircuits generated after cutting.
    gpu : bool, default=False
        If True, attempt to run circuit evaluations on GPU simulators when possible.
    gpu_min_qubits : int, optional
        Minimum number of qubits needed to offload execution to GPU simulators. Default is 0.
    qpu : bool, default=False
        If True, attempt to run circuit evaluations on real Quantum Processing Units (QPUs).
    qpu_dict : dict, optional
        Dictionary specifying QPU backends and their maximum supported qubits.
        Example: `{"IBM_Quantum": 153, "MN_Ona": 5}`.
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
        from mymodule import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)     # Gate "H_1"
        qc.cz(0, 1) # Gate "CZ_2"
        qc.h(1)     # Gate "H_3"

        # Define the gate to cut (e.g., the CZ gate)
        gates_to_cut = ["CZ_2"]  # You must identify the gate in order of definition

        reconstruction = gate_cutting(qc, gates_to_cut, shots=2048)

    """

    # A weighted sum of Pauli observables (a Hamiltonian) is evaluated term by
    # term: <sum_i c_i P_i> = sum_i c_i <P_i>. Each term reuses the single
    # Pauli-string path below. A plain string (or None) skips this and runs once.
    if observables is not None and not isinstance(observables, str):
        total = 0.0
        for coeff, pauli in normalize_observables(observables):
            base = dag.copy() if isinstance(dag, nx.Graph) else dag
            total += coeff * gate_cutting(
                base,
                gates_cut,
                observables=pauli,
                shots=shots,
                method=method,
                sync=sync,
                gpu=gpu,
                gpu_min_qubits=gpu_min_qubits,
                qpu=qpu,
                qpu_dict=qpu_dict,
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

        subcircuits = gate_cutting_subcircuits(
            dag, gates_cut, observables=observables, software=software
        )
        values = [
            expectation_value(s, shots=shots, cache=cache, cache_mode=cache_mode)
            for s in subcircuits
        ]
        return subcircuits.reconstruct(values)

    # Normalise optional dict/int args before they are indexed below.
    if qpu_dict is None:
        qpu_dict = {}
    if gpu_min_qubits is None:
        gpu_min_qubits = 0

    if qpu and "IBM_Quantum" in qpu_dict:
        batch, backend = _check_ibm_qc()
    else:
        batch, backend = None, None

    dag = convert_to_dag(dag, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy() for c in nx.connected_components(dag.to_undirected())
        ]
        results = []
        for s in S:
            # num_qubits = _max_qubit(s)
            tmp_cuts = []
            for c in gates_cut:
                if s.has_node(c):
                    tmp_cuts.append(c)
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
                    )
                )
        if sync:
            results = compss_wait_on(results)
        final_recons = 1 / (2 ** len(gates_cut)) * math.prod(results)
        logger.debug(
            "Reconstructed expected value %s from %s subcircuit values",
            final_recons,
            len(results),
        )
        return final_recons
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
            )
            if sync:
                results = compss_wait_on(results)
            final_recons = 1 / (2 ** len(gates_cut)) * sum(results)
            logger.debug(
                "Reconstructed expected value %s from %s subcircuit values",
                final_recons,
                len(results),
            )
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
            )
            if sync:
                final_recons = compss_wait_on(final_recons)
            logger.debug("No cut performed; expected value: %s", final_recons)
        return final_recons


def _generate_cut(dag, gates_cut):
    dag_copy = dag.copy()
    dag_copy = _remove_red_edges(dag_copy)
    for index, gate_name in enumerate(gates_cut, start=1):
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
            f"SUBS1_{index}", gate="S", qubits=(target_qubits[0],), parameters=()
        )
        dag_copy.add_node(
            f"SUBS2_{index}", gate="S", qubits=(target_qubits[1],), parameters=()
        )

        if pred_0:
            dag_copy.add_edge(pred_0[0], f"SUBS1_{index}", color="blue")

        if succ_0:
            dag_copy.add_edge(f"SUBS1_{index}", succ_0[0], color="blue")

        if pred_1:
            dag_copy.add_edge(pred_1[0], f"SUBS2_{index}", color="blue")

        if succ_1:
            dag_copy.add_edge(f"SUBS2_{index}", succ_1[0], color="blue")

    return dag_copy


# CPU-simulation backends whose expectation reconstruction decouples the internal
# cut-measurement sign from the final observable (see
# :func:`_expectation_from_counts_qiskit`), and so handle an identity observable
# on any cut qubit correctly. All three CPU backends do; the guard below still
# applies to GPU/QPU evaluation paths (not decoupled) and to the 3-step
# subcircuit workflow (which evaluates without the decoupled reconstruction).
_IDENTITY_CUT_DECOUPLED = {"qiskit", "qibo", "cudaq"}


def _assert_identity_observable_supported(cut_dag, software=None, gpu=False, qpu=False):
    """Reject identity-observable cut cases a path cannot yet reconstruct.

    A gate cut measures its cut qubits as part of the measure-and-prepare
    quasiprobability decomposition, and that measurement carries the
    reconstruction sign. If the cut qubit ends on single-qubit gates it stays
    isolated after the cut, so an identity observable on it is a trivial
    single-qubit subcircuit and always works. But if the cut qubit feeds a
    **multi-qubit** gate *after* the cut, an identity there would suppress the
    load-bearing measurement if it were read off the final measurement.

    All three CPU backends (:data:`_IDENTITY_CUT_DECOUPLED`) recover that sign
    from the cut-point measurement (decoupling it from the final observable) and
    handle this correctly, so the check is skipped for them on CPU. It still fires
    for GPU/QPU evaluation paths and the 3-step subcircuit workflow (``software``
    unset), which don't use the decoupled reconstruction -- raising rather than
    returning a silently-wrong value.

    :param cut_dag: DAG after :func:`_generate_cut` (carries the ``SUBS*`` nodes).
    :param software: CPU-simulation backend that will evaluate the subcircuits.
    :param gpu: Whether GPU evaluation is in use (not covered by the fix yet).
    :param qpu: Whether QPU evaluation is in use (not covered by the fix yet).
    :raises InvalidCutError: if an identity observable sits on a cut qubit that
        feeds a multi-qubit gate after the cut, on a not-yet-supported path.
    """
    if software in _IDENTITY_CUT_DECOUPLED and not gpu and not qpu:
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
                f"This gate-cut path cannot honour an identity ('I') observable on "
                f"cut qubit {cut_qubit}: it feeds a multi-qubit gate after the cut, "
                f"so the cut's measure-and-prepare sign cannot be dropped. It is "
                f"supported by the all-in-one gate_cutting() on CPU (qiskit, qibo, "
                f"cudaq), which recovers the sign from the cut-point measurement; "
                f"GPU/QPU and the subcircuit workflow do not yet. Use a CPU run, "
                f"wire cutting, or a cut whose qubits are not measured with identity."
            )


def _decimal_to_base6(num):
    if num == 0:
        return "0"

    base6 = ""
    while num > 0:
        base6 = str(num % 6) + base6
        num //= 6
    return base6


# Plain function (not a @task): the gate-cut substitution is cheap graph surgery,
# run on the master, populating graph_components in place with concrete DiGraphs.
# Keeping it task-free (no COLLECTION_OUT futures) lets the generated subcircuits
# be wrapped as serializable Subcircuit objects under PyCOMPSs.
def _generate_gate_cutting(updated_dag, gates_cut, index, graph_components):

    base6_rep = _decimal_to_base6(index)
    base6_rep = base6_rep.zfill(len(gates_cut))
    list_substitutions = list(map(int, base6_rep))

    # print(list_substitutions)
    for idx2, idx in enumerate(list_substitutions, start=1):
        list_succ1 = list(updated_dag.succ[f"SUBS1_{idx2}"])
        list_pred1 = list(updated_dag.pred[f"SUBS1_{idx2}"])

        list_succ2 = list(updated_dag.succ[f"SUBS2_{idx2}"])
        list_pred2 = list(updated_dag.pred[f"SUBS2_{idx2}"])

        # 1 - Rz(-pi/2) -- Rz(-pi/2)
        if idx == 0:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS1_{idx2}"]["parameters"] = [-np.pi / 2]
            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS2_{idx2}"]["parameters"] = [-np.pi / 2]

        # 2 - Z Rz(-pi/2) -- Z Rz(-pi/2)
        elif idx == 1:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "z"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "z"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        # 3 - MEASURE Rz(-pi/2) -- Rz(-pi)
        elif idx == 2:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS2_{idx2}"]["parameters"] = [-np.pi]

        # 4 - MEASURE Rz(-pi/2) -- RES
        elif idx == 3:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS11_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS1_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ1:
                succ1 = list_succ1[0]
                updated_dag.remove_edge(f"SUBS1_{idx2}", succ1)
                updated_dag.add_edge(f"SUBS11_{idx2}", succ1, color="blue")

            updated_dag.add_edge(f"SUBS1_{idx2}", f"SUBS11_{idx2}", color="blue")

            updated_dag.remove_node(f"SUBS2_{idx2}")

            if list_succ2 and list_pred2:
                succ2 = list_succ2[0]
                pred2 = list_pred2[0]

                updated_dag.add_edge(pred2, succ2)

        # 5 - Rz(-pi) -- MEASURE Rz(-pi/2)
        elif idx == 4:
            updated_dag.nodes[f"SUBS1_{idx2}"]["gate"] = "rz"
            updated_dag.nodes[f"SUBS1_{idx2}"]["parameters"] = [-np.pi]

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        # 6 - RES -- MEASURE Rz(-pi/2)
        elif idx == 5:
            updated_dag.remove_node(f"SUBS1_{idx2}")

            if list_succ1 and list_pred1:
                succ1 = list_succ1[0]
                pred1 = list_pred1[0]
                updated_dag.add_edge(pred1, succ1, color="blue")

            updated_dag.nodes[f"SUBS2_{idx2}"]["gate"] = "measure"
            updated_dag.add_node(
                f"SUBS22_{idx2}",
                gate="rz",
                qubits=updated_dag.nodes[f"SUBS2_{idx2}"].get("qubits"),
                parameters=([-np.pi / 2]),
            )

            if list_succ2:
                succ2 = list_succ2[0]
                updated_dag.remove_edge(f"SUBS2_{idx2}", succ2)
                updated_dag.add_edge(f"SUBS22_{idx2}", succ2, color="blue")

            updated_dag.add_edge(f"SUBS2_{idx2}", f"SUBS22_{idx2}", color="blue")

        else:
            raise ValueError(f"Unexpected gate-cut substitution index: {idx}")

    for i, c in enumerate(nx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")

    return updated_dag


# @constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=1, expectation_value=COLLECTION_IN)
def _change_sign_gate_cutting(expectation_value, index):
    expectation_value = [x for x in expectation_value if x is not None]
    expectation_value = math.prod(expectation_value)
    number = index

    change_sign = False

    while number != 0:
        digit = number % 6  # Get the last digit
        if digit in {3, 5}:  # Check if the digit is 3, 5
            change_sign = not change_sign  # Flip the sign change flag
        number //= 6  # Move to the next digit

    # If change_sign is True, we flip the sign of the original number
    if change_sign:
        return -expectation_value
    else:
        return expectation_value


def _measured_cut_qubits(dag):
    """Renumbered qubits that a gate-cut term *measures* (``SUBS`` -> ``measure``).

    Their measure-and-prepare sign must be recovered from the cut-point
    measurement when the user observes them with identity -- see
    :func:`_expectation_from_counts_qiskit`. Returns an empty set for a ``None``
    component (an empty gate-cut term, tolerated downstream).
    """
    if dag is None:
        return frozenset()
    return frozenset(
        (dag.nodes[n].get("qubits") or (None,))[0]
        for n in dag
        if dag.nodes[n].get("gate") == "measure"
    ) - {None}


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


def _sample_gate_cut_configs(dag, gates_cut, software, shots):
    """Sample every 6**n_cuts gate-cut configuration.

    For each configuration and each connected component, returns the component's
    qubit-0-first counts together with the component's *original* qubit indices,
    so the reconstruction can place each bit back on the right qubit (the piece
    the naive bitstring-concatenation reconstruction got wrong).
    """
    from Qdislib.core.backend.registry import build_execution_backend

    new_dag = _generate_cut(dag, gates_cut)
    configs = []
    for index in range(6 ** len(gates_cut)):
        copy_graph = _remove_red_edges(new_dag.copy())
        num_components = nx.number_connected_components(copy_graph.to_undirected())
        graph_components = [nx.DiGraph() for _ in range(num_components)]
        _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        per_component = []
        for comp in graph_components:
            original = _component_original_qubits(comp)
            s_new, highest_qubit = _update_qubits(comp)
            backend = build_execution_backend(software, num_qubits=highest_qubit)
            counts = backend.evaluate_samples(s_new, highest_qubit, shots)
            per_component.append((counts, original))
        configs.append(per_component)
    return configs


def _reconstruct_gate_cut_distribution(configs, num_cuts, num_qubits):
    """Reconstruct the output distribution from sampled gate-cut configurations.

    Each configuration contributes ``sign(index) / 2**n_cuts`` times the tensor
    product of its component distributions (built over *original* qubit indices).
    Gate cutting yields quasi-probabilities, so the signed sum is abs-normalised
    into a valid distribution. Validated against the exact statevector
    distribution (TVD ~ shot noise) for 1 and 2 cuts.
    """
    accumulated = {}
    for index, per_component in enumerate(configs):
        # Sign from the base-6 digits (3 or 5 flip), mirroring _change_sign_gate_cutting.
        sign, digits = 1, index
        while digits:
            if digits % 6 in (3, 5):
                sign = -sign
            digits //= 6

        combined = {(): 1.0}  # assignment (tuple of (qubit, bit)) -> probability
        for counts, original in per_component:
            total = sum(counts.values()) or 1
            component_dist = {}
            for bitstring, count in counts.items():
                assignment = tuple((original[i], bitstring[i]) for i in range(len(original)))
                component_dist[assignment] = component_dist.get(assignment, 0.0) + count / total
            merged = {}
            for prefix, p_prefix in combined.items():
                for assignment, p_assignment in component_dist.items():
                    key = prefix + assignment
                    merged[key] = merged.get(key, 0.0) + p_prefix * p_assignment
            combined = merged

        weight = sign / (2 ** num_cuts)
        for assignment, prob in combined.items():
            bits = dict(assignment)
            bitstring = "".join(bits.get(q, "0") for q in range(num_qubits))
            accumulated[bitstring] = accumulated.get(bitstring, 0.0) + weight * prob

    norm = sum(abs(v) for v in accumulated.values()) or 1
    return {k: abs(v) / norm for k, v in accumulated.items()}


def gate_cutting_sampling(circuit, gates_cut, shots=1024, software="qiskit"):
    """Reconstruct the output **probability distribution** of a gate-cut circuit.

    Unlike :func:`gate_cutting` (which returns an expectation value), this samples
    each subcircuit and reconstructs the full distribution over measurement
    bitstrings. Every subcircuit is sampled independently -- the reconstruction is
    embarrassingly parallel -- and combined via the gate-cut quasiprobability
    decomposition, then abs-normalised into a valid distribution.

    :param circuit: A Qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: The gates to cut (as in :func:`gate_cutting`).
    :param shots: Measurement shots per subcircuit. Sampling reconstruction has a
        quasiprobability overhead, so it needs more shots than an expectation.
    :param software: Execution backend, or the name of any registered backend.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    dag = convert_to_dag(circuit)
    num_qubits = _dag_num_qubits(dag)

    if not gates_cut:
        from Qdislib.core.backend.registry import build_execution_backend

        backend = build_execution_backend(software, num_qubits=num_qubits)
        counts = backend.evaluate_samples(dag, num_qubits, shots)
        total = sum(counts.values()) or 1
        return {k: v / total for k, v in counts.items()}

    configs = _sample_gate_cut_configs(dag, gates_cut, software, shots)
    return _reconstruct_gate_cut_distribution(configs, len(gates_cut), num_qubits)


@task(returns=1)
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
):
    """Build one gate-cut subcircuit and return its expectation value.

    A single task so that, under PyCOMPSs, the native circuit is built **and**
    executed in one worker process (only the picklable component DAG crosses the
    task boundary, never the fragile Qibo collapse-symbol circuit / CUDA-Q
    kernel); the inner ``dag_to_circuit_*`` / ``_expec_value_*`` calls run inline.

    Routes to a QPU or GPU when one is configured and large enough for the
    subcircuit; otherwise builds and runs on the chosen ``software`` backend
    (Qiskit Aer / Qibo / CUDA-Q).

    :param component: A connected-component DAG of the cut circuit.
    :param max_qubit: Qubit count of the component (for device-size checks).
    :param software: CPU-simulation backend ("qiskit", "qibo" or "cudaq").
    :return: The subcircuit's expectation value (a PyCOMPSs future under COMPSs).
    """
    qpu_dict = qpu_dict or {}
    s_new, highest_qubit = _update_qubits(component)
    prep_qubits = _measured_cut_qubits(s_new)

    # Hardware QPU paths (specialised; their backend extraction is a follow-up).
    if qpu and "MN_Ona" in qpu_dict and qpu_dict["MN_Ona"] >= max_qubit:
        subcirc = dag_to_circuit_qibo(s_new, highest_qubit)
        return _expec_value_qibo_qpu(subcirc, shots=shots, method=method)
    if qpu and "IBM_Quantum" in qpu_dict and qpu_dict["IBM_Quantum"] >= max_qubit:
        subcirc = dag_to_circuit_qiskit(s_new, highest_qubit)
        return _expec_value_qiskit_qpu(
            subcirc, shots=shots, batch=batch, backend=backend
        )

    # Simulator path: pick the backend (GPU-eligible Aer, else the ``software``
    # one) and delegate -- no per-backend branching in the algorithm itself.
    from Qdislib.core.backend.registry import build_execution_backend

    use_gpu = bool(gpu and gpu_min_qubits <= max_qubit and max_qubit <= 30)
    exec_backend = build_execution_backend(
        software, use_gpu=use_gpu, num_qubits=highest_qubit
    )
    return exec_backend.evaluate_expectation(
        s_new, highest_qubit, prep_qubits, shots, method
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
):
    new_dag = _generate_cut(dag, gates_cut)
    _assert_identity_observable_supported(
        new_dag, software=software, gpu=bool(gpu), qpu=bool(qpu)
    )

    S = [
        new_dag.subgraph(c).copy()
        for c in nx.connected_components(new_dag.to_undirected())
    ]

    max_qubit_list = []
    for s in S:
        _max_qubit = _max_qubits_graph(s)
        max_qubit_list.append(_max_qubit)

    reconstruction = []
    for index in range(6 ** len(gates_cut)):
        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(copy_graph.to_undirected())

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())

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
                )
            )

        exp = _change_sign_gate_cutting(exp_values, index)
        reconstruction.append(exp)

    return reconstruction


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
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expectation_from_counts_qibo,
    )

    result = submit_and_wait(subcirc)
    freq = result.counts()
    # Normalise by the actual shot count returned by the QPU (its num_bins may
    # differ from the classical-side ``shots`` argument).
    total_shots = sum(freq.values()) or 1
    return _expectation_from_counts_qibo(freq, observables, total_shots)


@constraint(processors=[{"processorType": "CPU", "computingUnits": "1"}])
@task(returns=1)
def _counts_qiskit(qc, shots=1024, method="automatic"):
    """Simulate a Qiskit circuit with Aer and return its measurement counts.

    Appends ``measure_all`` in place -- so the circuit is single-use (evaluate
    once), consistent with the Qibo evaluator -- then runs it on the CPU Aer
    simulator. Backend/method details are logged at DEBUG level.

    :param qc: A Qiskit ``QuantumCircuit``.
    :param shots: Number of measurement shots.
    :param method: Aer simulation method.
    :return: The Aer ``get_counts()`` dict (keys are ``"<measure-all reg> ..."``).
    """
    qc.measure_all()

    import subprocess

    try:
        subprocess.check_output("nvidia-smi")
        logger.debug("Nvidia GPU detected!")
    except Exception:
        logger.debug("No Nvidia GPU in system!")

    logger.debug("Simulation using CPUs; method: %s", method)
    simulator = AerSimulator(
        device="CPU",
        method=method,
        max_parallel_threads=1,
        mps_omp_threads=1,
        mps_parallel_threshold=1,
    )

    job = simulator.run(qc, shots=shots)
    result = job.result()
    logger.debug("Backend: %s", result.backend_name)

    counts = result.get_counts()
    logger.debug("Counts: %s", counts)
    return counts


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
    :param prep_qubits: Iterable of renumbered qubit indices measured for the cut.
    :raises ValueError: If an unsupported observable character is encountered.
    :return: The expectation value.
    """
    prep = set(prep_qubits or ())
    num_qubits = len(observables)
    expectation_value = 0
    for key, value in counts.items():
        # Qiskit get_counts() keys are "<measure-all reg> <cut-measure reg>", each
        # little-endian (qubit 0 is the last bit); index bit q as position
        # num_qubits-1-q to align with the qubit-0-first observable string.
        registers = key.split()
        measure_all = registers[0]
        cut_measure = registers[1] if len(registers) > 1 else ""
        contribution = 1
        for q, obs in enumerate(observables):
            if q in prep and obs == "I":
                # Identity on a measured cut qubit: recover the decomposition sign
                # from the cut-point outcome rather than the (traced-out) final bit.
                s = int(cut_measure[num_qubits - 1 - q]) if len(cut_measure) == num_qubits else 0
                contribution *= -((-1) ** s)
            elif obs == "Z":
                contribution *= (-1) ** int(measure_all[num_qubits - 1 - q])
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    return expectation_value


def _expec_value_qiskit(qc, shots=1024, method="automatic", prep_qubits=None):

    if qc is None:
        return None

    if type(qc) != qiskit.circuit.quantumcircuit.QuantumCircuit:
        observables = qc[1]
        qc = qc[0]
    else:
        observables = "Z" * qc.num_qubits

    counts = _counts_qiskit(qc, shots=shots, method=method)
    expectation_value = _expectation_from_counts_qiskit(
        counts, observables, shots, prep_qubits=prep_qubits
    )
    logger.debug("Expected value: %s", expectation_value)
    return expectation_value


@task(returns=1)
def _expec_value_qiskit_qpu(
    qc, shots=1024, method="automatic", batch=None, backend=None
):

    if qc is None:
        return None

    if type(qc) != qiskit.circuit.quantumcircuit.QuantumCircuit:
        observables = qc[1]
        subcirc = qc[0]
        qc = subcirc
    else:
        observables = ["Z"] * qc.num_qubits
        observables = "".join(observables)

    qc.measure_all()

    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    isa_circuit = pm.run(qc)

    sampler = Sampler(mode=batch)
    job = sampler.run([isa_circuit])
    logger.info(f"job id: {job.job_id()}")
    job_result = job.result()
    logger.info(job_result)
    pub_result = job_result[0].data.meas.get_counts()
    logger.info(pub_result)

    expectation_value = 0
    for key, value in pub_result.items():
        contribution = 1
        # Qiskit get_counts() keys are "<measure-all reg> <other regs>" with the
        # measure-all register first and little-endian (qubit 0 is the last bit);
        # take that register and reverse it to align with the qubit-0-first
        # observable string.
        for bit, obs in zip(reversed(key.split()[0]), observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    # print(expectation_value)
    return expectation_value


@constraint(
    processors=[
        {"processorType": "CPU", "computingUnits": "1"},
        {"processorType": "GPU", "computingUnits": "1"},
    ]
)
@task(returns=1)
def _expec_value_qiskit_gpu(qc, shots=1024, method="automatic"):
    if qc is None:
        return None

    if type(qc) != qiskit.circuit.quantumcircuit.QuantumCircuit:
        observables = qc[1]
        subcirc = qc[0]
        qc = subcirc
    else:
        observables = ["Z"] * qc.num_qubits
        observables = "".join(observables)

    qc.measure_all()

    import subprocess

    try:
        subprocess.check_output("nvidia-smi")
        logger.info("Nvidia GPU detected!")
    except Exception:
        logger.info("No Nvidia GPU in system!")

    simulator = AerSimulator(device="GPU", method=method, cuStateVec_enable=True)

    circ = transpile(qc, simulator)

    job = simulator.run(circ, shots=shots)
    result = job.result()

    counts = result.get_counts()
    logger.info(counts)

    expectation_value = 0
    for key, value in counts.items():
        contribution = 1
        # Qiskit get_counts() keys are "<measure-all reg> <other regs>" with the
        # measure-all register first and little-endian (qubit 0 is the last bit);
        # take that register and reverse it to align with the qubit-0-first
        # observable string.
        for bit, obs in zip(reversed(key.split()[0]), observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")

        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    return expectation_value


def _check_ibm_qc():
    urls = {"http": "http://localhost:44433", "https": "http://localhost:44433"}

    proxies = {"urls": urls}

    token = os.environ.get("IBM_QUANTUM_TOKEN", "")
    instance = os.environ.get("IBM_QUANTUM_INSTANCE", "")
    channel = os.environ.get("IBM_QUANTUM_CHANNEL", "")
    max_time_ibm = os.environ.get("IBM_QUANTUM_MAX_TIME", 60)
    quantum_chip = os.environ.get("IBM_QUANTUM_QPU_NAME", "")

    logger.info(quantum_chip)

    if channel == "ibm_quantum":
        service = QiskitRuntimeService(channel=channel, token=token, proxies=proxies)
    elif channel == "ibm_cloud":
        service = QiskitRuntimeService(
            channel=channel, token=token, instance=instance  # , proxies=proxies
        )
    else:
        raise ValueError(f"Unsupported IBM Runtime channel: {channel!r}")

    if quantum_chip != "":
        backend = service.backend(quantum_chip)
    else:
        backend = service.least_busy(operational=True, simulator=False)

    logger.info(backend.operation_names)

    batch = Batch(backend=backend, max_time=int(max_time_ibm))
    logger.info(batch)
    return batch, backend


def gate_cutting_subcircuits(
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

    observables : str, optional
        Per-qubit observable string (e.g. ``"ZZ"``). When given, the measurement
        basis is rotated accordingly before cutting; ``"I"`` entries become identity
        observables. Defaults to the all-``Z`` observable.

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

    dag = convert_to_dag(dag, observables)
    if observables:
        dag = change_basis_dag(dag, observables)

    # Serial -> native backend circuit; PyCOMPSs -> serializable Subcircuit.
    dispatch_map = {
        "qibo": lambda d, n: make_subcircuit(d, n, "qibo"),
        "qiskit": lambda d, n: make_subcircuit(d, n, "qiskit"),
        "cudaq": lambda d, n: make_subcircuit(d, n, "cudaq"),
    }

    number_cuts = len(gates_cut)

    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [
            dag.subgraph(c).copy() for c in nx.connected_components(dag.to_undirected())
        ]
        subcircuits = []
        for s in S:
            tmp_cuts = [c for c in gates_cut if s.has_node(c)]
            if tmp_cuts:
                subcircuits += _execute_gate_cutting_subcircuits(
                    s, tmp_cuts, software=software
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
    elif gates_cut:
        subcircuits = _execute_gate_cutting_subcircuits(
            dag, gates_cut, software=software
        )
    else:
        s_new, highest_qubit = _update_qubits(dag)
        try:
            subcircuits = [dispatch_map[software](s_new, highest_qubit)]
        except KeyError:
            raise UnsupportedBackendError(
                f"Unsupported software backend: '{software}'. "
                f"Expected one of {list(dispatch_map.keys())}."
            )

    number_components = max(1, len(subcircuits) // (6**number_cuts))
    return CutSubcircuits(
        subcircuits,
        number_cuts=number_cuts,
        number_components=number_components,
        reconstructor=gate_cutting_subcircuit_reconstruction,
    )


def _execute_gate_cutting_subcircuits(dag, gates_cut, software):
    new_dag = _generate_cut(dag, gates_cut)
    _assert_identity_observable_supported(new_dag)

    # Serial -> native backend circuit; PyCOMPSs -> serializable Subcircuit.
    dispatch_map = {
        "qibo": lambda d, n: make_subcircuit(d, n, "qibo"),
        "qiskit": lambda d, n: make_subcircuit(d, n, "qiskit"),
        "cudaq": lambda d, n: make_subcircuit(d, n, "cudaq"),
    }

    subcircuits = []
    for index in range(6 ** len(gates_cut)):

        copy_graph = new_dag.copy()
        copy_graph = _remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(copy_graph.to_undirected())

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy())

        # Populates graph_components in place (one DiGraph per component).
        _generate_gate_cutting(copy_graph, gates_cut, index, graph_components)

        for i, s in enumerate(graph_components):
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


def gate_cutting_subcircuit_reconstruction(
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
    >>> reconstructed_value = gate_cutting_subcircuit_reconstruction(raw_values, number_of_cuts, number_components=2)
    >>> print(reconstructed_value)
    0.0
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
    else:
        expected_values = [
            expected_values[i : i + number_components]
            for i in range(0, len(expected_values), number_components)
        ]
        results = []
        for index, exp_values in enumerate(expected_values):
            result = _change_sign_gate_cutting(exp_values, index)
            results.append(result)

        results = compss_wait_on(results)
        final_recons = 1 / (2**number_cuts) * sum(results)
        return final_recons
