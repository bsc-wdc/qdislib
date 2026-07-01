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

"""High-level one-call estimation: cut, run and reconstruct in a single step."""

import logging
import time

logger = logging.getLogger(__name__)

__all__ = ["estimate"]


def estimate(
    circuit,
    observables=None,
    max_qubits=None,
    shots=1024,
    software="qibo",
    cache=None,
    cache_mode="expectation",
    gate_cut=True,
    wire_cut=True,
    return_metrics=False,
):
    """Automatically cut, execute and reconstruct ``<O>`` for a circuit.

    The whole pipeline in one call: choose a cut with :func:`find_cut` (honouring
    ``max_qubits``), then gate- or wire-cut and reconstruct -- whichever
    :func:`find_cut` selected. For the step-by-step control use the
    ``*_subcircuits`` / :func:`expectation_value` / ``reconstruct`` workflow
    instead.

    :param circuit: Circuit to estimate (Qibo/Qiskit/CUDA-Q, or a DAG).
    :param observables: Observable spec -- a Pauli string, a ``{pauli: coeff}``
        mapping, or ``[(coeff, pauli)]`` weighted sum. Defaults to all-``Z``.
    :param max_qubits: Maximum qubits per subcircuit (guides the cut search).
    :param shots: Measurement shots per subcircuit.
    :param software: Backend used to execute the subcircuits.
    :param cache: Optional :class:`~Qdislib.utils.circuit_cache.CircuitCache` to
        reuse equivalent subcircuits (see
        :func:`~Qdislib.utils.circuit_cache.default_circuit_cache`).
    :param cache_mode: Cache granularity (``"expectation"``/``"counts"``/strategy).
    :param gate_cut: Allow :func:`find_cut` to consider gate cuts.
    :param wire_cut: Allow :func:`find_cut` to consider wire cuts.
    :param return_metrics: If ``True`` return ``(value, metrics)`` where
        ``metrics`` reports the chosen cut, its type, the number of cuts, the
        find-cut and run timings, and the cache statistics.
    :return: The reconstructed expectation value, or ``(value, metrics)``.

    Example:
        .. code-block:: python

            import Qdislib as qd

            value = qd.estimate(circuit, "ZZZZ", max_qubits=3)
            value, metrics = qd.estimate(circuit, max_qubits=3, return_metrics=True)
    """
    from Qdislib.core.find_cut.find_cut import find_cut
    from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
    from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

    start = time.perf_counter()
    cut = find_cut(circuit, max_qubits=max_qubits, gate_cut=gate_cut, wire_cut=wire_cut)
    find_cut_seconds = time.perf_counter() - start

    # find_cut returns gate identifiers (e.g. ``["CZ_3"]``) for a gate cut, or
    # endpoint pairs (e.g. ``[("H_1", "CZ_2")]``) for a wire cut; dispatch on the
    # element shape.
    is_wire_cut = bool(cut) and isinstance(cut[0], (tuple, list))

    start = time.perf_counter()
    if is_wire_cut:
        value = wire_cutting(
            circuit,
            cut,
            observables=observables,
            shots=shots,
            software=software,
            cache=cache,
            cache_mode=cache_mode,
        )
    else:
        value = gate_cutting(
            circuit,
            cut,
            observables=observables,
            shots=shots,
            software=software,
            cache=cache,
            cache_mode=cache_mode,
        )
    run_seconds = time.perf_counter() - start

    if not return_metrics:
        return value

    metrics = {
        "cut": cut,
        "cut_type": "wire" if is_wire_cut else "gate",
        "num_cuts": len(cut),
        "find_cut_seconds": find_cut_seconds,
        "run_seconds": run_seconds,
        "cache": cache.stats() if cache is not None else None,
    }
    return value, metrics
