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

from Qdislib.core.cutting_algorithms._api import standard_arguments

logger = logging.getLogger(__name__)

__all__ = ["estimate", "cut_cost"]


@standard_arguments
def estimate(
    circuit,
    observable=None,
    *,
    max_qubits=None,
    backend=None,
    shots=1024,
    seed=None,
    cache=None,
    cache_mode=None,
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

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param observable: Observable spec -- a Pauli string, a ``{pauli: coeff}``
        mapping, or ``[(coeff, pauli)]`` weighted sum. Defaults to all-``Z``.
    :param max_qubits: Maximum qubits per subcircuit (guides the cut search).
    :param shots: Measurement shots per subcircuit; ``None`` is exact.
    :param backend: Where the subcircuits run (see :func:`~Qdislib.gate_cutting`).
    :param seed: Seed for the simulators' sampling.
    :param cache: Reuse fragment circuits already run: what ``cache=`` takes on
        :func:`~Qdislib.circuit_cutting`.
    :param cache_mode: Deprecated and ignored, as on :func:`~Qdislib.gate_cutting`
        (a custom :class:`~Qdislib.utils.circuit_cache.CacheStrategy` still runs the
        1.x implementation).
    :param gate_cut: Allow :func:`find_cut` to consider gate cuts.
    :param wire_cut: Allow :func:`find_cut` to consider wire cuts.
    :param return_metrics: If ``True`` return ``(value, metrics)`` where
        ``metrics`` is ``{"cut", "cut_type", "num_cuts", "find_cut_seconds",
        "run_seconds", "cache"}``: the chosen cut, ``"gate"`` or ``"wire"``, the number
        of cuts, the timings and the store's :meth:`stats` (``None`` without a cache).
    :raises ValueError: For an invalid observable, or limits ``find_cut`` refuses.
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
            observable,
            backend=backend,
            shots=shots,
            seed=seed,
            cache=cache,
            cache_mode=cache_mode,
        )
    else:
        value = gate_cutting(
            circuit,
            cut,
            observable,
            backend=backend,
            shots=shots,
            seed=seed,
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
        "cache": _cache_stats(cache),
    }
    return value, metrics


def _cache_stats(cache):
    """The statistics of the store a call with this ``cache=`` uses, its own or the
    session's (:func:`~Qdislib.use_cache`); ``None`` when it is off."""
    from Qdislib.core.cutting_algorithms.engine_cache import resolve

    used = resolve(cache)
    return None if used is None else used.store.stats()


def cut_cost(cut, wire=None):
    """Deprecated: use :func:`~Qdislib.cutting_overhead`, which works for every cut
    type and mix and also reports the circuits to run and ``gamma``.

    Estimate the reconstruction cost of a cut *before* running it.

    Circuit cutting trades one circuit for many subcircuit evaluations: gate
    cutting expands each cut into 6 quasiprobability terms and wire cutting into
    8, so ``k`` cuts need ``base ** k`` term evaluations to reconstruct an
    expectation value (:func:`~Qdislib.gate_cutting` loops over ``range(6 ** k)``,
    :func:`~Qdislib.wire_cutting` over ``range(8 ** k)``). That factor -- not the
    subcircuit size -- dominates the cost and grows exponentially in the number
    of cuts, so this lets you size a run before launching it.

    :param cut: A cut from :func:`~Qdislib.find_cut`: a list of gate identifiers
        for a gate cut (``["CZ_3"]``) or of endpoint pairs for a wire cut
        (``[("H_1", "CZ_2")]``).
    :param wire: Force the cut kind instead of inferring it from the element shape
        (wire cuts are pairs, gate cuts are plain identifiers). Rarely needed.
    :return: A dict with ``num_cuts``, ``kind`` (``"gate"``/``"wire"``), ``base``
        (6 or 8) and ``terms`` (``base ** num_cuts`` -- the number of subcircuit
        configurations the reconstruction evaluates).

    Example:
        .. code-block:: python

            import Qdislib as qd

            cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)
            qd.cut_cost(cut)   # {'num_cuts': 2, 'kind': 'gate', 'base': 6, 'terms': 36}
    """
    import warnings

    warnings.warn(
        "cut_cost() is deprecated: use cutting_overhead(circuit, cuts); it is "
        "removed in the next release.",
        DeprecationWarning,
        stacklevel=2,
    )
    cut = list(cut or [])
    if wire is None:
        # Same shape test as estimate(): wire cuts are (source, target) pairs,
        # gate cuts are plain gate identifiers.
        wire = bool(cut) and isinstance(cut[0], (tuple, list))
    base = 8 if wire else 6
    num_cuts = len(cut)
    return {
        "num_cuts": num_cuts,
        "kind": "wire" if wire else "gate",
        "base": base,
        "terms": base**num_cuts,
    }
