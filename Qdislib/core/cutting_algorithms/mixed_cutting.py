#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Mixed wire + gate cutting on a single circuit (circuit knitting).

Cut some two-qubit gates (the measure-and-prepare gate cut) *and* some wires (the
wire cut) on the same circuit, in one integrated, fully parallel workflow. The two
quasiprobability decompositions are independent, so they **tensor**: the combined QPD
runs over the product of the two configuration spaces and the reconstruction is

    <O> = 1 / 2^(n_gate + n_wire)
          * sum_{gate cfg} sum_{wire cfg}
                gate_sign(cfg) * wire_sign(cfg) * prod_fragment <fragment>

Every fragment across the whole product space is independent, so they all fan out as
PyCOMPSs tasks and fold back through one reduction (single fan-in, or a balanced tree
when ``reduce_chunk`` is set) -- the same parallel shape as the standalone cutters.

The heavy lifting is reused unchanged: the gate-cut DAG surgery + base-6
measure-and-prepare substitution (:func:`_generate_cut`, :func:`_apply_gate_substitution`),
the wire-cut DAG surgery + base-8 basis substitution (:func:`_insert_wire_cut_markers`,
:func:`_apply_wire_substitution`), and the fragment builder/executor
(:func:`_build_and_run_gate_subcircuit`, which also handles the gate cut's
collapse-and-continue measures a mixed fragment can carry).
"""

from __future__ import annotations

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN

__all__ = ["gate_wire_cutting"]


def _config_sign(index, base, flip_digits):
    """The +/-1 quasiprobability sign of configuration ``index``: flip once per digit
    (in base ``base``) that lands in ``flip_digits``. Gate cut uses base 6 / {3, 5};
    wire cut uses base 8 / {3, 5, 7} -- matching the standalone ``_change_sign`` rules.
    """
    flip = False
    n = index
    while n:
        if n % base in flip_digits:
            flip = not flip
        n //= base
    return -1 if flip else 1


@task(returns=1, values=COLLECTION_IN)
def _signed_product(values, sign):
    """One configuration's contribution: ``sign * prod(fragment values)`` (the fan-in
    task; ``values`` is the ``COLLECTION_IN`` of that config's fragment futures)."""
    product = 1.0
    for value in values:
        if value is not None:
            product *= value
    return sign * product


def gate_wire_cutting(
    circuit,
    gates_cut,
    wire_cut,
    observable=None,
    software="qibo",
    shots=1024,
    method="automatic",
    sync=True,
    qpu=False,
    qpu_dict=None,
):
    """Reconstruct ``<observable>`` from a circuit cut by BOTH gate and wire cuts.

    :param circuit: A qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: Two-qubit gate node names to cut (CZ/CX), as in
        :func:`gate_cutting`.
    :param wire_cut: Wire cuts as gate-pair edges ``[(source, target), ...]``, as in
        :func:`wire_cutting` (each is expanded to one cut per shared qubit).
    :param observable: A Z/I (or general Pauli) string over the circuit's qubits.
    :param software: Fragment-execution backend (``"qibo"``/``"qiskit"``/``"cudaq"``/
        ``"pennylane"``).
    :param shots: Samples per fragment expectation.
    :param method: Backend evaluation method (forwarded to the gate-cut executor).
    :param sync: If ``True`` block and return the value; if ``False`` return the list
        of per-configuration futures (advanced overlap of several cuts).
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :return: The reconstructed expectation value.
    """
    import networkx as nx

    from Qdislib.utils.backends import convert_to_dag, change_basis_dag
    from Qdislib.utils.graph_qibo import _max_qubits_graph
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _generate_cut,
        _apply_gate_substitution,
        _build_and_run_gate_subcircuit,
        _expand_general_gate_cuts,
    )
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expand_wire_cuts,
        _insert_wire_cut_markers,
        _apply_wire_substitution,
        _remove_red_edges,
    )

    dag = convert_to_dag(circuit, observable)
    if observable:
        dag = change_basis_dag(dag, observable)
    num_qubits = _max_qubits_graph(dag)  # original qubit count (before wire additions)

    gates_cut = list(gates_cut or [])
    # A non-CZ/CX cut gate is decomposed into CZ + single-qubit gates and its CZs cut.
    if gates_cut:
        dag, gates_cut = _expand_general_gate_cuts(dag, gates_cut)
    # Gate surgery: insert the SUBS1/SUBS2 markers (also strips red edges internally).
    dag = _generate_cut(dag, gates_cut) if gates_cut else _remove_red_edges(dag.copy())
    # Wire surgery: expand gate-pair edges to single-wire cuts + insert O_i/PS_i markers.
    wire_cuts = _expand_wire_cuts(dag, wire_cut) if wire_cut else []
    if wire_cuts:
        _insert_wire_cut_markers(dag, wire_cuts, num_qubits)

    n_gate, n_wire = len(gates_cut), len(wire_cuts)
    total_qubits = num_qubits + n_wire  # fresh wire lines are numbered above the rest

    # Enumerate the PRODUCT of the two QPD configuration spaces; each configuration
    # applies both substitutions to a fresh DAG copy, splits into independent fragments,
    # fans them out, and contributes sign * prod(fragments).
    results = []
    for gate_idx in range(6**n_gate):
        for wire_idx in range(8**n_wire):
            config_dag = _remove_red_edges(dag.copy())
            if n_gate:
                _apply_gate_substitution(config_dag, gates_cut, gate_idx)
            if n_wire:
                _apply_wire_substitution(config_dag, total_qubits, wire_idx, wire_cuts)
            config_dag = _remove_red_edges(config_dag)

            fragment_values = [
                _build_and_run_gate_subcircuit(
                    config_dag.subgraph(component).copy(),
                    total_qubits,
                    software,
                    shots,
                    method,
                    False,  # gpu
                    0,  # gpu_min_qubits
                    qpu,
                    qpu_dict,
                    None,  # batch
                    None,  # ibm_backend
                )
                for component in nx.connected_components(config_dag.to_undirected())
            ]
            sign = _config_sign(gate_idx, 6, {3, 5}) * _config_sign(
                wire_idx, 8, {3, 5, 7}
            )
            results.append(_signed_product(fragment_values, sign))

    if not sync:
        return results
    results = compss_wait_on(results)
    return 1 / (2 ** (n_gate + n_wire)) * sum(results)
