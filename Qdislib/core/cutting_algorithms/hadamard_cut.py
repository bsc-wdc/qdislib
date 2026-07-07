#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Ancilla / double-Hadamard-test gate cutting (Harrow & Lowe, arXiv:2403.01018).

An alternative to the measure-and-prepare gate cut in :mod:`gate_cutting`: instead
of measuring a cut qubit mid-circuit and re-preparing on the other fragment, this
replaces a two-qubit gate ``U`` across a cut with an ensemble of **local** circuits,
each fragment carrying **one extra ancilla qubit** that runs a Hadamard test. The
fragments are fully independent (no shared ancilla, no mid-circuit feed-forward), so
they run as independent parallel subcircuits and their results are combined
classically -- a drop-in fit for Qdislib's parallel subcircuit model, and a good
match for superconducting hardware where mid-circuit measurement is costly.

Method (paper Fig. 3, Sec. 3.4). Write ``U = sum_i c_i V_i (x) W_i`` with ``c_i > 0``
and ``V_i, W_i`` local unitaries (an operator-Schmidt / Cartan decomposition). Per
sample, draw settings ``(i, j, g)`` from the PMF ``p(i,j,g) = c_i c_j / phi`` (zero if
``i == j and g == 1``), where ``phi = 2||c||_1^2 - ||c||_2^2`` is the 1-norm (the
sampling overhead; ``gamma(U) = 2||c||_1^2 - 1`` for an operator-Schmidt cut, e.g.
``gamma(CZ) = gamma(CNOT) = 3``, optimal). Each fragment runs
``|0> -H- G -ctrl(V_i/V_j)- H -measure b`` where ``G = I`` if ``g == 0`` else ``S``;
the unbiased estimator of ``<X>`` is ``phi * (-1)^(g + b1 + b2) * y``.

This module implements the **decomposition core** ``U -> [(c_i, V_i, W_i)]`` and the
``gamma``/``phi`` factors (milestone 1). The double-Hadamard-test subcircuit
construction, the sampling driver, and the reconstruction build on top of it.
"""

from __future__ import annotations

import cmath
import math

import numpy as np

__all__ = [
    "SUPPORTED_HADAMARD_GATES",
    "local_unitary_decomposition",
    "gamma_factor",
    "phi_factor",
    "gate_cutting_hadamard",
]

import itertools

# Single-qubit Paulis (the building blocks of the local unitaries below).
_I = np.eye(2, dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Two-qubit gates whose only non-local part is a ZZ interaction, so they cut into
# exactly two unitary terms with the optimal gamma. CX/CNOT reduce to CZ by
# conjugating the target with H (a local gate absorbed into the fragment).
SUPPORTED_HADAMARD_GATES = ("CZ", "CP", "CPHASE", "CU1", "RZZ", "CX", "CNOT")


def _diag_angles(diag):
    """The ``(II, IZ, ZI, ZZ)`` generator angles ``(alpha, beta, gamma, delta)`` of a
    diagonal 2-qubit unitary ``diag(e^{i phi_ab})`` -- i.e.
    ``D = e^{i alpha}(e^{i gamma Z}(x)e^{i beta Z}).exp(i delta Z(x)Z)``."""
    p00, p01, p10, p11 = np.angle(np.asarray(diag, dtype=complex))
    return (
        (p00 + p01 + p10 + p11) / 4,  # alpha: II (global phase)
        (p00 - p01 + p10 - p11) / 4,  # beta:  IZ (Z on the 2nd qubit)
        (p00 + p01 - p10 - p11) / 4,  # gamma: ZI (Z on the 1st qubit)
        (p00 - p01 - p10 + p11) / 4,  # delta: ZZ (the non-local part)
    )


def _diagonal_zz_decomposition(diag):
    """Decompose a *diagonal* 2-qubit unitary into local-unitary terms.

    A diagonal 2-qubit unitary ``D = diag(e^{i phi_ab})`` has generator in the
    ``{II, IZ, ZI, ZZ}`` basis; only the ``ZZ`` part is non-local. Writing
    ``D = e^{i alpha}(e^{i g Z} (x) e^{i b Z}) . exp(i d ZZ)`` and
    ``exp(i d ZZ) = cos d (I(x)I) + i sin d (Z(x)Z)`` gives two terms whose local
    single-qubit factors (the ``Z`` rotations and the global phase) are absorbed
    into ``V``/``W``. Returns ``[(c1, V1, W1), (c2, V2, W2)]`` with ``c_k > 0`` and
    ``V_k, W_k`` unitary, dropping a term whose coefficient is ~0.
    """
    alpha, beta, gamma, delta = _diag_angles(diag)
    rz_a = np.diag([cmath.exp(1j * gamma), cmath.exp(-1j * gamma)])  # e^{i gamma Z}
    rz_b = np.diag([cmath.exp(1j * beta), cmath.exp(-1j * beta)])  # e^{i beta Z}
    phase = cmath.exp(1j * alpha)

    terms = []
    c_cos, c_sin = abs(math.cos(delta)), abs(math.sin(delta))
    if c_cos > 1e-12:  # cos d * (I (x) I), with local rotations + global phase on V
        s = 1.0 if math.cos(delta) >= 0 else -1.0
        terms.append((c_cos, s * phase * (rz_a @ _I), rz_b @ _I))
    if c_sin > 1e-12:  # i sin d * (Z (x) Z)
        s = 1.0 if math.sin(delta) >= 0 else -1.0
        terms.append((c_sin, s * 1j * phase * (rz_a @ _Z), rz_b @ _Z))
    return terms


def local_unitary_decomposition(gate_name, params=()):
    """``U -> [(c_i, V_i, W_i)]``: the local-unitary decomposition of a cut gate.

    Supports the two-qubit gates Qdislib cuts whose non-local part is a single ``ZZ``
    interaction (:data:`SUPPORTED_HADAMARD_GATES`), giving two unitary terms with the
    optimal 1-norm. ``V_i`` acts on the control-side fragment qubit, ``W_i`` on the
    target-side. ``CX``/``CNOT`` reduce to ``CZ`` (the ``H`` on the target is local
    and applied inside the fragment by the caller). ``sum_i c_i V_i (x) W_i`` equals
    the gate up to a global phase.

    :param gate_name: Canonical gate name (upper-cased).
    :param params: Gate parameters (e.g. ``(theta,)`` for ``CP``/``RZZ``).
    :raises NotImplementedError: for a gate outside :data:`SUPPORTED_HADAMARD_GATES`.
    """
    name = str(gate_name).upper()
    theta = float(params[0]) if params else 0.0
    if name == "CZ" or (name in ("CX", "CNOT")):
        # CZ = diag(1, 1, 1, -1). CX is cut as CZ with H-conjugation on the target
        # (handled by the fragment builder); the decomposition is the CZ one.
        diag = [1, 1, 1, -1]
    elif name in ("CP", "CPHASE", "CU1"):  # CU1 is qibo's controlled-phase name
        diag = [1, 1, 1, cmath.exp(1j * theta)]
    elif name == "RZZ":
        e = cmath.exp(-1j * theta / 2)
        diag = [e, e.conjugate(), e.conjugate(), e]
    else:
        raise NotImplementedError(
            f"Hadamard-test cut not implemented for {name!r}; "
            f"supported: {SUPPORTED_HADAMARD_GATES}"
        )
    return _diagonal_zz_decomposition(diag)


def phi_factor(terms):
    """The sampling-overhead 1-norm ``phi = 2||c||_1^2 - ||c||_2^2`` of a decomposition."""
    c = np.array([t[0] for t in terms], dtype=float)
    return float(2 * c.sum() ** 2 - (c**2).sum())


def gamma_factor(terms):
    """The optimal gamma ``= 2||c||_1^2 - 1`` (holds for an operator-Schmidt cut)."""
    c = np.array([t[0] for t in terms], dtype=float)
    return float(2 * c.sum() ** 2 - 1)


# --------------------------------------------------------------------------- #
# Double-Hadamard-test cut: fragment construction + reconstruction
# --------------------------------------------------------------------------- #


def _generate_cut_hadamard(dag, gates_cut):
    """Insert SUBS1_k (control side) / SUBS2_k (target side) placeholder nodes for
    each cut gate -- like ``gate_cutting._generate_cut`` but for the Hadamard gate
    set (CZ/CP/CU1/RZZ and CX/CNOT), without the CZ/CX-only restriction. A CX/CNOT
    is cut as ``H CZ H`` (Hadamard conjugation on the target); the other supported
    gates are diagonal and need no conjugation. The SUBS nodes are later replaced by
    the double-Hadamard-test gadget in :func:`_component_hadamard_circuit`.
    """
    from Qdislib.core.find_cut.find_cut import _find_nodes_with_qubit
    from Qdislib.utils.graph_qibo import _remove_red_edges

    dag_copy = _remove_red_edges(dag.copy())
    for index, gate_name in enumerate(gates_cut, start=1):
        is_cx = (dag_copy.nodes[gate_name].get("gate") or "").lower() in ("cx", "cnot")
        q = dag_copy.nodes[gate_name]["qubits"]
        pred0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=q[0], direction="predecessor"
        )
        pred1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=q[1], direction="predecessor"
        )
        succ0 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=q[0], direction="successor"
        )
        succ1 = _find_nodes_with_qubit(
            dag_copy, gate_name, qubit=q[1], direction="successor"
        )
        dag_copy.remove_node(gate_name)
        dag_copy.add_node(f"SUBS1_{index}", gate="SUBS", qubits=(q[0],), parameters=())
        dag_copy.add_node(f"SUBS2_{index}", gate="SUBS", qubits=(q[1],), parameters=())
        if pred0:
            dag_copy.add_edge(pred0[0], f"SUBS1_{index}", color="blue")
        if succ0:
            dag_copy.add_edge(f"SUBS1_{index}", succ0[0], color="blue")
        if is_cx:  # H . SUBS2 . H on the target (CX = H CZ H)
            hpre, hpost = f"HPRE_{index}", f"HPOST_{index}"
            dag_copy.add_node(hpre, gate="h", qubits=(q[1],), parameters=())
            dag_copy.add_node(hpost, gate="h", qubits=(q[1],), parameters=())
            if pred1:
                dag_copy.add_edge(pred1[0], hpre, color="blue")
            dag_copy.add_edge(hpre, f"SUBS2_{index}", color="blue")
            dag_copy.add_edge(f"SUBS2_{index}", hpost, color="blue")
            if succ1:
                dag_copy.add_edge(hpost, succ1[0], color="blue")
        else:
            if pred1:
                dag_copy.add_edge(pred1[0], f"SUBS2_{index}", color="blue")
            if succ1:
                dag_copy.add_edge(f"SUBS2_{index}", succ1[0], color="blue")
    return dag_copy


def _cut_settings(num_terms):
    """The valid ``(i, j, g)`` settings: all i,j in [m], g in {0,1}, minus i==j,g==1.

    Each carries weight ``c_i c_j (-1)^g`` in the reconstruction (the PMF is this
    over ``phi``); ``g`` selects the Hadamard test that reads the real (g=0) or
    imaginary (g=1) part.
    """
    settings = []
    for i in range(num_terms):
        for j in range(num_terms):
            for g in (0, 1):
                if i == j and g == 1:
                    continue
                settings.append((i, j, g))
    return settings


def _subs_gadget(circuit, gates, qubit, ancilla, u_i, u_j, g):
    """Append the double-Hadamard-test gadget for one cut side to ``circuit``.

    Implements the select gate ``|0><0| (x) U_i + |1><1| (x) U_j`` (``ancilla`` the
    control), wrapped in the Hadamard test ``H . G(g) . [select] . H``. ``G = I`` if
    ``g == 0`` else ``S`` (real- vs imaginary-part test). ``U_i, U_j`` are diagonal
    for the supported gates, so the select gate is a diagonal 2-qubit unitary and is
    emitted in **standard gates** (RZ + CX for its ZZ part) -- no arbitrary ``Unitary``
    matrix -- so the fragment is expressible in the neutral DAG and runnable on any
    backend. The ancilla is read out as a Z expectation by the caller.
    """
    # Select gate diagonal on (ancilla, qubit): anc=0 -> U_i, anc=1 -> U_j.
    gsel = [u_i[0, 0], u_i[1, 1], u_j[0, 0], u_j[1, 1]]
    _, beta, gamma, delta = _diag_angles(gsel)  # global phase (alpha) is irrelevant
    circuit.add(gates.H(ancilla))
    if g == 1:
        circuit.add(gates.S(ancilla))
    circuit.add(gates.RZ(ancilla, -2 * gamma))  # e^{i gamma Z} on the control
    circuit.add(gates.RZ(qubit, -2 * beta))  # e^{i beta Z} on the target
    circuit.add(gates.CNOT(ancilla, qubit))  # exp(i delta ZZ) = CX . RZ . CX
    circuit.add(gates.RZ(qubit, -2 * delta))
    circuit.add(gates.CNOT(ancilla, qubit))
    circuit.add(gates.H(ancilla))


def _component_hadamard_circuit(component, assignments):
    """Build one fragment's Hadamard-test circuit from its DAG component.

    Walks the component in topological order, renumbering its qubits to ``0..k-1``;
    each ``SUBS`` node is replaced by :func:`_subs_gadget` on a fresh ancilla
    (indices ``k, k+1, ...``). ``assignments`` maps a SUBS node id to its
    ``(u_i, u_j, g)`` (already the V- or W-side unitaries). Returns ``(circuit,
    qubit_map, ancilla_indices)`` where ``qubit_map`` is original-qubit -> new-index.
    """
    import networkx as nx
    from qibo import gates, models

    from Qdislib.utils.graph_qibo import _resolve_qibo_gate

    orig_qubits = sorted(
        {q for _, d in component.nodes(data=True) for q in (d.get("qubits") or ())}
    )
    qubit_map = {q: i for i, q in enumerate(orig_qubits)}
    n = len(orig_qubits)
    num_subs = sum(1 for nid in component.nodes if str(nid).startswith("SUBS"))
    circuit = models.Circuit(n + num_subs)
    ancillas = []
    next_anc = n
    for node in nx.topological_sort(component):
        data = component.nodes[node]
        name = (data.get("gate") or "").upper()
        qubits = tuple(qubit_map[q] for q in (data.get("qubits") or ()))
        if str(node).startswith("SUBS"):
            u_i, u_j, g = assignments[node]
            _subs_gadget(circuit, gates, qubits[0], next_anc, u_i, u_j, g)
            ancillas.append(next_anc)
            next_anc += 1
        elif name in ("OBSERVABLE I", "MEASURE"):
            continue
        else:
            gate_cls = _resolve_qibo_gate(data.get("gate") or "")
            params = data.get("parameters")
            if params:
                circuit.add(gate_cls(*qubits, *params))
            else:
                circuit.add(gate_cls(*qubits))
    return circuit, qubit_map, ancillas


def _expectation_z(circuit, z_qubits):
    """``<prod_q Z_q>`` for the fragment state -- a diagonal operator, so a cheap
    signed sum over ``|amplitude|^2`` (qibo big-endian: qubit q is bit ``n-1-q``)."""
    state = np.asarray(circuit().state())
    probs = np.abs(state) ** 2
    n = circuit.nqubits
    masks = [1 << (n - 1 - q) for q in z_qubits]
    total = 0.0
    for x, p in enumerate(probs):
        if p == 0.0:
            continue
        sign = 1
        for m in masks:
            if x & m:
                sign = -sign
        total += p * sign
    return float(total)


def _expectation_z_sampled(circuit, z_qubits, shots, seed=None):
    """Shot-based estimate of ``<prod_q Z_q>`` for the fragment (measures all qubits,
    then averages the parity over ``z_qubits`` across the shot frequencies)."""
    from qibo import gates

    circ = circuit.copy(deep=True)
    circ.add(gates.M(*range(circ.nqubits)))
    freqs = circ(nshots=shots).frequencies(binary=True)  # bitstring -> count
    total = 0.0
    for bitstring, count in freqs.items():
        sign = 1
        for q in z_qubits:  # bitstring[q] is qubit q (M added in 0..n-1 order)
            if bitstring[q] == "1":
                sign = -sign
        total += sign * count
    return total / shots


def gate_cutting_hadamard(
    circuit, gates_cut, observable, software="qibo", shots=None, seed=None
):
    """Estimate ``<observable>`` via the ancilla / double-Hadamard-test gate cut.

    Cuts the two-qubit ``gates_cut`` (names like ``"CZ_2"``) with the Harrow-Lowe
    space-like cut (arXiv:2403.01018): each fragment carries one local ancilla
    running a Hadamard test, the fragments are independent, and the result is
    combined classically. This is the **exact** (statevector, full-enumeration)
    reconstruction -- the correctness anchor for the sampling variant.

    :param circuit: A qibo circuit.
    :param gates_cut: List of cut gate names (each a two-qubit gate in
        :data:`SUPPORTED_HADAMARD_GATES`).
    :param observable: A Z/I Pauli string over the circuit's qubits.
    :param software: Only ``"qibo"`` today.
    :param shots: If given, estimate each fragment value from ``shots`` measurement
        samples (unbiased, with shot noise); if ``None``, use the exact state vector.
    :param seed: Optional RNG seed for the sampled path.
    :return: The reconstructed expectation value.
    """
    if software != "qibo":
        raise NotImplementedError("gate_cutting_hadamard supports software='qibo'")
    import networkx as nx

    from Qdislib.utils.backends import convert_to_dag

    dag = convert_to_dag(circuit)
    # Per-cut local-unitary decompositions [(c, V, W), ...] and the settings.
    decomps, settings_per_cut = [], []
    for name in gates_cut:
        gate = (dag.nodes[name].get("gate") or "").upper()
        params = dag.nodes[name].get("parameters") or ()
        terms = local_unitary_decomposition(gate, params)
        decomps.append(terms)
        settings_per_cut.append(_cut_settings(len(terms)))

    cut_dag = _generate_cut_hadamard(dag, gates_cut)  # SUBS1_k / SUBS2_k (+ CX H-conj)
    num_qubits = circuit.nqubits

    total = 0.0
    for config in itertools.product(*settings_per_cut):  # one (i,j,g) per cut
        weight = 1.0
        assignments = {}
        for k, (i, j, g) in enumerate(config, start=1):
            c_i, v_i, w_i = decomps[k - 1][i]
            c_j, v_j, w_j = decomps[k - 1][j]
            weight *= c_i * c_j * (-1) ** g
            assignments[f"SUBS1_{k}"] = (v_i, v_j, g)  # control side -> V
            assignments[f"SUBS2_{k}"] = (w_i, w_j, g)  # target side -> W
        if weight == 0.0:
            continue
        # Split into independent fragments and multiply their signed Z-expectations.
        undirected = _remove_red_free(cut_dag).to_undirected()
        product = 1.0
        for comp_nodes in nx.connected_components(undirected):
            comp = cut_dag.subgraph(comp_nodes).copy()
            frag, qmap, ancillas = _component_hadamard_circuit(comp, assignments)
            z_qubits = list(ancillas)
            for q, ch in enumerate(observable):
                if ch == "Z" and q in qmap:
                    z_qubits.append(qmap[q])
            if shots is None:
                product *= _expectation_z(frag, z_qubits)
            else:
                product *= _expectation_z_sampled(frag, z_qubits, shots, seed)
        total += weight * product
    return total


def _remove_red_free(dag):
    """Drop the red (cut-dependency) edges, leaving the blue circuit-flow edges that
    define the independent fragments (mirrors ``_remove_red_edges``)."""
    from Qdislib.utils.graph_qibo import _remove_red_edges

    return _remove_red_edges(dag.copy())
