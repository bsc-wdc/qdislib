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

The two-qubit gates whose only non-local part is a ``ZZ`` interaction (``CZ``, ``CX``,
``CP``/``CU1``, ``RZZ``) take a cheap closed-form decomposition with the optimal
1-norm; **any other two-qubit gate** (``iSWAP``, ``SWAP``, ``fSim``, an arbitrary
``U``) is decomposed via the general KAK / Cartan (Weyl) route, which writes the gate
as ``sum_P alpha_P (P (x) P)`` sandwiched by single-qubit gates. Both feed the same
double-Hadamard-test fragment builder, which emits only standard named gates
(``H``/``S``/``U3``/``CU3``/``U1``/``CNOT``) so the fragments stay expressible in the
neutral DAG and runnable on any backend.
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
    "gate_cutting_hadamard_subcircuits",
    "gate_cutting_hadamard_subcircuit_reconstruction",
    "gate_cutting_hadamard_sampling",
    "gate_cutting_hadamard_sampling_subcircuits",
    "HadamardCutSubcircuits",
]

import itertools

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN
from Qdislib.utils.graph_qibo import circuit_qibo_to_dag

# Single-qubit Paulis (the building blocks of the local unitaries below).
_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULI = {"I": _I, "X": _X, "Y": _Y, "Z": _Z}

# Two-qubit gates whose only non-local part is a ZZ interaction: they cut into two
# unitary terms with the optimal gamma via a cheap closed form (no KAK). CX/CNOT
# reduce to CZ by conjugating the target with H (a local gate absorbed into the
# fragment). Every *other* two-qubit gate is still cuttable -- via the general KAK /
# Cartan path (:func:`_general_local_terms`) -- just without this shortcut.
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


def _u3_angles(u):
    """Decompose a single-qubit unitary ``u`` (any ``U(2)``) into qibo ``U3`` angles
    plus a global phase: returns ``(theta, phi, lam, phase)`` with
    ``u = e^{i phase} * U3(theta, phi, lam)`` and
    ``U3(theta, phi, lam) = RZ(phi) . RY(theta) . RZ(lam)`` (qibo's convention, an
    SU(2) matrix whose entries are ``e^{-iA}c, -e^{-iB}s, e^{iB}s, e^{iA}c`` with
    ``c = cos(theta/2)``, ``s = sin(theta/2)``, ``A = (phi+lam)/2``, ``B = (phi-lam)/2``).

    The half-angles ``A``, ``B`` are read straight off the entry phases -- *not* by
    halving a phase *difference*, which a 2*pi wrap would turn into a spurious pi error
    and corrupt the SU(2) structure. Handles the ``theta = 0`` / ``theta = pi``
    degeneracies (then only one of ``A``, ``B`` is defined; the other is free)."""
    u = np.asarray(u, dtype=complex)
    # Direct 2x2 determinant: numpy's LU-based det spuriously warns "divide by zero"
    # on matrices with exact off-diagonal zeros (e.g. diag(1, i)).
    det = u[0, 0] * u[1, 1] - u[0, 1] * u[1, 0]
    phase = cmath.phase(det) / 2  # det(U3) == 1, so det(u) == e^{2i phase}
    su = cmath.exp(-1j * phase) * u  # now in SU(2), entries as documented above
    theta = 2 * math.atan2(abs(su[1, 0]), abs(su[0, 0]))
    a = cmath.phase(su[1, 1]) if abs(su[1, 1]) > 1e-12 else 0.0  # A = (phi+lam)/2
    b = cmath.phase(su[1, 0]) if abs(su[1, 0]) > 1e-12 else 0.0  # B = (phi-lam)/2
    return theta, a + b, a - b, phase


def _cartan_pauli_coefficients(a, b, c):
    """The ``{II, XX, YY, ZZ}`` coefficients of the canonical non-local gate
    ``N(a,b,c) = exp(i(a XX + b YY + c ZZ)) = sum_P alpha_P (P (x) P)``. Since ``XX``,
    ``YY``, ``ZZ`` commute, ``N`` factorises into three commuting rotations whose
    product expands onto the four commuting Pauli pairs."""
    ca, sa = math.cos(a), math.sin(a)
    cb, sb = math.cos(b), math.sin(b)
    cc, sc = math.cos(c), math.sin(c)
    return {
        "I": ca * cb * cc + 1j * sa * sb * sc,
        "X": ca * sb * sc + 1j * sa * cb * cc,
        "Y": sa * cb * sc + 1j * ca * sb * cc,
        "Z": sa * sb * cc + 1j * ca * cb * sc,
    }


def _terms_from_factors(coeffs, k_left_a, k_left_b, k_right_a, k_right_b):
    """Assemble ``[(c_P, V_P, W_P)]`` from Cartan Pauli coefficients and the KAK local
    gates: ``U = sum_P coeffs[P] (k_left_a P k_right_a) (x) (k_left_b P k_right_b)``.
    Folds each coefficient's phase into ``V_P`` (so ``c_P = |coeff| > 0``) and drops
    near-zero terms. ``V`` acts on the control-side qubit, ``W`` on the target side."""
    terms = []
    for name, coeff in coeffs.items():
        if abs(coeff) < 1e-12:
            continue
        p = _PAULI[name]
        v = k_left_a @ p @ k_right_a
        w = k_left_b @ p @ k_right_b
        terms.append((abs(coeff), (coeff / abs(coeff)) * v, w))
    return terms


def _kak_local_terms(u):
    """General ``U -> [(c_i, V_i, W_i)]`` via the KAK / Cartan (Weyl) decomposition.

    Any two-qubit unitary is ``U = e^{i phi}(k1l (x) k1r) N(a,b,c) (k2l (x) k2r)`` with
    single-qubit ``k`` and ``N(a,b,c) = sum_P alpha_P (P (x) P)`` (four unitary Pauli
    pairs). This yields the operator-Schmidt / local-unitary decomposition the
    double-Hadamard-test cut needs, with the optimal 1-norm (e.g. ``gamma(iSWAP) =
    gamma(SWAP) = 7``, ``gamma(CZ) = 3``). Uses qiskit's ``TwoQubitWeylDecomposition``
    as the (well-tested) Cartan engine; qiskit's qubit ordering is the reverse of
    qibo's, so both control/target assignments are tried and the one that reconstructs
    ``U`` (to 1e-9) is returned.

    :raises RuntimeError: if neither assignment reconstructs ``U`` (caller falls back).
    """
    from qiskit.synthesis import TwoQubitWeylDecomposition

    u = np.asarray(u, dtype=complex)
    weyl = TwoQubitWeylDecomposition(u)
    coeffs = _cartan_pauli_coefficients(weyl.a, weyl.b, weyl.c)
    gphase = cmath.exp(1j * weyl.global_phase)
    coeffs = {name: gphase * c for name, c in coeffs.items()}
    # qiskit is little-endian vs qibo big-endian: try both control/target orderings.
    for kla, klb, kra, krb in (
        (weyl.K1l, weyl.K1r, weyl.K2l, weyl.K2r),
        (weyl.K1r, weyl.K1l, weyl.K2r, weyl.K2l),
    ):
        terms = _terms_from_factors(coeffs, kla, klb, kra, krb)
        recon = sum(c * np.kron(v, w) for c, v, w in terms)
        if np.max(np.abs(recon - u)) < 1e-9:
            return terms
    raise RuntimeError("KAK decomposition failed to reconstruct the gate")


def _pauli_local_terms(u):
    """Backend-free fallback ``U -> [(c_i, V_i, W_i)]``: expand ``U`` in the two-qubit
    Pauli basis ``U = sum_{P,Q} c_{PQ} (P (x) Q)`` (``c_{PQ} = Tr[(P (x) Q)^dag U]/4``).
    Always works and every factor is a (phased) Pauli -- hence unitary -- but with up
    to 16 terms it is *not* gamma-optimal. Used only when qiskit (the KAK engine) is
    unavailable."""
    u = np.asarray(u, dtype=complex)
    terms = []
    for pn, p in _PAULI.items():
        for qn, q in _PAULI.items():
            kron = np.kron(p, q)
            coeff = np.trace(kron.conj().T @ u) / 4
            if abs(coeff) < 1e-12:
                continue
            terms.append((abs(coeff), (coeff / abs(coeff)) * p, q))
    return terms


def _general_local_terms(u):
    """``U -> [(c_i, V_i, W_i)]`` for an arbitrary two-qubit unitary: the gamma-optimal
    KAK decomposition when qiskit is available, else the (correct, suboptimal)
    Pauli-basis fallback."""
    try:
        return _kak_local_terms(u)
    except Exception as exc:  # noqa: BLE001 -- any failure => use the robust fallback
        import logging

        logging.getLogger(__name__).debug(
            "KAK unavailable (%s); using Pauli basis", exc
        )
        return _pauli_local_terms(u)


def _two_qubit_matrix(gate_name, params=()):
    """The 4x4 unitary of a named two-qubit qibo gate (for the general cut path)."""
    import qibo

    from Qdislib.utils.graph_qibo import _resolve_qibo_gate

    try:
        gate = _resolve_qibo_gate(str(gate_name))(0, 1, *params)
        matrix = np.asarray(gate.matrix(qibo.get_backend()))
    except Exception as exc:  # noqa: BLE001
        raise NotImplementedError(
            f"cannot build a two-qubit matrix for {gate_name!r}: {exc}"
        )
    if matrix.shape != (4, 4):
        raise NotImplementedError(
            f"Hadamard-test cut needs a two-qubit gate; {gate_name!r} is {matrix.shape}"
        )
    return matrix


def local_unitary_decomposition(gate_name, params=()):
    """``U -> [(c_i, V_i, W_i)]``: the local-unitary decomposition of a cut gate.

    ``V_i`` acts on the control-side fragment qubit, ``W_i`` on the target-side, every
    factor is unitary, ``c_i > 0``, and ``sum_i c_i V_i (x) W_i`` equals the gate up to
    a global phase. The gates in :data:`SUPPORTED_HADAMARD_GATES` (whose only non-local
    part is a single ``ZZ`` interaction) take a cheap closed-form path giving two terms
    with the optimal 1-norm; ``CX``/``CNOT`` reduce to ``CZ`` (the ``H`` on the target
    is local, applied inside the fragment by the caller). **Any other two-qubit gate**
    falls through to the general KAK / Cartan decomposition
    (:func:`_general_local_terms`), so arbitrary non-diagonal gates (``iSWAP``,
    ``SWAP``, ``fSim``, ...) are cuttable too.

    :param gate_name: Canonical gate name (case-insensitive).
    :param params: Gate parameters (e.g. ``(theta,)`` for ``CP``/``RZZ``).
    :raises NotImplementedError: if ``gate_name`` is not a resolvable two-qubit gate.
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
    else:  # general non-diagonal two-qubit gate: KAK / Cartan (optimal) or Pauli basis
        return _general_local_terms(_two_qubit_matrix(gate_name, params))
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
    ``g == 0`` else ``S`` (real- vs imaginary-part test). ``U_i, U_j`` are *general*
    single-qubit unitaries (the KAK path produces non-diagonal ones), so the select
    gate is realised as ``U_i`` unconditionally followed by the controlled
    ``U_i^dag U_j``: ``U3(U_i)`` then ``CU3(U_i^dag U_j)`` plus a ``U1`` on the ancilla
    for that operator's global phase. ``U_i``'s own global phase is a fragment-global
    phase (irrelevant to the Z read-out) and is dropped. Emitting only standard named
    gates (H/S/U3/CU3/U1/CNOT) keeps the fragment expressible in the neutral DAG and
    runnable on any backend. The ancilla is read out as a Z expectation by the caller.
    """
    theta_i, phi_i, lam_i, _ = _u3_angles(u_i)  # drop U_i's fragment-global phase
    # U_i is applied unconditionally first, so the controlled correction that turns the
    # anc==1 branch into U_j must act *after* U_i: M @ U_i = U_j => M = U_j U_i^dag.
    # (For diagonal U_i, U_j this equals U_i^dag U_j since they commute -- but not in
    # general, so the order matters for non-diagonal KAK terms.)
    m = np.asarray(u_j) @ np.asarray(u_i).conj().T
    theta_m, phi_m, lam_m, gphase_m = _u3_angles(m)
    circuit.add(gates.H(ancilla))
    if g == 1:
        circuit.add(gates.S(ancilla))
    circuit.add(gates.U3(qubit, theta_i, phi_i, lam_i))  # U_i on the data qubit
    circuit.add(
        gates.CU3(ancilla, qubit, theta_m, phi_m, lam_m)
    )  # controlled U_j U_i^dag
    if abs(gphase_m) > 1e-12:  # that operator's global phase, on the anc==1 branch
        circuit.add(gates.U1(ancilla, gphase_m))
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


def _z_observable(num_qubits, z_qubits):
    """The per-qubit Pauli string with ``Z`` on ``z_qubits`` and ``I`` elsewhere --
    the fragment observable (joint ``Z`` over the ancillas + the observed qubits)."""
    obs = ["I"] * num_qubits
    for q in z_qubits:
        obs[q] = "Z"
    return "".join(obs)


def _hadamard_setup(circuit, gates_cut):
    """Shared preamble for every Hadamard-cut driver: the per-cut local-unitary
    decompositions, the cut DAG (SUBS1/SUBS2 placeholders), and its fixed fragment
    partition. Returns ``(cut_dag, components, decomps)``."""
    import networkx as nx

    from Qdislib.utils.backends import convert_to_dag

    dag = convert_to_dag(circuit)
    decomps = []
    for name in gates_cut:
        gate = (dag.nodes[name].get("gate") or "").upper()
        params = dag.nodes[name].get("parameters") or ()
        decomps.append(local_unitary_decomposition(gate, params))

    cut_dag = _generate_cut_hadamard(dag, gates_cut)  # SUBS1_k / SUBS2_k (+ CX H-conj)
    # The fragment partition (connected components) is fixed by the cut graph; only the
    # per-node gadget assignment changes per configuration, so compute it once.
    components = list(
        nx.connected_components(_remove_red_free(cut_dag).to_undirected())
    )
    return cut_dag, components, decomps


def _config_assignments(decomps, config):
    """SUBS-node -> ``(V_i, V_j, g)`` / ``(W_i, W_j, g)`` gadget assignments for one
    ``config`` (a list of ``(i, j, g)`` per cut)."""
    assignments = {}
    for k, (i, j, g) in enumerate(config, start=1):
        _, v_i, w_i = decomps[k - 1][i]
        _, v_j, w_j = decomps[k - 1][j]
        assignments[f"SUBS1_{k}"] = (v_i, v_j, g)  # control side -> V
        assignments[f"SUBS2_{k}"] = (w_i, w_j, g)  # target side -> W
    return assignments


def _config_fragments(cut_dag, components, assignments, observable):
    """Build the ``[(fragment_circuit, z_qubits), ...]`` for one configuration's gadget
    assignments; ``z_qubits`` = ancillas + observed (``Z``) qubits, renumbered."""
    fragments = []
    for comp_nodes in components:
        comp = cut_dag.subgraph(comp_nodes).copy()
        frag, qmap, ancillas = _component_hadamard_circuit(comp, assignments)
        z_qubits = list(ancillas)
        for q, ch in enumerate(observable):
            if ch == "Z" and q in qmap:
                z_qubits.append(qmap[q])
        fragments.append((frag, z_qubits))
    return fragments


def _cut_pmf(terms):
    """The sampling PMF for one cut: ``(settings, probs, phi)`` where ``settings`` are
    the valid ``(i, j, g)``, ``probs[t] = c_i c_j / phi`` (the Harrow-Lowe distribution
    ``p(i,j,g) = c_i c_j / phi``), and ``phi = sum_{valid} c_i c_j = 2||c||_1^2 -
    ||c||_2^2`` is the per-cut 1-norm. The sign ``(-1)^g`` is *not* in the PMF -- it
    rides in the estimator -- so the weights are non-negative and sum to ``phi``."""
    settings = _cut_settings(len(terms))
    weights = [terms[i][0] * terms[j][0] for (i, j, g) in settings]
    phi = sum(weights)
    probs = [w / phi for w in weights]
    return settings, probs, phi


def _hadamard_fragment_specs(circuit, gates_cut, observable):
    """Enumerate the double-Hadamard-test fragments for a set of cuts.

    Returns a list with one record per non-zero enumeration setting:
    ``(weight, [(fragment_circuit, z_qubits), ...])``. Each fragment is an independent
    qibo circuit carrying its local ancilla(s); ``z_qubits`` are the qubits read out as
    a joint ``Z`` (ancillas + observed qubits). The reconstruction of ``<observable>``
    is ``sum_record weight * prod_fragment <Z_{z_qubits}>`` -- the shared core of the
    all-in-one and three-step *enumeration* drivers.
    """
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut)
    settings_per_cut = [_cut_settings(len(terms)) for terms in decomps]

    specs = []
    for config in itertools.product(*settings_per_cut):  # one (i, j, g) per cut
        weight = 1.0
        for k, (i, j, g) in enumerate(config, start=1):
            weight *= decomps[k - 1][i][0] * decomps[k - 1][j][0] * (-1) ** g
        if weight == 0.0:
            continue
        assignments = _config_assignments(decomps, config)
        fragments = _config_fragments(cut_dag, components, assignments, observable)
        specs.append((weight, fragments))
    return specs


@task(returns=1)
def _evaluate_hadamard_fragment(frag_dag, num_qubits, software, shots, qpu, qpu_dict):
    """Evaluate one fragment's ``<Z_{z_qubits}>`` on any backend (a COMPSs task).

    The fragment rides as a neutral DAG whose non-read-out qubits are marked
    "Observable I", so :meth:`build_execution_backend(...).evaluate_expectation`
    returns exactly the joint ``Z`` over the ancillas + observed qubits. This is the
    per-fragment unit of work fanned out across PyCOMPSs and routed (qibo / Qiskit /
    CUDA-Q / PennyLane / a QPU) by ``software``/``qpu``.
    """
    from Qdislib.core.backend.registry import build_execution_backend

    backend = build_execution_backend(
        software, num_qubits=num_qubits, qpu=bool(qpu), qpu_dict=qpu_dict
    )
    return backend.evaluate_expectation(frag_dag, num_qubits, shots=shots)


@task(returns=1, values=COLLECTION_IN)
def _reduce_hadamard(values, weights, config_of):
    """Fold the per-fragment values into ``<observable>`` (the fan-in COMPSs task).

    Takes the whole fan-out (``values`` a ``COLLECTION_IN`` of the fragment futures)
    and returns ``sum_config weight * prod_fragment value`` -- so the task graph is a
    clean fan-out of fragment tasks into this single reduction, and the reduction runs
    off the master.
    """
    return gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of)


def gate_cutting_hadamard(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    shots=None,
    seed=None,
    qpu=False,
    qpu_dict=None,
):
    """Estimate ``<observable>`` via the ancilla / double-Hadamard-test gate cut.

    Cuts the two-qubit ``gates_cut`` (names like ``"CZ_2"``, ``"ISWAP_3"``) with the
    Harrow-Lowe space-like cut (arXiv:2403.01018): each fragment carries one local
    ancilla running a Hadamard test, the fragments are independent, and the result is
    combined classically as ``sum_config weight * prod_fragment <Z>``.

    With ``software="qibo"`` and ``shots=None`` this is the **exact** (statevector,
    full-enumeration) reconstruction -- the correctness anchor. Otherwise each fragment
    is evaluated by :func:`build_execution_backend` on the chosen backend (qibo /
    Qiskit / CUDA-Q / PennyLane / a QPU) from ``shots`` samples, fanned out as
    independent PyCOMPSs tasks (one per fragment per configuration).

    :param circuit: A qibo circuit.
    :param gates_cut: List of cut gate names (each a two-qubit gate; the ZZ-diagonal
        gates use the optimal fast path, others the general KAK decomposition).
    :param observable: A Z/I Pauli string over the circuit's qubits.
    :param software: Execution backend for the sampled path (``"qibo"``, ``"qiskit"``,
        ``"cudaq"``, ``"pennylane"``).
    :param shots: Samples per fragment; ``None`` selects the exact statevector path
        (``software="qibo"`` only).
    :param seed: Unused placeholder kept for signature compatibility.
    :param qpu: Route fragment execution to a QPU via ``build_execution_backend``.
    :param qpu_dict: QPU connection settings forwarded to the backend.
    :return: The reconstructed expectation value (resolved, not a future).
    """
    exact = software == "qibo" and shots is None
    specs = _hadamard_fragment_specs(circuit, gates_cut, observable)

    # Flat fan-out: one value (exact) or future (@task) per fragment, with parallel
    # weights / config_of bookkeeping so a single reduction folds them back.
    values, config_of, weights = [], [], []
    for cfg_idx, (weight, fragments) in enumerate(specs):
        weights.append(weight)
        for frag, z_qubits in fragments:
            if exact:
                values.append(_expectation_z(frag, z_qubits))
            else:
                obs = _z_observable(frag.nqubits, z_qubits)
                frag_dag = circuit_qibo_to_dag(frag, obs_I=obs)
                values.append(
                    _evaluate_hadamard_fragment(
                        frag_dag, frag.nqubits, software, shots or 1024, qpu, qpu_dict
                    )
                )
            config_of.append(cfg_idx)

    if exact:  # in-process, no tasks to synchronise
        return gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
    # Fan the per-fragment futures into one reduction task (fan-out -> fan-in).
    return compss_wait_on(_reduce_hadamard(values, weights, config_of))


class HadamardCutSubcircuits(list):
    """The Hadamard-cut fragments plus the metadata to recombine them.

    Behaves like the ``list`` of subcircuits it is (iterate/index/``len``), and adds
    :meth:`reconstruct` so the three-step workflow needs no manual bookkeeping::

        subs = qd.gate_cutting_hadamard_subcircuits(circuit, ["CZ_2"], "ZZ")
        values = [qd.expectation_value(s) for s in subs]
        result = subs.reconstruct(values)

    Each item is a ``[fragment_circuit, observable]`` pair; every fragment belongs to
    one enumeration configuration (``config_of[i]``) whose signed coefficient is
    ``weights[config_of[i]]``. Reconstruction multiplies the values within a
    configuration and sums the weighted products.
    """

    def __init__(self, subcircuits, *, weights, config_of, number_cuts):
        super().__init__(subcircuits)
        #: Signed coefficient of each configuration (``c_i c_j (-1)^g`` product).
        self.weights = weights
        #: Configuration index each subcircuit (in order) belongs to.
        self.config_of = config_of
        #: Number of gates cut.
        self.number_cuts = number_cuts

    def reconstruct(self, values):
        """Recombine per-subcircuit expectation values into ``<observable>``.

        :param values: One value per subcircuit, in this object's order.
        :return: The reconstructed expectation value of the uncut circuit.
        """
        return gate_cutting_hadamard_subcircuit_reconstruction(
            values, self.weights, self.config_of
        )


def _materialize_subcircuit(frag, z_qubits, software):
    """Turn a fragment qibo circuit into a ``[circuit, observable]`` subcircuit for the
    requested backend, so :func:`expectation_value` can evaluate it directly."""
    obs = _z_observable(frag.nqubits, z_qubits)
    if software == "qibo":
        return [frag, obs]
    # Convert through the neutral DAG (the non-read-out qubits become "Observable I",
    # which each backend converter turns back into the observable string).
    frag_dag = circuit_qibo_to_dag(frag, obs_I=obs)
    n = frag.nqubits
    if software == "qiskit":
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        return dag_to_circuit_qiskit(frag_dag, n)
    if software in ("cudaq", "cuda-q", "cuda_q"):
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq

        return dag_to_circuit_cudaq(frag_dag, n)
    if software == "pennylane":
        from Qdislib.utils.graph_pennylane import dag_to_circuit_pennylane

        return dag_to_circuit_pennylane(frag_dag, n)
    raise NotImplementedError(f"unknown software backend {software!r}")


def gate_cutting_hadamard_subcircuits(circuit, gates_cut, observable, software="qibo"):
    """Build the Hadamard-cut fragments as standalone subcircuits (step one of three).

    Returns a :class:`HadamardCutSubcircuits` -- a list of ``[fragment, observable]``
    subcircuits that the caller evaluates independently (e.g. with
    :func:`expectation_value`, on any backend / QPU, in parallel) and then feeds to
    :meth:`~HadamardCutSubcircuits.reconstruct`. This is the explicit three-step form
    of :func:`gate_cutting_hadamard`.

    :param circuit: A qibo circuit.
    :param gates_cut: List of cut gate names.
    :param observable: A Z/I Pauli string over the circuit's qubits.
    :param software: Backend the subcircuits are materialised for.
    :return: A :class:`HadamardCutSubcircuits`.
    """
    specs = _hadamard_fragment_specs(circuit, gates_cut, observable)
    subcircuits, weights, config_of = [], [], []
    for cfg_idx, (weight, fragments) in enumerate(specs):
        weights.append(weight)
        for frag, z_qubits in fragments:
            subcircuits.append(_materialize_subcircuit(frag, z_qubits, software))
            config_of.append(cfg_idx)
    return HadamardCutSubcircuits(
        subcircuits,
        weights=weights,
        config_of=config_of,
        number_cuts=len(gates_cut),
    )


def gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of):
    """Recombine per-fragment values into ``<observable>`` (step three of three).

    :param values: One expectation value per subcircuit, in the order produced by
        :func:`gate_cutting_hadamard_subcircuits`.
    :param weights: Signed coefficient of each configuration.
    :param config_of: Configuration index of each subcircuit (parallel to ``values``).
    :return: ``sum_config weight * prod_fragment value``.
    """
    products = [1.0] * len(weights)
    for value, cfg in zip(values, config_of):
        if value is None:  # a skipped/empty fragment contributes the identity factor
            continue
        products[cfg] *= value
    return sum(w * p for w, p in zip(weights, products))


# --------------------------------------------------------------------------- #
# Monte-Carlo sampling driver (the native Harrow-Lowe estimator)
# --------------------------------------------------------------------------- #


def _sample_configs(decomps, num_samples, seed=None):
    """Draw ``num_samples`` cut configurations from the Harrow-Lowe PMF.

    Returns ``(configs, phi_total)`` where each config is ``([(i, j, g), ...], sign)``
    with ``sign = prod_k (-1)^{g_k}`` and ``phi_total = prod_k phi_k`` the total 1-norm.
    Each cut ``k`` draws ``(i, j, g)`` from ``p = c_i c_j / phi_k`` independently.
    """
    pmfs = [_cut_pmf(terms) for terms in decomps]  # (settings, probs, phi) per cut
    phi_total = 1.0
    for _, _, phi in pmfs:
        phi_total *= phi
    rng = np.random.default_rng(seed)
    configs = []
    for _ in range(num_samples):
        config, sign = [], 1
        for settings, probs, _phi in pmfs:
            i, j, g = settings[int(rng.choice(len(settings), p=probs))]
            config.append((i, j, g))
            sign *= (-1) ** g
        configs.append((config, sign))
    return configs, phi_total


def _sampling_weight(phi_total, sign, num_samples):
    """The per-sample coefficient that turns the shared weighted-product reconstruction
    into the Monte-Carlo estimator: with ``weight_s = phi_total * sign_s / N`` the sum
    ``sum_s weight_s * prod_frag value`` equals ``phi_total * mean_s(sign_s * prod
    value)`` -- the unbiased Harrow-Lowe estimator of ``<observable>``."""
    return phi_total * sign / num_samples


def gate_cutting_hadamard_sampling(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    num_samples=1000,
    shots=1024,
    seed=None,
    qpu=False,
    qpu_dict=None,
):
    """Estimate ``<observable>`` by **Monte-Carlo sampling** the Hadamard-test cut.

    The native Harrow-Lowe estimator (arXiv:2403.01018): instead of enumerating all
    ``(2m^2-m)^{n_cuts}`` configurations (as :func:`gate_cutting_hadamard`), draw
    ``num_samples`` configurations from the quasiprobability PMF ``p(i,j,g) = c_i c_j /
    phi`` and average ``phi_total * sign * prod_fragment <Z>`` (``sign = prod (-1)^g``).
    This bounds the cost by ``num_samples`` regardless of the number of cuts -- the
    enumeration is exponential in ``n_cuts``, so sampling is what makes many-cut
    circuits tractable (at the cost of a ``phi_total^2 / num_samples`` variance).

    Fragments fan out as independent PyCOMPSs tasks and fold through the same
    :func:`_reduce_hadamard` fan-in as the enumeration driver.

    :param circuit: A qibo circuit.
    :param gates_cut: List of cut gate names.
    :param observable: A Z/I Pauli string over the circuit's qubits.
    :param software: Fragment-execution backend (``"qibo"``/``"qiskit"``/``"cudaq"``/
        ``"pennylane"``).
    :param num_samples: Number of configurations drawn from the PMF.
    :param shots: Samples per fragment; ``None`` uses the exact statevector per
        fragment (``software="qibo"`` only) -- then only the *configuration* sampling
        is stochastic.
    :param seed: Seed for the configuration draws (reproducible).
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :return: The estimated expectation value.
    """
    exact = software == "qibo" and shots is None
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut)
    configs, phi_total = _sample_configs(decomps, num_samples, seed)

    values, config_of, weights = [], [], []
    for s, (config, sign) in enumerate(configs):
        weights.append(_sampling_weight(phi_total, sign, num_samples))
        assignments = _config_assignments(decomps, config)
        for frag, z_qubits in _config_fragments(
            cut_dag, components, assignments, observable
        ):
            if exact:
                values.append(_expectation_z(frag, z_qubits))
            else:
                obs = _z_observable(frag.nqubits, z_qubits)
                frag_dag = circuit_qibo_to_dag(frag, obs_I=obs)
                values.append(
                    _evaluate_hadamard_fragment(
                        frag_dag, frag.nqubits, software, shots or 1024, qpu, qpu_dict
                    )
                )
            config_of.append(s)

    if exact:  # in-process, no tasks to synchronise
        return gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
    return compss_wait_on(_reduce_hadamard(values, weights, config_of))


def gate_cutting_hadamard_sampling_subcircuits(
    circuit, gates_cut, observable, software="qibo", num_samples=1000, seed=None
):
    """Monte-Carlo-sampled Hadamard-cut fragments as standalone subcircuits.

    The sampling counterpart of :func:`gate_cutting_hadamard_subcircuits`: draws
    ``num_samples`` configurations from the PMF and returns their fragments as a
    :class:`HadamardCutSubcircuits` (already carrying the ``phi_total * sign /
    num_samples`` weights), so the usual ``evaluate -> reconstruct`` workflow yields the
    Monte-Carlo estimate::

        subs = qd.gate_cutting_hadamard_sampling_subcircuits(c, ["CZ_2"], "ZZ", num_samples=2000)
        values = [qd.expectation_value(s) for s in subs]
        estimate = subs.reconstruct(values)

    :param num_samples: Number of configurations drawn from the PMF.
    :param seed: Seed for the configuration draws.
    :return: A :class:`HadamardCutSubcircuits`.
    """
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut)
    configs, phi_total = _sample_configs(decomps, num_samples, seed)

    subcircuits, weights, config_of = [], [], []
    for s, (config, sign) in enumerate(configs):
        weights.append(_sampling_weight(phi_total, sign, num_samples))
        assignments = _config_assignments(decomps, config)
        for frag, z_qubits in _config_fragments(
            cut_dag, components, assignments, observable
        ):
            subcircuits.append(_materialize_subcircuit(frag, z_qubits, software))
            config_of.append(s)
    return HadamardCutSubcircuits(
        subcircuits, weights=weights, config_of=config_of, number_cuts=len(gates_cut)
    )


def _remove_red_free(dag):
    """Drop the red (cut-dependency) edges, leaving the blue circuit-flow edges that
    define the independent fragments (mirrors ``_remove_red_edges``)."""
    from Qdislib.utils.graph_qibo import _remove_red_edges

    return _remove_red_edges(dag.copy())
