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
    "HadamardSampledCutSubcircuits",
]

import itertools

from Qdislib.core.cutting_algorithms._reconstruction import next_seed
from Qdislib.utils.backends import Subcircuit
from Qdislib.core.cutting_algorithms._api import (
    standard_arguments,
    resolve_backend,
    resolve_software,
    default_observable,
)
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


class _ProductTerms:
    """The decomposition of a joint group: the products of its gates' own terms,
    ``(c_s, [(V, W) per gate])`` with ``c_s = prod_m c_{s_m}``, built on demand.

    A group of k gates has ``prod_m r_m`` terms (``2**k`` for CZ/CX/RZZ), so they are
    never listed: ``terms[s]`` builds term ``s`` from its per-gate indices (first gate
    most significant), and the norms and counts come from the gates' own terms.
    """

    def __init__(self, per_gate):
        self.per_gate = [list(t) for t in per_gate]
        self.radix = [len(t) for t in self.per_gate]
        self._len = int(np.prod(self.radix)) if self.radix else 1
        #: ``||c||_1`` and ``||c||_2**2`` of the product, and its non-zero terms.
        self.l1 = float(np.prod([sum(abs(c) for c, _, _ in t) for t in self.per_gate]))
        self.l2sq = float(
            np.prod([sum(abs(c) ** 2 for c, _, _ in t) for t in self.per_gate])
        )
        self.nonzero = int(
            np.prod([sum(1 for c, _, _ in t if c != 0) for t in self.per_gate])
        )

    def __len__(self):
        return self._len

    def indices(self, s):
        """Term ``s`` as one index per gate."""
        out = []
        for r in reversed(self.radix):
            out.append(s % r)
            s //= r
        return out[::-1]

    def __getitem__(self, s):
        if not 0 <= s < self._len:
            raise IndexError(s)
        c, factors = 1.0, []
        for terms, i in zip(self.per_gate, self.indices(s)):
            c_m, v, w = terms[i]
            c *= c_m
            factors.append((v, w))
        return c, factors

    def __iter__(self):
        return (self[s] for s in range(self._len))


def phi_factor(terms):
    """The sampling-overhead 1-norm ``phi = 2||c||_1^2 - ||c||_2^2`` of a decomposition."""
    if isinstance(terms, _ProductTerms):
        return float(2 * terms.l1**2 - terms.l2sq)
    c = np.array([t[0] for t in terms], dtype=float)
    return float(2 * c.sum() ** 2 - (c**2).sum())


def _num_settings(terms):
    """How many ``(s, t, g)`` settings carry a non-zero weight: ``2 r**2 - r`` over the
    ``r`` terms with a non-zero coefficient (``s = t`` has only ``g = 0``)."""
    if isinstance(terms, _ProductTerms):
        r = terms.nonzero
    else:
        r = sum(1 for t in terms if t[0] != 0)
    return 2 * r * r - r


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


def _open_ancilla(circuit, gates, ancilla, g):
    """Start a Hadamard test on ``ancilla``: ``H``, then ``S`` for the imaginary part."""
    circuit.add(gates.H(ancilla))
    if g == 1:
        circuit.add(gates.S(ancilla))


def _select(circuit, gates, qubit, ancilla, u_i, u_j):
    """The select gate ``|0><0| (x) U_i + |1><1| (x) U_j`` on ``qubit``, ``ancilla``
    the control: ``U3(U_i)`` then the controlled ``U_j U_i^dag`` (see
    :func:`_subs_gadget`)."""
    theta_i, phi_i, lam_i, _ = _u3_angles(u_i)  # drop U_i's fragment-global phase
    # U_i is applied unconditionally first, so the controlled correction that turns the
    # anc==1 branch into U_j must act *after* U_i: M @ U_i = U_j => M = U_j U_i^dag.
    # (For diagonal U_i, U_j this equals U_i^dag U_j since they commute -- but not in
    # general, so the order matters for non-diagonal KAK terms.)
    m = np.asarray(u_j) @ np.asarray(u_i).conj().T
    theta_m, phi_m, lam_m, gphase_m = _u3_angles(m)
    circuit.add(gates.U3(qubit, theta_i, phi_i, lam_i))  # U_i on the data qubit
    circuit.add(
        gates.CU3(ancilla, qubit, theta_m, phi_m, lam_m)
    )  # controlled U_j U_i^dag
    if abs(gphase_m) > 1e-12:  # that operator's global phase, on the anc==1 branch
        circuit.add(gates.U1(ancilla, gphase_m))


def _close_ancilla(circuit, gates, ancilla):
    """End a Hadamard test on ``ancilla`` (read out later as a ``Z``)."""
    circuit.add(gates.H(ancilla))


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

    It is :func:`_open_ancilla`, :func:`_select` and :func:`_close_ancilla` in a row;
    the joint cut spreads those over several cut gates sharing one ancilla.
    """
    _open_ancilla(circuit, gates, ancilla, g)
    _select(circuit, gates, qubit, ancilla, u_i, u_j)
    _close_ancilla(circuit, gates, ancilla)


def _component_hadamard_circuit(component, assignments):
    """Build one fragment's Hadamard-test circuit from its DAG component.

    Walks the component in topological order, renumbering its qubits to ``0..k-1``.
    ``assignments`` maps a SUBS node id to ``(u_i, u_j, g, ancilla_key)`` (the V- or
    W-side unitaries). SUBS nodes with the same ``ancilla_key`` share one ancilla
    (indices ``k, k+1, ...``): it is opened just before the first of them, gets a
    select at each, and is closed just after the last -- for a single node, exactly
    :func:`_subs_gadget`. A 3-tuple (no key) gets its own ancilla, as before.
    Returns ``(circuit, qubit_map, ancilla_indices)`` where ``qubit_map`` is
    original-qubit -> new-index.
    """
    import networkx as nx
    from qibo import gates, models

    from Qdislib.utils.graph_qibo import _resolve_qibo_gate

    orig_qubits = sorted(
        {q for _, d in component.nodes(data=True) for q in (d.get("qubits") or ())}
    )
    qubit_map = {q: i for i, q in enumerate(orig_qubits)}
    n = len(orig_qubits)
    order = list(nx.topological_sort(component))

    def ancilla_key(node):
        assignment = assignments[node]
        return assignment[3] if len(assignment) > 3 else ("node", node)

    subs = [node for node in order if str(node).startswith("SUBS")]
    keys = list(dict.fromkeys(ancilla_key(node) for node in subs))  # first-use order
    anc_of = {key: n + index for index, key in enumerate(keys)}
    last_of = {ancilla_key(node): node for node in subs}  # last SUBS node per ancilla
    circuit = models.Circuit(n + len(keys))
    opened = set()
    for node in order:
        data = component.nodes[node]
        name = (data.get("gate") or "").upper()
        qubits = tuple(qubit_map[q] for q in (data.get("qubits") or ()))
        if str(node).startswith("SUBS"):
            u_i, u_j, g = assignments[node][:3]
            key = ancilla_key(node)
            ancilla = anc_of[key]
            if key not in opened:
                _open_ancilla(circuit, gates, ancilla, g)
                opened.add(key)
            _select(circuit, gates, qubits[0], ancilla, u_i, u_j)
            if last_of[key] == node:
                _close_ancilla(circuit, gates, ancilla)
        elif name in ("OBSERVABLE I", "MEASURE"):
            continue
        else:
            gate_cls = _resolve_qibo_gate(data.get("gate") or "")
            params = data.get("parameters")
            if params:
                circuit.add(gate_cls(*qubits, *params))
            else:
                circuit.add(gate_cls(*qubits))
    return circuit, qubit_map, [anc_of[key] for key in keys]


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


def _validate_hadamard_cuts_partition(gates_cut, cut_graph):
    """Raise unless every Hadamard gate cut separates the circuit (see
    :func:`~Qdislib.core.cutting_algorithms._reconstruction.check_cuts_separate`).

    :param gates_cut: Gate node ids to cut.
    :param cut_graph: The cut DAG without its red edges.
    :raises InvalidCutError: If a cut gate's two sides stay in one fragment.
    """
    from Qdislib.core.cutting_algorithms._reconstruction import check_cuts_separate

    check_cuts_separate(
        cut_graph,
        {
            gate: (f"SUBS1_{k}", f"SUBS2_{k}")
            for k, gate in enumerate(gates_cut, start=1)
        },
    )


def _observable_terms(observable):
    """``[(coefficient, pauli)]`` for a Pauli string or a weighted sum of them
    (anything :func:`~Qdislib.utils.backends.normalize_observables` accepts)."""
    from Qdislib.utils.backends import normalize_observables

    terms = (
        [(1.0, observable)]
        if isinstance(observable, str)
        else list(normalize_observables(observable))
    )
    for _, pauli in terms:
        if not isinstance(pauli, str) or set(pauli) - set("IXYZ"):
            raise ValueError(
                f"Observable terms must be strings over I, X, Y, Z, got {pauli!r}."
            )
    return terms


def _is_single_term(terms):
    return len(terms) == 1 and terms[0][0] == 1.0


def _z_basis(observable):
    """The observable after :func:`_hadamard_setup` rotated X/Y into Z."""
    return observable.replace("X", "Z").replace("Y", "Z")


@task(returns=1, values=COLLECTION_IN)
def _weighted_sum(values, coefficients):
    """``sum_k c_k * value_k`` -- folds a Pauli sum's per-term results."""
    return sum(c * v for c, v in zip(coefficients, values))


def _hadamard_setup(circuit, gates_cut, observable=None):
    """Shared preamble for every Hadamard-cut driver: the per-cut local-unitary
    decompositions, the cut DAG (SUBS1/SUBS2 placeholders), and its fixed fragment
    partition. Returns ``(cut_dag, components, decomps)``.

    An X or Y in ``observable`` is measured by rotating that qubit into the Z basis
    at its end (``H``, or ``SDG`` then ``H``), as ``gate_cutting`` does; read the
    fragments with :func:`_z_basis` of it."""
    import networkx as nx

    from Qdislib.utils.backends import change_basis_dag, convert_to_dag

    dag = convert_to_dag(circuit)
    if observable and set(observable) & {"X", "Y"}:
        dag = change_basis_dag(dag, observable)
    decomps = []
    for name in gates_cut:
        gate = (dag.nodes[name].get("gate") or "").upper()
        params = dag.nodes[name].get("parameters") or ()
        decomps.append(local_unitary_decomposition(gate, params))

    cut_dag = _generate_cut_hadamard(dag, gates_cut)  # SUBS1_k / SUBS2_k (+ CX H-conj)
    # The fragment partition (connected components) is fixed by the cut graph; only the
    # per-node gadget assignment changes per configuration, so compute it once.
    cut_graph = _remove_red_free(cut_dag)
    components = list(nx.connected_components(cut_graph.to_undirected()))
    _validate_hadamard_cuts_partition(gates_cut, cut_graph)
    return cut_dag, components, decomps


def _cut_groups(cut_dag, components, decomps, joint=False):
    """Group the cuts that are decomposed (and Hadamard-tested) together.

    Returns ``[(cut_indices, terms), ...]``: ``cut_indices`` are 1-based positions in
    ``gates_cut``; ``terms`` the group's decomposition ``[(c_s, [(V, W) per gate in
    cut_indices order]), ...]``.

    ``joint=False``: one group per cut, today's independent cut (``6**k`` settings,
    ``gamma = 3**k`` for k CZ/CX). ``joint=True``: the cuts crossing the same pair of
    fragments form one group, whose terms are the products of the gates' own terms,
    ``c_s = prod_m c_{s_m}`` (Schmitt, Piveteau & Sutter, Quantum 9, 1634 (2025)). The
    product is a decomposition of the whole boundary even when its gates sit at
    different times, since everything between them is local to one side. Sampled
    with one ancilla per fragment it costs ``2**k (2**(k+1) - 1)`` settings and
    ``gamma = 2**(k+1) - 1`` for k CZ/CX (see ``docs/design/joint_hadamard_cut.md``).
    A cut whose two sides land in the same fragment stays a group of its own.
    """
    component_of = {}
    for index, nodes in enumerate(components):
        for node in nodes:
            component_of[node] = index
    keys = {}
    for k in range(1, len(decomps) + 1):
        a, b = component_of[f"SUBS1_{k}"], component_of[f"SUBS2_{k}"]
        key = frozenset((a, b)) if joint and a != b else ("cut", k)
        keys.setdefault(key, []).append(k)
    groups = []
    for cut_indices in keys.values():  # in order of each group's first cut
        # the product of the gates' terms, built term by term when needed
        terms = _ProductTerms([decomps[k - 1] for k in cut_indices])
        groups.append((cut_indices, terms))
    return groups


def _config_assignments(groups, config):
    """SUBS-node -> ``(V_s, V_t, g, ancilla_key)`` / ``(W_s, W_t, g, ancilla_key)``
    gadget assignments for one ``config`` (an ``(s, t, g)`` per group). Every gate of
    a group gets its own factor of the group's terms ``s`` and ``t`` and the group's
    ``g``. A group of one cut keeps an ancilla per side (keyed by node); a joint
    group shares one per fragment (keyed by group)."""
    assignments = {}
    for index, ((cut_indices, terms), (s, t, g)) in enumerate(zip(groups, config)):
        _, factors_s = terms[s]
        _, factors_t = terms[t]
        joint = len(cut_indices) > 1
        for position, k in enumerate(cut_indices):
            v_s, w_s = factors_s[position]
            v_t, w_t = factors_t[position]
            for node, u_s, u_t in ((f"SUBS1_{k}", v_s, v_t), (f"SUBS2_{k}", w_s, w_t)):
                key = ("group", index) if joint else ("node", node)
                assignments[node] = (u_s, u_t, g, key)  # SUBS1 -> V, SUBS2 -> W
    return assignments


def _hadamard_cut_overhead(circuit, gates_cut, joint=False):
    """Sampling overhead of a Hadamard cut: ``(phi_total, [(gates, phi), ...])``.

    ``phi`` is the 1-norm of each group's decomposition (:func:`phi_factor`), the
    shots grow as ``phi_total**2``. For k CZ/CX across one boundary: ``3**k`` cut
    separately, ``2**(k+1) - 1`` with ``joint=True``.
    """
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut)
    per_group = [
        ([gates_cut[k - 1] for k in cut_indices], phi_factor(terms))
        for cut_indices, terms in _cut_groups(cut_dag, components, decomps, joint)
    ]
    total = 1.0
    for _, phi in per_group:
        total *= phi
    return total, per_group


def _config_fragments(cut_dag, components, assignments, observable):
    """Build the ``[(fragment_circuit, z_qubits), ...]`` for one configuration's gadget
    assignments; ``z_qubits`` = ancillas + observed (``Z``) qubits, renumbered."""
    # Only Z and I can be read off the ancilla-joint Z measurement: X/Y must have
    # been rotated into Z (_hadamard_setup) and a Pauli sum split into terms by the
    # caller. Anything else would be silently read as I.
    if not isinstance(observable, str) or set(observable) - {"Z", "I"}:
        raise ValueError(
            "Fragments are read in the Z basis: expected a Z/I string after the "
            f"basis rotation, got {observable!r}."
        )
    fragments = []
    for frag, qmap, ancillas in _config_fragment_circuits(
        cut_dag, components, assignments
    ):
        z_qubits = list(ancillas)
        for q, ch in enumerate(observable):
            if ch == "Z" and q in qmap:
                z_qubits.append(qmap[q])
        fragments.append((frag, z_qubits))
    return fragments


def _config_fragment_circuits(cut_dag, components, assignments):
    """``[(fragment_circuit, qubit_map, ancillas), ...]`` for one configuration.

    ``qubit_map`` maps each original qubit to its index in the fragment (the data
    qubits are ``0..k-1`` in original order, the ancillas follow).
    """
    return [
        _component_hadamard_circuit(cut_dag.subgraph(nodes).copy(), assignments)
        for nodes in components
    ]


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


def _hadamard_fragment_specs(circuit, gates_cut, observable, joint=False):
    """Enumerate the double-Hadamard-test fragments for a set of cuts.

    Returns a list with one record per non-zero enumeration setting:
    ``(weight, [(fragment_circuit, z_qubits), ...])``. Each fragment is an independent
    qibo circuit carrying its local ancilla(s); ``z_qubits`` are the qubits read out as
    a joint ``Z`` (ancillas + observed qubits). The reconstruction of ``<observable>``
    is ``sum_record weight * prod_fragment <Z_{z_qubits}>`` -- the shared core of the
    all-in-one and three-step *enumeration* drivers.
    """
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut, observable)
    groups = _cut_groups(cut_dag, components, decomps, joint)
    settings_per_group = [_cut_settings(len(terms)) for _, terms in groups]

    specs = []
    for config in itertools.product(*settings_per_group):  # one (s, t, g) per group
        weight = 1.0
        for (_, terms), (i, j, g) in zip(groups, config):
            weight *= terms[i][0] * terms[j][0] * (-1) ** g
        if weight == 0.0:
            continue
        assignments = _config_assignments(groups, config)
        fragments = _config_fragments(
            cut_dag, components, assignments, _z_basis(observable)
        )
        specs.append((weight, fragments))
    return specs


@task(returns=1)
def _evaluate_hadamard_fragment(
    frag_dag, num_qubits, software, shots, qpu, qpu_dict, seed=None
):
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
    return backend.evaluate_expectation(frag_dag, num_qubits, shots=shots, seed=seed)


@task(returns=1, values=COLLECTION_IN)
def _reduce_hadamard(values, weights, config_of):
    """Fold the per-fragment values into ``<observable>`` (the fan-in COMPSs task).

    Takes the whole fan-out (``values`` a ``COLLECTION_IN`` of the fragment futures)
    and returns ``sum_config weight * prod_fragment value`` -- so the task graph is a
    clean fan-out of fragment tasks into this single reduction, and the reduction runs
    off the master.
    """
    return _gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of)


@task(returns=1, partials=COLLECTION_IN)
def _sum_reduce(partials):
    """Sum a batch of partial reductions -- one node of the reduction tree."""
    return sum(p for p in partials if p is not None)


def _reduce_hadamard_tree(values, weights, config_of, chunk):
    """Hierarchical fan-in reduction (returns the final future).

    Funnelling every fragment value into one :func:`_reduce_hadamard` becomes a
    data-transfer + serial-compute bottleneck once the fragment count is large (e.g.
    the sampling driver with many ``num_samples``). Instead, reduce ``chunk``
    configurations at a time with independent partial :func:`_reduce_hadamard` tasks,
    then combine those partial sums through a balanced tree of :func:`_sum_reduce`
    tasks (arity ``chunk``) -- so the reduction itself parallelises and no single task
    sees all values. Configurations partition cleanly (a fragment product never spans
    configs), so the result is identical to the single reduction, just regrouped.
    """
    from collections import defaultdict

    frags_by_config = defaultdict(list)
    for value, cfg in zip(values, config_of):
        frags_by_config[cfg].append(value)

    n_configs = len(weights)
    if n_configs == 0:
        return _reduce_hadamard([], [], [])

    partials = []  # one partial-sum future per chunk of configurations
    for start in range(0, n_configs, chunk):
        chunk_values, chunk_weights, chunk_config_of = [], [], []
        for local, cfg in enumerate(range(start, min(start + chunk, n_configs))):
            for value in frags_by_config[cfg]:
                chunk_values.append(value)
                chunk_config_of.append(local)
            chunk_weights.append(weights[cfg])
        partials.append(_reduce_hadamard(chunk_values, chunk_weights, chunk_config_of))

    # Combine the partial sums up a balanced tree. The arity must be >= 2 or the count
    # never shrinks (chunk == 1 would loop forever).
    arity = max(2, chunk)
    while len(partials) > 1:
        partials = [
            _sum_reduce(partials[i : i + arity]) for i in range(0, len(partials), arity)
        ]
    return partials[0]


def _dispatch_reduce(values, weights, config_of, reduce_chunk, sync):
    """Reduce the fan-out either with a single fan-in task or, when ``reduce_chunk`` is
    set and there are more configurations than that, a balanced reduction tree. Returns
    the value (``sync=True``) or the reduction future (``sync=False``)."""
    if reduce_chunk and len(weights) > reduce_chunk:
        result = _reduce_hadamard_tree(values, weights, config_of, reduce_chunk)
    else:
        result = _reduce_hadamard(values, weights, config_of)
    return compss_wait_on(result) if sync else result


def _gate_cutting_hadamard(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    shots=None,
    seed=None,
    qpu=False,
    qpu_dict=None,
    sync=True,
    reduce_chunk=None,
    joint=False,
    num_samples=None,
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
    :param observable: A Pauli string over the circuit's qubits (I/X/Y/Z; X and Y
        are rotated into Z), or a weighted sum ``[(coef, "PAULI"), ...]``.
    :param software: Execution backend for the sampled path (``"qibo"``, ``"qiskit"``,
        ``"cudaq"``, ``"pennylane"``).
    :param shots: Samples per fragment; ``None`` selects the exact statevector path
        (``software="qibo"`` only).
    :param seed: Unused placeholder kept for signature compatibility.
    :param qpu: Route fragment execution to a QPU via ``build_execution_backend``.
    :param qpu_dict: QPU connection settings forwarded to the backend.
    :param sync: If ``True`` (default) block on the reduction and return the value; if
        ``False`` return the (unresolved) reduction future so several independent cuts
        can be launched back-to-back and awaited together (no barrier between them --
        the caller resolves with ``compss_wait_on``). Ignored on the exact path.
    :param reduce_chunk: If set, reduce the fan-out through a balanced tree that folds
        this many configurations per partial reduction (and combines partials with the
        same arity) instead of a single fan-in task -- keeps the reduction parallel
        when the fragment count is large. ``None`` (default) uses one reduction task.
    :param joint: Cut the gates crossing the same boundary together, with one shared
        ancilla per fragment: ``gamma = 2**(k+1) - 1`` instead of ``3**k`` for k
        CZ/CX, and fewer settings (``2**k (2**(k+1) - 1)`` instead of ``6**k``). The
        shared ancilla must stay coherent from the boundary's first cut gate to its
        last. ``False`` (default) cuts each gate on its own.
    :param num_samples: Draw this many configurations from the quasiprobability PMF
        ``p(i, j, g) = c_i c_j / phi`` and average ``phi * sign * prod <Z>`` (the
        native Harrow-Lowe estimator), instead of enumerating every configuration.
        Bounds the cost by ``num_samples`` whatever the number of cuts, at a
        ``phi**2 / num_samples`` variance. ``seed`` seeds the draws.
    :return: The reconstructed expectation value, or its future when ``sync=False``.
    """
    if num_samples is not None:
        return _gate_cutting_hadamard_montecarlo(
            circuit,
            gates_cut,
            observable,
            software,
            num_samples,
            shots,
            seed,
            qpu,
            qpu_dict,
            sync,
            reduce_chunk,
            joint,
        )
    terms = _observable_terms(observable)
    if not _is_single_term(terms):
        values = [
            _gate_cutting_hadamard(
                circuit,
                gates_cut,
                pauli,
                software,
                shots,
                seed,
                qpu,
                qpu_dict,
                sync,
                reduce_chunk,
                joint,
            )
            for _, pauli in terms
        ]
        coefficients = [c for c, _ in terms]
        if sync:
            return sum(c * v for c, v in zip(coefficients, values))
        return _weighted_sum(values, coefficients)
    observable = terms[0][1]
    exact = shots is None  # statevector, whatever the backend
    specs = _hadamard_fragment_specs(circuit, gates_cut, observable, joint)

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
                        frag_dag,
                        frag.nqubits,
                        software,
                        shots,
                        qpu,
                        qpu_dict,
                        next_seed(seed, len(values)),
                    )
                )
            config_of.append(cfg_idx)

    if exact:  # in-process, no tasks to synchronise
        return _gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
    # Fan the per-fragment futures into the reduction (single task, or a tree when
    # reduce_chunk is set); sync=False returns the future so callers can overlap cuts.
    return _dispatch_reduce(values, weights, config_of, reduce_chunk, sync)


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

        :param values: One value (or future) per subcircuit, in this object's order.
        :return: The reconstructed expectation value of the uncut circuit.
        """
        return _dispatch_reduce(list(values), self.weights, self.config_of, None, True)


def _materialize_subcircuit(frag, z_qubits, software):
    """A fragment qibo circuit as a :class:`~Qdislib.utils.backends.Subcircuit`
    written in ``software``: the qubits it does not read become "Observable I"."""
    from Qdislib.utils.backends import Subcircuit

    obs = _z_observable(frag.nqubits, z_qubits)
    return Subcircuit(circuit_qibo_to_dag(frag, obs_I=obs), frag.nqubits, software)


def _gate_cutting_hadamard_subcircuits(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    joint=False,
    num_samples=None,
    seed=None,
):
    """Build the Hadamard-cut fragments as standalone subcircuits (step one of three).

    Returns a :class:`HadamardCutSubcircuits` -- a list of ``[fragment, observable]``
    subcircuits that the caller evaluates independently (e.g. with
    :func:`expectation_value`, on any backend / QPU, in parallel) and then feeds to
    :meth:`~HadamardCutSubcircuits.reconstruct`. This is the explicit three-step form
    of :func:`gate_cutting_hadamard`.

    :param circuit: A qibo circuit.
    :param gates_cut: List of cut gate names.
    :param observable: A Pauli string over the circuit's qubits (I/X/Y/Z; X and Y
        are rotated into Z), or a weighted sum ``[(coef, "PAULI"), ...]``.
    :param software: Backend the subcircuits are materialised for.
    :param joint: Cut the gates crossing each boundary together (see
        :func:`gate_cutting_hadamard`).
    :param num_samples: Draw this many configurations instead of enumerating them
        (see :func:`gate_cutting_hadamard`); ``seed`` seeds the draws.
    :return: A :class:`HadamardCutSubcircuits`.
    """
    if num_samples is not None:
        return _gate_cutting_hadamard_montecarlo_subcircuits(
            circuit, gates_cut, observable, software, num_samples, seed, joint
        )
    terms = _observable_terms(observable)
    if not _is_single_term(terms):
        from Qdislib.utils.backends import CombinedCutSubcircuits

        parts = [
            _gate_cutting_hadamard_subcircuits(circuit, gates_cut, p, software, joint)
            for _, p in terms
        ]
        return CombinedCutSubcircuits(parts, [c for c, _ in terms])
    observable = terms[0][1]
    specs = _hadamard_fragment_specs(circuit, gates_cut, observable, joint)
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


def _gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of):
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
    ``decomps`` may be a group's terms (:func:`_cut_groups`) as well as a cut's: only
    each term's leading coefficient is read.
    """
    rng = np.random.default_rng(seed)
    phi_total = 1.0
    draws = []  # per cut / group: arrays s, t, g over the samples
    for terms in decomps:
        phi_total *= phi_factor(terms)
        draws.append(_draw_settings(terms, num_samples, rng))
    configs = []
    for n in range(num_samples):
        config = [(int(s[n]), int(t[n]), int(g[n])) for s, t, g in draws]
        configs.append((config, (-1) ** sum(g for _, _, g in config)))
    return configs, phi_total


def _draw_settings(terms, num_samples, rng):
    """``num_samples`` settings ``(s, t, g)`` with probability ``c_s c_t / phi``,
    without listing the settings: ``s`` and ``t`` are drawn ``∝ c`` (gate by gate for
    a joint group) and ``s = t``, which has one setting where ``s != t`` has two, is
    kept half the time."""
    if isinstance(terms, _ProductTerms):
        per_gate = [np.array([abs(c) for c, _, _ in t]) for t in terms.per_gate]
    else:
        per_gate = [np.array([abs(t[0]) for t in terms])]

    def term():
        index = np.zeros(num_samples, dtype=np.int64)
        for c in per_gate:
            index = index * len(c) + rng.choice(len(c), size=num_samples, p=c / c.sum())
        return index

    s, t = term(), term()
    redo = (s == t) & (rng.random(num_samples) < 0.5)
    while redo.any():
        s_new, t_new = term(), term()
        s[redo], t[redo] = s_new[redo], t_new[redo]
        redo = redo & (s == t) & (rng.random(num_samples) < 0.5)
    g = np.where(s == t, 0, rng.integers(0, 2, size=num_samples))
    return s, t, g


def _sampling_weight(phi_total, sign, num_samples):
    """The per-sample coefficient that turns the shared weighted-product reconstruction
    into the Monte-Carlo estimator: with ``weight_s = phi_total * sign_s / N`` the sum
    ``sum_s weight_s * prod_frag value`` equals ``phi_total * mean_s(sign_s * prod
    value)`` -- the unbiased Harrow-Lowe estimator of ``<observable>``."""
    return phi_total * sign / num_samples


def _gate_cutting_hadamard_montecarlo(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    num_samples=1000,
    shots=1024,
    seed=None,
    qpu=False,
    qpu_dict=None,
    sync=True,
    reduce_chunk=None,
    joint=False,
):
    """Estimate ``<observable>`` by **Monte-Carlo sampling** the Hadamard-test cut.

    Reached through ``gate_cutting_hadamard(..., num_samples=N)``.

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
    :param observable: A Pauli string over the circuit's qubits (I/X/Y/Z; X and Y
        are rotated into Z), or a weighted sum ``[(coef, "PAULI"), ...]``.
    :param software: Fragment-execution backend (``"qibo"``/``"qiskit"``/``"cudaq"``/
        ``"pennylane"``).
    :param num_samples: Number of configurations drawn from the PMF.
    :param shots: Samples per fragment; ``None`` uses the exact statevector per
        fragment (``software="qibo"`` only) -- then only the *configuration* sampling
        is stochastic.
    :param seed: Seed for the configuration draws (reproducible).
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :param sync: If ``True`` (default) block on the reduction and return the value; if
        ``False`` return the (unresolved) reduction future so several cuts can run
        concurrently without a barrier between them. Ignored on the exact path.
    :param reduce_chunk: If set, reduce through a balanced tree folding this many
        configurations per partial reduction (recommended for large ``num_samples`` so
        the reduction stays parallel); ``None`` (default) uses one reduction task.
    :param joint: Cut the gates crossing each boundary together (see
        :func:`gate_cutting_hadamard`). Here the variance, ``phi_total**2 /
        num_samples``, is what shrinks: ``(2**(k+1) - 1)**2`` instead of ``9**k``.
    :return: The estimated expectation value, or its future when ``sync=False``.
    """
    terms = _observable_terms(observable)
    if not _is_single_term(terms):
        values = [
            _gate_cutting_hadamard_montecarlo(
                circuit,
                gates_cut,
                pauli,
                software,
                num_samples,
                shots,
                seed,
                qpu,
                qpu_dict,
                sync,
                reduce_chunk,
                joint,
            )
            for _, pauli in terms
        ]
        coefficients = [c for c, _ in terms]
        if sync:
            return sum(c * v for c, v in zip(coefficients, values))
        return _weighted_sum(values, coefficients)
    observable = terms[0][1]
    exact = shots is None  # statevector, whatever the backend
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut, observable)
    groups = _cut_groups(cut_dag, components, decomps, joint)
    configs, phi_total = _sample_configs([t for _, t in groups], num_samples, seed)

    values, config_of, weights = [], [], []
    for s, (config, sign) in enumerate(configs):
        weights.append(_sampling_weight(phi_total, sign, num_samples))
        assignments = _config_assignments(groups, config)
        for frag, z_qubits in _config_fragments(
            cut_dag, components, assignments, _z_basis(observable)
        ):
            if exact:
                values.append(_expectation_z(frag, z_qubits))
            else:
                obs = _z_observable(frag.nqubits, z_qubits)
                frag_dag = circuit_qibo_to_dag(frag, obs_I=obs)
                values.append(
                    _evaluate_hadamard_fragment(
                        frag_dag,
                        frag.nqubits,
                        software,
                        shots,
                        qpu,
                        qpu_dict,
                        next_seed(seed, "run", len(values)),
                    )
                )
            config_of.append(s)

    if exact:  # in-process, no tasks to synchronise
        return _gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
    return _dispatch_reduce(values, weights, config_of, reduce_chunk, sync)


def _gate_cutting_hadamard_montecarlo_subcircuits(
    circuit,
    gates_cut,
    observable,
    software="qibo",
    num_samples=1000,
    seed=None,
    joint=False,
):
    """Monte-Carlo-sampled Hadamard-cut fragments as standalone subcircuits.

    Reached through ``gate_cutting_hadamard_subcircuits(..., num_samples=N)``.

    The sampling counterpart of :func:`gate_cutting_hadamard_subcircuits`: draws
    ``num_samples`` configurations from the PMF and returns their fragments as a
    :class:`HadamardCutSubcircuits` (already carrying the ``phi_total * sign /
    num_samples`` weights), so the usual ``evaluate -> reconstruct`` workflow yields the
    Monte-Carlo estimate::

        subs = qd.gate_cutting_hadamard_subcircuits(c, ["CZ_2"], "ZZ", num_samples=2000)
        values = [qd.expectation_value(s) for s in subs]
        estimate = subs.reconstruct(values)

    :param num_samples: Number of configurations drawn from the PMF.
    :param seed: Seed for the configuration draws.
    :param joint: Cut the gates crossing each boundary together (see
        :func:`gate_cutting_hadamard`).
    :return: A :class:`HadamardCutSubcircuits`.
    """
    terms = _observable_terms(observable)
    if not _is_single_term(terms):
        from Qdislib.utils.backends import CombinedCutSubcircuits

        parts = [
            _gate_cutting_hadamard_montecarlo_subcircuits(
                circuit, gates_cut, p, software, num_samples, seed, joint
            )
            for _, p in terms
        ]
        return CombinedCutSubcircuits(parts, [c for c, _ in terms])
    observable = terms[0][1]
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut, observable)
    groups = _cut_groups(cut_dag, components, decomps, joint)
    configs, phi_total = _sample_configs([t for _, t in groups], num_samples, seed)

    subcircuits, weights, config_of = [], [], []
    for s, (config, sign) in enumerate(configs):
        weights.append(_sampling_weight(phi_total, sign, num_samples))
        assignments = _config_assignments(groups, config)
        for frag, z_qubits in _config_fragments(
            cut_dag, components, assignments, _z_basis(observable)
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


# --------------------------------------------------------------------------- #
# Output distribution
# --------------------------------------------------------------------------- #
#
# Every fragment is read in the Z basis; an observable only chooses which data bits
# join the ancillas' parity. So one run of a fragment gives <Z_ancillas Z^a> for
# EVERY subset a of its data qubits. Weighted and multiplied across fragments
# exactly as for an expectation value, these give every <Z^A> of the full circuit,
# and the inverse Walsh-Hadamard transform P(x) = 2**-N sum_A (-1)**(A.x) <Z^A>
# turns them into the output distribution -- the Hadamard-cut counterpart of
# gate_cutting_sampling / wire_cutting_sampling.


def _fragment_correlators(
    frag_dag, nqubits, num_data, software, shots, seed, qpu, qpu_dict
):
    """``<Z_ancillas Z^a>`` for every subset ``a`` of a fragment's data qubits.

    The fragment's data qubits are ``0..num_data-1`` and its ancillas follow. With
    ``shots=None`` (``software="qibo"``) the values are exact (statevector);
    otherwise they come from ``shots`` samples on the chosen backend.

    :return: An array indexed by the data-qubit mask (bit ``i`` = data qubit ``i``).
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard

    probs = np.zeros(2**nqubits)  # index bit q = fragment qubit q
    if shots is None:
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo

        built = dag_to_circuit_qibo(frag_dag, nqubits)
        circuit = built[0] if isinstance(built, (list, tuple)) else built
        state = np.asarray(circuit().state())
        for index, amplitude in enumerate(state):  # qibo index: qubit 0 is the MSB
            x = sum(((index >> (nqubits - 1 - q)) & 1) << q for q in range(nqubits))
            probs[x] = abs(amplitude) ** 2
    else:
        from Qdislib.core.backend.registry import build_execution_backend

        backend = build_execution_backend(
            software, num_qubits=nqubits, qpu=bool(qpu), qpu_dict=qpu_dict
        )
        counts = backend.evaluate_samples(frag_dag, nqubits, shots, seed=seed)
        total = sum(counts.values()) or 1
        for bitstring, count in counts.items():  # qubit-0-first
            probs[sum(int(b) << q for q, b in enumerate(bitstring))] += count / total
    ancillas = sum(1 << q for q in range(num_data, nqubits))
    signs = np.array([(-1) ** bin(x & ancillas).count("1") for x in range(2**nqubits)])
    return _fast_walsh_hadamard(probs * signs)[: 2**num_data]


@task(returns=1)
def _fragment_correlators_task(
    frag_dag, nqubits, num_data, original, software, shots, seed, qpu, qpu_dict
):
    """:func:`_fragment_correlators` as a COMPSs task; returns ``(corr, original)``."""
    corr = _fragment_correlators(
        frag_dag, nqubits, num_data, software, shots, seed, qpu, qpu_dict
    )
    return corr, original


def _hadamard_distribution_units(circuit, gates_cut, joint, num_samples, seed):
    """``(weights, [(config, frag_dag, nqubits, num_data, original), ...])`` for the
    distribution: every configuration (or ``num_samples`` drawn from the PMF)."""
    cut_dag, components, decomps = _hadamard_setup(circuit, gates_cut)
    groups = _cut_groups(cut_dag, components, decomps, joint)
    if num_samples is None:
        configs = []
        for config in itertools.product(*(_cut_settings(len(t)) for _, t in groups)):
            weight = 1.0
            for (_, terms), (i, j, g) in zip(groups, config):
                weight *= terms[i][0] * terms[j][0] * (-1) ** g
            if weight != 0.0:
                configs.append((config, weight))
    else:
        drawn, phi_total = _sample_configs([t for _, t in groups], num_samples, seed)
        configs = [
            (config, _sampling_weight(phi_total, sign, num_samples))
            for config, sign in drawn
        ]
    weights, units = [], []
    for index, (config, weight) in enumerate(configs):
        weights.append(weight)
        assignments = _config_assignments(groups, config)
        for frag, qmap, _ancillas in _config_fragment_circuits(
            cut_dag, components, assignments
        ):
            original = sorted(qmap, key=qmap.get)
            units.append(
                (
                    index,
                    circuit_qibo_to_dag(frag),
                    frag.nqubits,
                    len(original),
                    original,
                )
            )
    return weights, units


def _num_circuit_qubits(circuit):
    from Qdislib.utils.backends import convert_to_dag, num_qubits

    try:
        return num_qubits(circuit)
    except TypeError:  # a DAG
        dag = convert_to_dag(circuit)
        return 1 + max(q for _, d in dag.nodes(data=True) for q in d.get("qubits", ()))


def _gate_cutting_hadamard_sampling_reconstruction(
    sampled, weights, config_of, num_qubits
):
    """Recombine per-fragment correlators into the output distribution.

    :param sampled: One ``(corr, original)`` per fragment (from
        :func:`~Qdislib.sample_subcircuit` or the all-in-one driver), in order.
    :param weights: Signed coefficient of each configuration.
    :param config_of: Configuration index of each fragment (parallel to ``sampled``).
    :param num_qubits: Width of the original circuit.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings), the signed
        estimate abs-normalised like :func:`gate_cutting_sampling`.
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard

    sampled = compss_wait_on(sampled)
    dim = 2**num_qubits
    masks = np.arange(dim)
    contributions = [np.full(dim, float(w)) for w in weights]
    for (corr, original), cfg in zip(sampled, config_of):
        local = np.zeros(dim, dtype=np.int64)
        for i, qubit in enumerate(original):
            local |= ((masks >> qubit) & 1) << i
        contributions[cfg] = contributions[cfg] * np.asarray(corr)[local]
    z_correlators = np.sum(contributions, axis=0) if contributions else np.ones(dim)
    probabilities = _fast_walsh_hadamard(z_correlators) / dim
    keys = ["".join(str((x >> q) & 1) for q in range(num_qubits)) for x in range(dim)]
    norm = float(np.abs(probabilities).sum()) or 1.0
    return {key: abs(float(p)) / norm for key, p in zip(keys, probabilities)}


def _gate_cutting_hadamard_sampling(
    circuit,
    gates_cut,
    software="qibo",
    shots=None,
    seed=None,
    qpu=False,
    qpu_dict=None,
    joint=False,
    num_samples=None,
    observable=None,
    sync=True,
):
    """Reconstruct the output **probability distribution** with the Hadamard cut.

    The Hadamard-cut counterpart of :func:`gate_cutting_sampling` /
    :func:`wire_cutting_sampling`. Each fragment of each configuration is run once
    and read in the Z basis, which gives its correlators for every subset of its
    output qubits; their weighted products give every ``<Z^A>`` of the circuit, and
    the inverse Walsh-Hadamard transform gives ``P(x)``. Same circuits and cost as
    one :func:`gate_cutting_hadamard` call. Fragments run as parallel PyCOMPSs tasks.

    :param circuit: A qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: Two-qubit gate names to cut.
    :param software: Fragment backend (``"qibo"``/``"qiskit"``/``"cudaq"``/
        ``"pennylane"``).
    :param shots: Samples per fragment; ``None`` gives the exact distribution
        (``software="qibo"``, statevector).
    :param seed: Seed for the fragment sampling (``seed + i`` per fragment) and for
        the configuration draws when ``num_samples`` is set.
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :param joint: Cut the gates crossing each boundary together (see
        :func:`gate_cutting_hadamard`).
    :param num_samples: Draw this many configurations from the quasiprobability
        PMF instead of enumerating them all (as ``gate_cutting_hadamard(...,
        num_samples=N)`` does for an expectation value).
    :param observable: Not accepted: this returns a distribution. (Up to 1.0.1 this
        name estimated an expectation value; that is now
        ``gate_cutting_hadamard(..., num_samples=N)``.)
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    if observable is not None:
        raise TypeError(
            "gate_cutting_hadamard_sampling* now returns the output distribution and "
            "takes no observable (it is the Hadamard counterpart of "
            "gate_cutting_sampling). For the Monte-Carlo estimate of an expectation "
            "value use gate_cutting_hadamard(circuit, gates_cut, observable, "
            "num_samples=N) or gate_cutting_hadamard_subcircuits(..., num_samples=N)."
        )
    exact = shots is None  # statevector (qibo), whatever the backend
    weights, units = _hadamard_distribution_units(
        circuit, gates_cut, joint, num_samples, seed
    )
    sampled, config_of = [], []
    for n, (cfg, frag_dag, nqubits, num_data, original) in enumerate(units):
        unit_seed = next_seed(seed, n)
        if exact:
            corr = _fragment_correlators(
                frag_dag, nqubits, num_data, software, None, None, qpu, qpu_dict
            )
            sampled.append((corr, original))
        else:
            sampled.append(
                _fragment_correlators_task(
                    frag_dag,
                    nqubits,
                    num_data,
                    original,
                    software,
                    shots,
                    unit_seed,
                    qpu,
                    qpu_dict,
                )
            )
        config_of.append(cfg)
    if sync or exact:
        return _gate_cutting_hadamard_sampling_reconstruction(
            sampled, weights, config_of, _num_circuit_qubits(circuit)
        )
    return _hadamard_sampling_reconstruction_task(
        sampled, weights, config_of, _num_circuit_qubits(circuit)
    )


@task(returns=1, sampled=COLLECTION_IN)
def _hadamard_sampling_reconstruction_task(sampled, weights, config_of, num_qubits):
    """:func:`gate_cutting_hadamard_sampling_reconstruction` as a task (sync=False)."""
    return _gate_cutting_hadamard_sampling_reconstruction(
        sampled, weights, config_of, num_qubits
    )


class _HadamardSampleUnit(Subcircuit):
    """One Hadamard-cut fragment to sample for the distribution (three-step):
    :func:`~Qdislib.sample_subcircuit` returns its ``(corr, original)`` correlators,
    the ancillas' parity folded in. A renumbered DAG plus metadata (picklable)."""

    kind = "distribution"

    def __init__(self, frag_dag, nqubits, num_data, original, software, qpu, qpu_dict):
        super().__init__(frag_dag, nqubits, software)
        self.num_data = num_data
        self.original = original
        self.qpu = qpu
        self.qpu_dict = qpu_dict

    @property
    def frag_dag(self):
        return self.dag

    @property
    def nqubits(self):
        return self.num_qubits

    @property
    def observable(self):
        return None

    def _expectation(self, shots, backend, seed, cache, cache_mode):
        raise TypeError(
            "This fragment comes from a *_sampling_subcircuits list (distribution): "
            "run it with sample_subcircuit()."
        )

    def _sample(self, shots, backend, seed):
        if shots is None:
            corr = _fragment_correlators(
                self.dag,
                self.num_qubits,
                self.num_data,
                "qibo",
                None,
                None,
                False,
                None,
            )
            return corr, self.original
        return _fragment_correlators_task(
            self.dag,
            self.num_qubits,
            self.num_data,
            self.original,
            self.software if backend is None else backend,
            shots,
            seed,
            self.qpu,
            self.qpu_dict,
        )

    def sample(self, shots=1024, seed=None):
        """``(corr, original)`` (the 1.x method; see :func:`~Qdislib.sample_subcircuit`)."""
        return self._sample(shots, None, seed)


class HadamardSampledCutSubcircuits(list):
    """Hadamard-cut fragments to sample, plus what reconstructs the distribution::

    subs = qd.gate_cutting_hadamard_sampling_subcircuits(circuit, ["CZ_3"])
    sampled = [qd.sample_subcircuit(s, shots=8192) for s in subs]
    distribution = subs.reconstruct(sampled)
    """

    def __init__(self, units, *, weights, config_of, num_qubits, number_cuts):
        super().__init__(units)
        self.weights = weights
        self.config_of = config_of
        self.num_qubits = num_qubits
        self.number_cuts = number_cuts

    def reconstruct(self, sampled_results):
        """``{bitstring: probability}`` from one sample result per subcircuit."""
        return _gate_cutting_hadamard_sampling_reconstruction(
            sampled_results, self.weights, self.config_of, self.num_qubits
        )


def _gate_cutting_hadamard_sampling_subcircuits(
    circuit,
    gates_cut,
    software="qibo",
    qpu=False,
    qpu_dict=None,
    joint=False,
    num_samples=None,
    seed=None,
    observable=None,
):
    """The Hadamard-cut fragments for the output distribution (three-step workflow).

    Sample each with :func:`~Qdislib.sample_subcircuit` (any backend, your loop),
    then call ``reconstruct`` on the list of results. Arguments as in
    :func:`gate_cutting_hadamard_sampling`.

    :return: A :class:`HadamardSampledCutSubcircuits`.
    """
    if observable is not None:
        raise TypeError(
            "gate_cutting_hadamard_sampling* now returns the output distribution and "
            "takes no observable (it is the Hadamard counterpart of "
            "gate_cutting_sampling). For the Monte-Carlo estimate of an expectation "
            "value use gate_cutting_hadamard(circuit, gates_cut, observable, "
            "num_samples=N) or gate_cutting_hadamard_subcircuits(..., num_samples=N)."
        )
    weights, units = _hadamard_distribution_units(
        circuit, gates_cut, joint, num_samples, seed
    )
    return HadamardSampledCutSubcircuits(
        [
            _HadamardSampleUnit(dag, n, k, original, software, qpu, qpu_dict)
            for _cfg, dag, n, k, original in units
        ],
        weights=weights,
        config_of=[cfg for cfg, *_ in units],
        num_qubits=_num_circuit_qubits(circuit),
        number_cuts=len(gates_cut),
    )


# --------------------------------------------------------------------------- #
# Public API (the calling convention shared by every cut type, see _api.py)
# --------------------------------------------------------------------------- #


@standard_arguments
def gate_cutting_hadamard(
    circuit,
    cuts,
    observable=None,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    joint=True,
    num_samples=None,
    reduce_chunk=None,
    _legacy=None,
):
    """Reconstruct ``<observable>`` with the ancilla / double-Hadamard-test gate cut
    (Harrow & Lowe, arXiv:2403.01018): any two-qubit gate, ``gamma = 3`` per CZ/CX,
    no mid-circuit measurement.

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: Two-qubit gate names to cut.
    :param observable: A Pauli string or a weighted Pauli sum; ``None`` is all ``Z``.
    :param backend: Where the fragments run (see :func:`gate_cutting`).
    :param shots: Shots per fragment; ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling and for the ``num_samples`` draws.
    :param sync: ``False`` returns the value as a future (PyCOMPSs).
    :param joint: Cut the gates crossing the same boundary together, with one shared
        ancilla per fragment: ``gamma = 2**(k+1) - 1`` instead of ``3**k`` for k
        CZ/CX. The shared ancilla must stay coherent from the boundary's first cut
        gate to its last; ``False`` cuts each gate on its own.
    :param num_samples: Draw this many configurations from the quasiprobability
        distribution instead of running them all.
    :param reduce_chunk: Deprecated, removed in the next release: the reconstruction
        is one cheap task. Still folds a shot run's results through a balanced
        tree of this arity.
    :return: The reconstructed expectation value.
    """
    if reduce_chunk is not None:
        import warnings

        warnings.warn(
            "reduce_chunk is deprecated and removed in the next release: the "
            "reconstruction is one task, and a cheap one.",
            DeprecationWarning,
            stacklevel=3,
        )
    if shots is None or num_samples is not None:
        # Exact and Monte-Carlo runs go through circuit_cutting, as gate_cutting's
        # do: every fragment circuit is built and run in its own task.
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _route_to_circuit_cutting,
        )
        from Qdislib.core.cutting_algorithms.mixed_cutting import circuit_cutting

        _route_to_circuit_cutting(_legacy)
        return circuit_cutting(
            circuit,
            {"hadamard": cuts},
            observable,
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            joint=joint,
            num_samples=num_samples,
        )
    legacy = _legacy or {}
    return _gate_cutting_hadamard(
        circuit,
        cuts,
        default_observable(observable, circuit),
        software=resolve_backend(backend, circuit),
        shots=shots,
        seed=seed,
        qpu=legacy.get("qpu", False),
        qpu_dict=legacy.get("qpu_dict"),
        sync=sync,
        reduce_chunk=reduce_chunk,
        joint=joint,
        num_samples=num_samples,
    )


@standard_arguments
def gate_cutting_hadamard_sampling(
    circuit,
    cuts,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    joint=True,
    num_samples=None,
    _legacy=None,
):
    """Reconstruct the output **probability distribution** with the Hadamard cut.
    Arguments as in :func:`gate_cutting_hadamard`.

    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings).
    """
    if shots is None or num_samples is not None:
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _route_to_circuit_cutting,
        )
        from Qdislib.core.cutting_algorithms.mixed_cutting import (
            circuit_cutting_sampling,
        )

        _route_to_circuit_cutting(_legacy)
        return circuit_cutting_sampling(
            circuit,
            {"hadamard": cuts},
            backend=backend,
            shots=shots,
            seed=seed,
            sync=sync,
            joint=joint,
            num_samples=num_samples,
        )
    legacy = _legacy or {}
    return _gate_cutting_hadamard_sampling(
        circuit,
        cuts,
        software=resolve_backend(backend, circuit),
        shots=shots,
        seed=seed,
        qpu=legacy.get("qpu", False),
        qpu_dict=legacy.get("qpu_dict"),
        joint=joint,
        num_samples=num_samples,
        sync=sync,
    )


@standard_arguments
def gate_cutting_hadamard_subcircuits(
    circuit,
    cuts,
    observable=None,
    *,
    software=None,
    joint=True,
    num_samples=None,
    seed=None,
):
    """The fragments of :func:`gate_cutting_hadamard`, to run yourself (three-step).

    :param software: The SDK the fragments are written in (``None``: the
        circuit's own); where they run is chosen when you run them.
    :param seed: Seed for the ``num_samples`` draws.
    :return: A list of subcircuits with ``reconstruct(values)``.
    """
    return _gate_cutting_hadamard_subcircuits(
        circuit,
        cuts,
        default_observable(observable, circuit),
        software=resolve_software(software, circuit),
        joint=joint,
        num_samples=num_samples,
        seed=seed,
    )


@standard_arguments
def gate_cutting_hadamard_sampling_subcircuits(
    circuit, cuts, *, software=None, joint=True, num_samples=None, seed=None
):
    """The fragments of :func:`gate_cutting_hadamard_sampling`, to sample yourself.

    :return: A list of subcircuits with ``reconstruct(results)``.
    """
    return _gate_cutting_hadamard_sampling_subcircuits(
        circuit,
        cuts,
        software=resolve_software(software, circuit),
        joint=joint,
        num_samples=num_samples,
        seed=seed,
    )


# Deprecated in 2.0.0: ``subcircuits.reconstruct(results)`` is the one way.


def gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of):
    """Deprecated: use ``subcircuits.reconstruct(results)``."""
    import warnings

    warnings.warn(
        "gate_cutting_hadamard_subcircuit_reconstruction() is deprecated: use subcircuits.reconstruct(results); it is "
        "removed in the next release.",
        DeprecationWarning,
        stacklevel=2,
    )
    return _gate_cutting_hadamard_subcircuit_reconstruction(values, weights, config_of)
