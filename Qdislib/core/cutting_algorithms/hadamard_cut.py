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
]

# Single-qubit Paulis (the building blocks of the local unitaries below).
_I = np.eye(2, dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Two-qubit gates whose only non-local part is a ZZ interaction, so they cut into
# exactly two unitary terms with the optimal gamma. CX/CNOT reduce to CZ by
# conjugating the target with H (a local gate absorbed into the fragment).
SUPPORTED_HADAMARD_GATES = ("CZ", "CP", "CPHASE", "RZZ", "CX", "CNOT")


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
    phi = np.angle(np.asarray(diag, dtype=complex))  # phi_00, phi_01, phi_10, phi_11
    p00, p01, p10, p11 = phi
    alpha = (p00 + p01 + p10 + p11) / 4  # II (global phase)
    beta = (p00 - p01 + p10 - p11) / 4  # IZ (local on B)
    gamma = (p00 + p01 - p10 - p11) / 4  # ZI (local on A)
    delta = (p00 - p01 - p10 + p11) / 4  # ZZ (the only non-local part)

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
    elif name in ("CP", "CPHASE"):
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
