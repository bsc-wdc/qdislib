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

"""A small numpy statevector simulator for Qdislib DAGs (exact mode).

The DAG's gates are OpenQASM's (:mod:`Qdislib.utils.gates`), so their matrices are
fixed here and applied straight to the state, with no SDK circuit in between. For the
small fragment circuits of a cut this is many times faster than building and
simulating an SDK circuit, whose per-gate overhead dominates at that size.
"""

import cmath
import math

import networkx as nx
import numpy as np

from Qdislib.utils.gates import canonical

__all__ = ["dag_probabilities", "gate_unitary"]

_I2 = np.eye(2, dtype=complex)


def _u(theta, phi, lam):
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array(
        [
            [c, -cmath.exp(1j * lam) * s],
            [cmath.exp(1j * phi) * s, cmath.exp(1j * (phi + lam)) * c],
        ]
    )


def _controlled(u):
    out = np.eye(4, dtype=complex)
    out[2:, 2:] = u
    return out


def _two_pauli(theta, pauli):
    return math.cos(theta / 2) * np.eye(4) - 1j * math.sin(theta / 2) * np.kron(
        pauli, pauli
    )


_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]])
_Z = np.diag([1.0 + 0j, -1.0])
_H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)

_FIXED = {
    "id": _I2,
    "x": _X,
    "y": _Y,
    "z": _Z,
    "h": _H,
    "s": np.diag([1, 1j]),
    "sdg": np.diag([1, -1j]),
    "t": np.diag([1, cmath.exp(1j * math.pi / 4)]),
    "tdg": np.diag([1, cmath.exp(-1j * math.pi / 4)]),
    "sx": np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2,
    "cx": _controlled(_X),
    "cy": _controlled(_Y),
    "cz": _controlled(_Z),
    "ch": _controlled(_H),
    "swap": np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex),
    "iswap": np.array(
        [[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]], dtype=complex
    ),
}
_CCX = np.eye(8, dtype=complex)
_CCX[6:, 6:] = _X
_CSWAP = np.eye(8, dtype=complex)
_CSWAP[[5, 6], :] = _CSWAP[[6, 5], :]
_FIXED["ccx"] = _CCX
_FIXED["cswap"] = _CSWAP


def gate_unitary(name, params=()):
    """The matrix of a gate (any vocabulary, see :func:`canonical`), its first qubit
    the most significant; ``None`` if this simulator does not know it."""
    name, params = canonical(name, params)
    if name in _FIXED:
        return _FIXED[name]
    if name == "unitary":  # a dense gate: its matrix, first qubit most significant
        return np.asarray(params[0], dtype=complex)
    p = [float(v) for v in params]
    if name == "rx":
        return _u(p[0], -math.pi / 2, math.pi / 2)
    if name == "ry":
        return _u(p[0], 0.0, 0.0)
    if name == "rz":
        return np.diag([cmath.exp(-0.5j * p[0]), cmath.exp(0.5j * p[0])])
    if name == "p":
        return np.diag([1, cmath.exp(1j * p[0])])
    if name == "u":
        return _u(*p[:3])
    if name == "cp":
        return np.diag([1, 1, 1, cmath.exp(1j * p[0])])
    if name in ("crx", "cry", "crz"):
        return _controlled(gate_unitary(name[1:], p))
    if name == "cu":
        return _controlled(cmath.exp(1j * p[3]) * _u(*p[:3]))
    if name == "rzz":
        e = cmath.exp(-0.5j * p[0])
        return np.diag([e, e.conjugate(), e.conjugate(), e])
    if name == "rxx":
        return _two_pauli(p[0], _X)
    if name == "ryy":
        return _two_pauli(p[0], _Y)
    return None


def dag_probabilities(dag, num_qubits):
    """Basis-state probabilities of a DAG on qubits ``0..num_qubits-1``, index bit
    ``q`` = qubit ``q``; ``None`` if a gate is unknown here (or a measurement is
    left), so the caller can use an SDK instead."""
    state = np.zeros((2,) * num_qubits, dtype=complex)
    state[(0,) * num_qubits] = 1.0
    for node in nx.topological_sort(dag):
        data = dag.nodes[node]
        name = str(data.get("gate") or "")
        if name.upper() == "OBSERVABLE I":
            continue
        matrix = gate_unitary(name, data.get("parameters") or ())
        if matrix is None:
            return None
        qubits = list(data["qubits"])
        k = len(qubits)
        tensor = np.asarray(matrix).reshape((2,) * (2 * k))
        state = np.tensordot(tensor, state, axes=(list(range(k, 2 * k)), qubits))
        state = np.moveaxis(state, list(range(k)), qubits)
    # axis q is qubit q; index bit q = qubit q means reversing the axes
    probs = np.abs(state) ** 2
    return probs.transpose(list(range(num_qubits))[::-1]).reshape(-1)
