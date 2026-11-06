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

"""One gate vocabulary for every Qdislib DAG: OpenQASM 3's standard gates.

Every SDK names and parameterises gates its own way (Qiskit ``cx``/``u``/``p``,
qibo ``CNOT``/``U3``/``U1``, CUDA-Q ``x.ctrl``/``r1``, PennyLane ``PauliX``/
``PhaseShift``). A circuit's DAG stores its gates in one vocabulary -- the gate names
and definitions of OpenQASM 3's ``stdgates.inc``, plus ``rxx``, ``ryy``, ``rzz`` and
``iswap`` -- so the cutting code, the cut finder and every converter speak one
language. :func:`canonical` converts any SDK's gate to it, parameters included (qibo's
``CU3`` is a controlled SU(2) gate, so it becomes ``cu`` with the phase that makes it
so); :func:`lowered` turns a canonical gate back into the upper-case names the qibo,
CUDA-Q and PennyLane builders use. Qiskit's names are the OpenQASM ones.

Measurements, barriers and Qdislib's own markers (``Observable I``, cut markers) are
not gates and pass through unchanged.
"""

import math
import re

__all__ = ["GATES", "canonical", "lowered", "resolve_cuts"]

#: The canonical gates: ``name -> (qubits, parameters)``. OpenQASM 3 ``stdgates.inc``
#: (``u`` is OpenQASM's ``U``; ``cu`` has its fourth, phase parameter) plus the
#: two-qubit rotations ``rxx``, ``ryy``, ``rzz`` and ``iswap``.
GATES = {
    "id": (1, 0),
    "x": (1, 0),
    "y": (1, 0),
    "z": (1, 0),
    "h": (1, 0),
    "s": (1, 0),
    "sdg": (1, 0),
    "t": (1, 0),
    "tdg": (1, 0),
    "sx": (1, 0),
    "rx": (1, 1),
    "ry": (1, 1),
    "rz": (1, 1),
    "p": (1, 1),
    "u": (1, 3),
    "cx": (2, 0),
    "cy": (2, 0),
    "cz": (2, 0),
    "ch": (2, 0),
    "cp": (2, 1),
    "crx": (2, 1),
    "cry": (2, 1),
    "crz": (2, 1),
    "cu": (2, 4),
    "swap": (2, 0),
    "iswap": (2, 0),
    "rxx": (2, 1),
    "ryy": (2, 1),
    "rzz": (2, 1),
    "ccx": (3, 0),
    "cswap": (3, 0),
}

# SDKs whose U3 / CU3 are OpenQASM's U and controlled U.
_OPENQASM_SDKS = ("qiskit", "qilisdk")

# Not gates: left as they are.
_NOT_GATES = {"M", "MEASURE", "BARRIER", "OBSERVABLE I", "SUBS", "RESET"}

# Other names for a canonical gate with the same parameters.
_RENAMES = {
    "I": "id",
    "ID": "id",
    "IDENTITY": "id",
    "CNOT": "cx",
    "TOFFOLI": "ccx",
    "FREDKIN": "cswap",
    "HADAMARD": "h",
    "PAULIX": "x",
    "PAULIY": "y",
    "PAULIZ": "z",
    "ADJOINT(S)": "sdg",
    "ADJOINT(T)": "tdg",
    "SDAG": "sdg",
    "TDAG": "tdg",
    "U1": "p",
    "R1": "p",
    "PHASE": "p",
    "PHASESHIFT": "p",
    "CU1": "cp",
    "CPHASE": "cp",
    "CR1": "cp",
    "CONTROLLEDPHASESHIFT": "cp",
    "U3": "u",  # = U up to a global phase (qibo's U3 is the SU(2) form)
    "ISWAP": "iswap",
    "SWAP": "swap",
}


def reverse_qubit_order(matrix):
    """A gate's matrix with its qubits' significance reversed: Qiskit's (first
    qubit least significant) to Qdislib's DAG convention (first qubit most
    significant), and back. A dense ``unitary`` gate keeps its matrix in the DAG's
    convention, so every SDK can read it."""
    import numpy as np

    matrix = np.asarray(matrix)
    n = int(round(math.log2(matrix.shape[0])))
    axes = list(range(n))[::-1]
    tensor = matrix.reshape((2,) * (2 * n))
    return tensor.transpose(axes + [n + a for a in axes]).reshape(matrix.shape)


def canonical(name, params=(), sdk=None):
    """``(name, params)`` of a gate in the canonical vocabulary.

    :param name: The gate's name in any SDK (any case).
    :param params: Its parameters.
    :param sdk: ``"qiskit"`` and ``"qilisdk"`` read ``CU3`` with its OpenQASM meaning
        (a controlled ``U``); otherwise ``CU3`` is qibo's controlled SU(2) gate, which
        is what Qdislib's own cut gadgets emit.
    :return: The canonical name (lower case) and parameters; a name that is not a
        gate (a measurement, barrier or marker) unchanged; an unknown gate lower-cased.
    """
    raw = str(name)
    upper = raw.upper()
    params = tuple(params or ())
    if upper in _NOT_GATES:
        return raw, params
    lower = raw.lower()
    if lower in GATES:
        return lower, params
    if upper == "U2":  # u2(phi, lam) = U(pi/2, phi, lam)
        return "u", (math.pi / 2,) + params
    if upper == "CU3":
        if not params:  # only the name is asked for (resolve_cuts)
            return "cu", params
        theta, phi, lam = (float(p) for p in params[:3])
        if sdk in _OPENQASM_SDKS:
            return "cu", (theta, phi, lam, 0.0)
        # qibo's CU3 controls the SU(2) U3 = exp(-i(phi+lam)/2) U
        return "cu", (theta, phi, lam, -(phi + lam) / 2)
    if upper in _RENAMES:
        return _RENAMES[upper], params
    return lower, params


# canonical -> the upper-case names the qibo / CUDA-Q / PennyLane builders take
_LOWER = {
    "id": "I",
    "cx": "CNOT",
    "ccx": "TOFFOLI",
    "p": "U1",
    "cp": "CU1",
    "u": "U3",
    "iswap": "ISWAP",
}


def lowered(name, params=()):
    """A canonical gate as ``[(NAME, qubit positions, params), ...]`` in the
    upper-case names the qibo, CUDA-Q and PennyLane builders use (``CNOT``, ``U3``,
    ``U1``, ``CU1``, ``CU3``...). ``cu`` becomes ``CU3`` plus, when its phase is not
    the SU(2) one, a ``U1`` on the control. Not-gates come back unchanged."""
    name, params = canonical(name, params)
    if name.upper() in _NOT_GATES:
        return [(name, None, params)]
    if name == "unitary":
        raise NotImplementedError(
            "a dense unitary gate runs on Qiskit (Aer), qibo and Qdislib's exact "
            "simulator, not on this SDK yet: decompose it first"
        )
    if name == "cu":
        theta, phi, lam, gamma = (float(p) for p in params)
        ops = [("CU3", None, (theta, phi, lam))]
        phase = gamma + (phi + lam) / 2
        if abs(math.remainder(phase, 2 * math.pi)) > 1e-12:
            ops.append(("U1", (0,), (phase,)))
        return ops
    return [(_LOWER.get(name, name.upper()), None, params)]


_INDEXED = re.compile(r"^(.+)_(\d+)$")


def resolve_cuts(circuit, cuts):
    """The cut names as the circuit's DAG names them, accepting any SDK's name for
    the gate: ``CX_7``, ``CNOT_7`` and ``cx_7`` all name gate 7 when it is a CX.

    A name already in the circuit is kept. A name whose gate does not match gate
    ``N`` raises a ``ValueError`` naming what gate ``N`` is. Works on a list of gate
    names, ``(source, target)`` wire edges, or a ``{"gate", "hadamard", "wire"}``
    dict.
    """
    if not cuts:
        return cuts
    names = _gate_names(circuit)
    if names is None:
        return cuts

    def one(cut):
        if not isinstance(cut, str) or cut in names["ids"]:
            return cut
        match = _INDEXED.match(cut)
        if not match:
            return cut
        index = int(match.group(2))
        node = names["by_index"].get(index)
        if node is None:
            return cut
        want = canonical(match.group(1))[0]
        have = canonical(_INDEXED.match(node).group(1), sdk=names["sdk"])[0]
        if want != have:
            raise ValueError(
                f"No gate {cut!r}: gate {index} of this circuit is {node!r} ({have})."
            )
        return node

    def walk(cut):
        if isinstance(cut, str):
            return one(cut)
        if isinstance(cut, tuple):
            return tuple(walk(c) for c in cut)
        if isinstance(cut, list):
            return [walk(c) for c in cut]
        return cut

    if isinstance(cuts, dict):
        # a string where a list belongs is left for the caller to reject (list()
        # would split it into characters)
        return {
            kind: value if isinstance(value, str) else walk(list(value or []))
            for kind, value in cuts.items()
        }
    return walk(cuts) if isinstance(cuts, (list, tuple, str)) else cuts


def _gate_names(circuit):
    """``{"ids": set, "by_index": {index: node id}, "sdk": ...}`` for a circuit or
    DAG, or ``None`` when its gates cannot be listed."""
    try:
        import networkx as nx

        if isinstance(circuit, nx.Graph):
            ids = set(circuit.nodes)
            sdk = None
        else:
            from Qdislib.utils.backends import iter_gates
            from Qdislib.utils.sdks import sdk_of

            ids = {
                f"{name}_{index}"
                for index, (name, qubits) in enumerate(iter_gates(circuit), 1)
                if qubits
            }
            sdk = sdk_of(circuit).name
    except Exception:  # noqa: BLE001 -- an unknown circuit type: leave names alone
        return None
    by_index = {}
    for node in ids:
        match = _INDEXED.match(str(node))
        if match:
            by_index.setdefault(int(match.group(2)), str(node))
    return {"ids": ids, "by_index": by_index, "sdk": sdk}
