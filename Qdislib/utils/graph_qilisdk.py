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

"""QiliSDK (Qilimanjaro) circuits: to and from a Qdislib DAG, and sampling.

QiliSDK's gates follow OpenQASM's definitions (``U3`` is OpenQASM's ``U``, ``U1``
its ``p``), so converting is renaming: ``CNOT`` is ``cx``, ``Adjoint(S)`` (``S†``)
is ``sdg``, ``Controlled(..., basic_gate=G)`` is ``c`` + ``G``. Gates QiliSDK has no
class for (``rzz``, ``iswap``, ``sx``, ...) are built from ones it has. Samples are
qubit-0-first bitstrings, Qdislib's convention.
"""

import math

from Qdislib.utils.circuit_dag import ops_to_dag
from Qdislib.utils.gates import canonical

__all__ = [
    "circuit_qilisdk_to_dag",
    "dag_to_circuit_qilisdk",
    "qilisdk_gates",
    "sample_qilisdk",
]


def _label(gate):
    """A QiliSDK gate's name, as a node label (``S†`` -> ``SDG``)."""
    return gate.name.upper().replace("†", "DG")


def _parameters(gate):
    try:
        return tuple(float(v) for v in gate.get_parameter_values())
    except Exception:  # noqa: BLE001 -- a gate without parameters
        return ()


def qilisdk_gates(circuit):
    """``[(label, qubits, is_measurement), ...]`` per gate, in circuit order."""
    return [
        (_label(gate), tuple(gate.qubits), gate.name == "M") for gate in circuit.gates
    ]


def circuit_qilisdk_to_dag(circuit, obs_I=None):
    """A Qdislib DAG from a QiliSDK ``Circuit`` (a final measurement is left out,
    as for every SDK)."""
    from Qdislib.utils.graph_qibo import _final_measurements

    final = _final_measurements(
        [(gate.name == "M", list(gate.qubits)) for gate in circuit.gates]
    )
    ops = [
        (
            f"{_label(gate)}_{index}",
            "sdg" if gate.name == "S†" else "tdg" if gate.name == "T†" else gate.name,
            tuple(gate.qubits),
            _parameters(gate),
        )
        for index, gate in enumerate(circuit.gates, start=1)
        if index - 1 not in final
    ]
    return ops_to_dag(ops, obs_I, sdk="qilisdk")


def _emit(circuit, name, qubits, params):
    """Add one canonical gate to a QiliSDK circuit (decomposed when QiliSDK has no
    class for it)."""
    from qilisdk.digital import (
        CNOT,
        CZ,
        RX,
        RY,
        RZ,
        SWAP,
        U1,
        U3,
        Adjoint,
        Controlled,
        H,
        I,
        S,
        T,
        X,
        Y,
        Z,
    )

    name, params = canonical(name, params)
    p = [float(v) for v in params]
    one = {"id": I, "x": X, "y": Y, "z": Z, "h": H, "s": S, "t": T}
    if name in one:
        circuit.add(one[name](qubits[0]))
    elif name == "sdg":
        circuit.add(Adjoint(S(qubits[0])))
    elif name == "tdg":
        circuit.add(Adjoint(T(qubits[0])))
    elif name == "sx":  # sqrt(X) = RX(pi/2) up to a global phase
        circuit.add(RX(qubits[0], theta=math.pi / 2))
    elif name == "rx":
        circuit.add(RX(qubits[0], theta=p[0]))
    elif name == "ry":
        circuit.add(RY(qubits[0], theta=p[0]))
    elif name == "rz":
        circuit.add(RZ(qubits[0], phi=p[0]))
    elif name == "p":
        circuit.add(U1(qubits[0], phi=p[0]))
    elif name == "u":
        circuit.add(U3(qubits[0], theta=p[0], phi=p[1], gamma=p[2]))
    elif name == "cx":
        circuit.add(CNOT(*qubits))
    elif name == "cz":
        circuit.add(CZ(*qubits))
    elif name == "swap":
        circuit.add(SWAP(*qubits))
    elif name in ("cy", "ch", "cp", "crx", "cry", "crz", "cu"):
        target = {
            "cy": lambda q: Y(q),
            "ch": lambda q: H(q),
            "cp": lambda q: U1(q, phi=p[0]),
            "crx": lambda q: RX(q, theta=p[0]),
            "cry": lambda q: RY(q, theta=p[0]),
            "crz": lambda q: RZ(q, phi=p[0]),
            "cu": lambda q: U3(q, theta=p[0], phi=p[1], gamma=p[2]),
        }[name]
        circuit.add(Controlled(qubits[0], basic_gate=target(qubits[1])))
        if name == "cu" and p[3]:  # cu's phase: a phase gate on the control
            circuit.add(U1(qubits[0], phi=p[3]))
    elif name == "ccx":
        circuit.add(Controlled(qubits[0], qubits[1], basic_gate=X(qubits[2])))
    elif name == "rzz":  # exp(-i t/2 ZZ) = CX . RZ(t) on the target . CX
        a, b = qubits
        circuit.add(CNOT(a, b))
        circuit.add(RZ(b, phi=p[0]))
        circuit.add(CNOT(a, b))
    elif name in ("rxx", "ryy"):  # rotate into Z, rzz, rotate back
        a, b = qubits
        to_z = (H,) if name == "rxx" else None
        for q in (a, b):
            if to_z:
                circuit.add(H(q))
            else:
                circuit.add(RX(q, theta=math.pi / 2))
        _emit(circuit, "rzz", qubits, p)
        for q in (a, b):
            if to_z:
                circuit.add(H(q))
            else:
                circuit.add(RX(q, theta=-math.pi / 2))
    elif name == "iswap":  # S S . H(a) . CX(a,b) . CX(b,a) . H(b)
        a, b = qubits
        circuit.add(S(a))
        circuit.add(S(b))
        circuit.add(H(a))
        circuit.add(CNOT(a, b))
        circuit.add(CNOT(b, a))
        circuit.add(H(b))
    else:
        raise ValueError(f"QiliSDK: no conversion for the gate {name!r}")


def dag_to_circuit_qilisdk(dag, num_qubits=None):
    """``[circuit, observables]``: a QiliSDK ``Circuit`` on qubits ``0..n-1`` and the
    Z/I string its ``Observable I`` markers give (``None`` without markers).

    A mid-circuit ``measure`` (the measure-and-prepare gate cut) is not supported
    here: on a simulator Qdislib defers it to the end (see
    :func:`~Qdislib.core.cutting_algorithms.gate_cutting._defer_cut_measurements`).
    """
    import networkx as nx
    from qilisdk.digital import Circuit

    from Qdislib.utils.graph_qibo import _update_qubits_serie

    dag = dag.copy()
    identity = []
    for node in list(nx.topological_sort(dag)):
        if (dag.nodes[node].get("gate") or "").upper() == "OBSERVABLE I":
            identity.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)
    dag, width, _ = _update_qubits_serie(dag)
    circuit = Circuit(width)
    for node in nx.topological_sort(dag):
        data = dag.nodes[node]
        if (data.get("gate") or "").lower() in ("measure", "m"):
            raise NotImplementedError(
                "QiliSDK circuits from Qdislib carry no mid-circuit measurement; "
                "run the gate cut on a simulator (it is deferred) or another SDK."
            )
        _emit(
            circuit, data["gate"], tuple(data["qubits"]), data.get("parameters") or ()
        )
    if not identity:
        return [circuit, None]
    observables = ["Z"] * width
    for qubit in identity:
        observables[qubit] = "I"
    return [circuit, "".join(observables)]


def sample_qilisdk(circuit, shots, seed=None, gpu=False):
    """``{qubit-0-first bitstring: count}`` from QiliSim (every qubit read)."""
    from qilisdk.backends import ExecutionConfig, QiliSim
    from qilisdk.functionals import DigitalPropagation
    from qilisdk.readout import Readout

    config = ExecutionConfig(seed=None if seed is None else int(seed) % 2**31, gpu=gpu)
    result = QiliSim(execution_config=config).execute(
        DigitalPropagation(circuit), Readout().with_sampling(nshots=int(shots))
    )
    return {str(bits): int(count) for bits, count in result.sampling.samples.items()}


def qilisdk_num_qubits(circuit):
    """A QiliSDK circuit's width (idle qubits included)."""
    return int(circuit.nqubits)


def qilisdk_statevector(circuit):
    """The exact state from QiliSim, qubit 0 the most significant bit."""
    import numpy as np
    from qilisdk.backends import QiliSim
    from qilisdk.functionals import DigitalPropagation
    from qilisdk.readout import Readout

    state = (
        QiliSim()
        .execute(DigitalPropagation(circuit), Readout().with_state_tomography())
        .get_state()
    )
    dense = state.dense() if callable(getattr(state, "dense", None)) else state.dense
    return np.asarray(dense, dtype=complex).reshape(-1)
