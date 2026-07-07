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

"""
Qdislib PennyLane graph utils.

Converts a PennyLane circuit (a :class:`~pennylane.tape.QuantumScript` /
``QuantumTape``, or a :class:`~pennylane.QNode`) to and from Qdislib's
backend-neutral ``networkx`` DAG, mirroring :mod:`~Qdislib.utils.graph_qibo` and
:mod:`~Qdislib.utils.graph_cudaq`.

PennyLane names its operations differently from the other SDKs (``Hadamard``
rather than ``H``, ``PauliX`` rather than ``X``), so :func:`circuit_pennylane_to_dag`
normalises them to the same canonical vocabulary the other converters use
(``H``/``X``/``Y``/``Z``/``S``/``T``/``SDG``/``TDG``/``RX``/``RY``/``RZ``/``CNOT``/``CZ``);
:func:`dag_to_circuit_pennylane` accepts that whole vocabulary back (including the
``CX`` alias Qiskit/CUDA-Q emit and the gate-cut ``MEASURE`` node), so a cut found
on a Qibo/Qiskit/CUDA-Q circuit can still be executed with ``software="pennylane"``.
"""

from __future__ import annotations  # defer annotation eval so pennylane can be optional

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_pennylane_to_dag",
    "dag_to_circuit_pennylane",
    "PennylaneCircuit",
]

import networkx as nx

try:
    import pennylane as qml

    pennylane_available = True
except ImportError:  # pennylane is optional; functions that need it fail when called
    qml = None
    pennylane_available = False

from Qdislib.utils.runtime import task
from Qdislib.utils.graph_qibo import _update_qubits_serie

#: PennyLane operation name -> Qdislib's canonical DAG gate name. Kept in step with
#: :func:`_apply_pennylane_gate` (the inverse) and with the other ``circuit_*_to_dag``
#: converters so the ``f"{name}_{index}"`` node identifiers line up across SDKs.
_PENNYLANE_TO_CANONICAL = {
    "Hadamard": "H",
    "PauliX": "X",
    "PauliY": "Y",
    "PauliZ": "Z",
    "S": "S",
    "T": "T",
    "Adjoint(S)": "SDG",
    "Adjoint(T)": "TDG",
    "SX": "SX",
    "CNOT": "CNOT",
    "CZ": "CZ",
    "CY": "CY",
    "RX": "RX",
    "RY": "RY",
    "RZ": "RZ",
}


class PennylaneCircuit:
    """A PennyLane subcircuit as an ordered, picklable operation list.

    The native PennyLane execution objects -- a ``QNode`` (bound to a device and a
    Python quantum function) or a ``QuantumScript`` carrying mid-circuit
    ``MeasurementValue`` handles -- do not survive a PyCOMPSs task boundary, so, as
    with the CUDA-Q ``(dag, num_qubits)`` payload, Qdislib carries the plain data
    (``ops`` + ``num_qubits``) and rebuilds a QNode **in the process that executes
    it** (see :func:`~Qdislib.core.cutting_algorithms.wire_cutting._counts_pennylane`).

    :param ops: Ordered ``(canonical_gate_name, qubits_tuple, parameters_list)``
        triples; the gate-cut collapse point is the pseudo-gate ``"MEASURE"``.
    :param num_qubits: Width of the (renumbered 0..n-1) subcircuit.
    """

    def __init__(self, ops, num_qubits):
        self.ops = ops
        self.num_qubits = num_qubits

    def __repr__(self):
        return f"PennylaneCircuit(num_qubits={self.num_qubits}, gates={len(self.ops)})"


def _apply_pennylane_gate(gate, qubits, parameters):
    """Apply one canonical DAG gate inside a PennyLane QNode / recording context.

    The inverse of :data:`_PENNYLANE_TO_CANONICAL`, extended with the ``CX`` alias
    (Qiskit/CUDA-Q) and the gate-cut ``MEASURE`` node.

    :param gate: Canonical (upper-cased) gate name.
    :param qubits: Tuple of integer wire indices.
    :param parameters: List of gate parameters (rotation angles).
    :raises ValueError: For a gate outside the supported set.
    :return: The ``MeasurementValue`` for a ``MEASURE`` (so the caller can sample
        the collapse outcome), otherwise ``None``.
    """
    g = gate.upper()
    if g == "H":
        qml.Hadamard(qubits[0])
    elif g == "X":
        qml.PauliX(qubits[0])
    elif g == "Y":
        qml.PauliY(qubits[0])
    elif g == "Z":
        qml.PauliZ(qubits[0])
    elif g == "S":
        qml.S(qubits[0])
    elif g == "T":
        qml.T(qubits[0])
    elif g == "SDG":
        qml.adjoint(qml.S)(qubits[0])
    elif g == "TDG":
        qml.adjoint(qml.T)(qubits[0])
    elif g == "SX":
        qml.SX(qubits[0])
    elif g in ("RX", "RY", "RZ"):
        angle = float(parameters[0]) if parameters else 0.0
        getattr(qml, g)(angle, wires=qubits[0])
    elif g == "U1":  # phase gate diag(1, e^{i lam})
        qml.PhaseShift(float(parameters[0]), wires=qubits[0])
    elif g == "U2":
        qml.U2(float(parameters[0]), float(parameters[1]), wires=qubits[0])
    elif g == "U3":
        qml.U3(*[float(p) for p in parameters], wires=qubits[0])
    elif g in ("CNOT", "CX"):
        qml.CNOT(wires=[qubits[0], qubits[1]])
    elif g == "CZ":
        qml.CZ(wires=[qubits[0], qubits[1]])
    elif g == "CY":
        qml.CY(wires=[qubits[0], qubits[1]])
    elif g in ("CU1", "CPHASE"):  # controlled phase on |11>
        qml.ControlledPhaseShift(float(parameters[0]), wires=[qubits[0], qubits[1]])
    elif g == "CU3":
        # Controlled SU(2) U3. qml.U3 carries an e^{i(phi+lam)/2} global phase, which
        # under control becomes a relative phase on the control-1 block; cancel it with
        # a phase on the control qubit so the block is exactly qibo's SU(2) U3.
        theta, phi, lam = (float(p) for p in parameters)
        qml.ctrl(qml.U3(theta, phi, lam, wires=qubits[1]), control=qubits[0])
        qml.PhaseShift(-(phi + lam) / 2, wires=qubits[0])
    elif g == "SWAP":
        qml.SWAP(wires=[qubits[0], qubits[1]])
    elif g == "ISWAP":
        qml.ISWAP(wires=[qubits[0], qubits[1]])
    elif g == "RZZ":  # qibo RZZ(theta) == qml.IsingZZ(theta)
        qml.IsingZZ(float(parameters[0]), wires=[qubits[0], qubits[1]])
    elif g in ("TOFFOLI", "CCX"):
        qml.Toffoli(wires=[qubits[0], qubits[1], qubits[2]])
    elif g == "MEASURE":
        # Gate-cut measure: a projective Z measurement that collapses the qubit to
        # its outcome |s> and lets it CONTINUE (the CZ cut's (I +/- Z)/2 projector,
        # Mitarai-Fujii arXiv:1909.07534). qml.measure does exactly this; the
        # returned MeasurementValue is sampled for the reconstruction sign.
        return qml.measure(qubits[0])
    elif g in ("M", "OBSERVABLE I"):
        # A plain terminal measurement / identity-observable marker: not a unitary,
        # handled by the caller (which measures every qubit itself).
        return None
    else:
        raise ValueError(f"Unsupported gate type: {gate}")
    return None


def _build_pennylane_qnode(qfunc, num_qubits, shots, device="default.qubit", seed=None):
    """Build a shot-based QNode, robust across PennyLane releases.

    Uses the ``qml.set_shots`` transform on PennyLane >= 0.40; on older releases it
    sets the shots on the device instead (the transform did not exist yet, and
    setting shots on the device is not deprecated there). A ``seed`` is forwarded
    to the device when given (supported by ``default.qubit``).

    :param qfunc: A PennyLane quantum function (applies the gates, returns the
        measurement(s)).
    :param num_qubits: Number of device wires.
    :param shots: Number of measurement shots.
    :param device: PennyLane device short name.
    :param seed: Optional RNG seed for reproducible sampling.
    :return: A callable ``QNode`` bound to ``shots`` shots.
    """
    device_kwargs = {"wires": num_qubits}
    if seed is not None:
        device_kwargs["seed"] = seed
    if hasattr(qml, "set_shots"):
        dev = qml.device(device, **device_kwargs)
        return qml.set_shots(qml.QNode(qfunc, dev), shots=shots)
    device_kwargs["shots"] = shots
    dev = qml.device(device, **device_kwargs)
    return qml.QNode(qfunc, dev)


def _to_tape(circuit):
    """Return the ``QuantumScript`` for a PennyLane circuit or QNode.

    :param circuit: A ``QuantumScript``/``QuantumTape`` (returned as-is) or a
        parameterless ``QNode`` (traced into its tape).
    :raises TypeError: If ``circuit`` is neither.
    """
    if isinstance(circuit, qml.tape.QuantumScript):
        return circuit
    if isinstance(circuit, qml.QNode):
        return qml.workflow.construct_tape(circuit)()
    raise TypeError(
        f"Unsupported PennyLane circuit type: {type(circuit).__name__}; "
        "expected a QuantumScript/QuantumTape or a QNode."
    )


def circuit_pennylane_to_dag(circuit, obs_I=None) -> nx.DiGraph:
    """Convert a PennyLane circuit into a Qdislib DAG.

    :param circuit: A PennyLane ``QuantumScript``/``QuantumTape`` or ``QNode``.
    :param obs_I: Optional observable string; each ``"I"`` appends an identity
        observable sink node on that qubit (internal use when cutting).
    :return: A :class:`networkx.DiGraph` whose nodes carry ``gate`` (canonical
        upper-cased name), ``qubits`` (tuple of ints) and ``parameters`` (list of
        floats), with blue (nearest-predecessor) and red (other-predecessor) edges.
    """
    tape = _to_tape(circuit)
    dag = nx.DiGraph()

    for gate_idx, op in enumerate(tape.operations, start=1):
        name = _PENNYLANE_TO_CANONICAL.get(op.name)
        if name is None:
            raise ValueError(
                f"Unsupported PennyLane operation '{op.name}'; supported: "
                f"{sorted(_PENNYLANE_TO_CANONICAL)}."
            )
        qubits = tuple(int(w) for w in op.wires)
        parameters = [float(p) for p in op.parameters]
        gate_name = f"{name}_{gate_idx}"
        dag.add_node(gate_name, gate=name, qubits=qubits, parameters=parameters)

        # Blue edges: nearest preceding gate on each shared qubit.
        for qubit in qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                ):
                    dag.add_edge(pred_gate, gate_name, color="blue")
                    break

        # Red edges: other earlier gates on a shared qubit, if not already linked.
        for qubit in qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                    and not dag.has_edge(pred_gate, gate_name)
                ):
                    dag.add_edge(pred_gate, gate_name, color="red")

    # Append identity-observable sink nodes when requested, after the LAST gate
    # acting on each identity qubit (see circuit_qibo_to_dag for the rationale).
    if obs_I:
        for idx, observable in enumerate(obs_I):
            if observable != "I":
                continue
            last = None
            for node, data in dag.nodes(data=True):
                if (
                    idx in (data.get("qubits") or ())
                    and data.get("gate") != "Observable I"
                ):
                    last = node
            if last is not None:
                obs_node = f"OBSI_{idx}"
                dag.add_node(
                    obs_node, gate="Observable I", qubits=(idx,), parameters=[]
                )
                dag.add_edge(last, obs_node, color="blue")

    return dag


@task(returns=1)
def dag_to_circuit_pennylane(dag: nx.DiGraph, num_qubits: int):
    """Rebuild a PennyLane subcircuit from a DAG.

    Mirrors :func:`~Qdislib.utils.graph_qibo.dag_to_circuit_qibo`: strips the
    ``Observable I`` nodes into a per-qubit observable string, compacts the qubit
    indices to a contiguous ``0..n-1`` range, and returns the ordered operation
    list wrapped in a :class:`PennylaneCircuit` (a QNode is built from it, on the
    device, only at execution time).

    :param dag: A connected-component DAG of the cut circuit.
    :param num_qubits: Qubit count of the original circuit (unused; the compacted
        width is recomputed, matching the other ``dag_to_circuit_*`` helpers).
    :return: ``[PennylaneCircuit, observables]``, or ``None`` for an empty DAG.
    """
    if dag is None:  # empty component (e.g. a qubit touched only by the cut gate)
        return None

    # Extract identity-observable markers before compaction (their positions become
    # "I" in the returned observable string).
    obs_I = []
    for node in list(nx.topological_sort(dag)):
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)

    dag, highest_qubit, _smallest = _update_qubits_serie(dag)

    ops = []
    for node in nx.topological_sort(dag):
        node_data = dag.nodes[node]
        gate = node_data["gate"].upper()
        if gate == "OBSERVABLE I":
            continue
        ops.append(
            (gate, tuple(node_data["qubits"]), list(node_data["parameters"] or []))
        )

    observables = ["Z"] * highest_qubit
    for qubit in obs_I:
        observables[qubit] = "I"
    observables = "".join(observables)

    return [PennylaneCircuit(ops, highest_qubit), observables]
