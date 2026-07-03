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

"""PennyLane backend integration: conversion, find_cut, gate/wire cutting."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

PENNYLANE_AVAILABLE = importlib.util.find_spec("pennylane") is not None
QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None

# A small circuit with a clearly non-zero <ZZZZ> (small rotations keep it near |0>),
# so a matching reconstruction is meaningful rather than trivially ~0. The Qibo and
# PennyLane builders below apply the SAME gates, so Qibo's analytical value is the
# reference for both.
_GATES = [
    ("RY", 0, np.pi / 6),
    ("CZ", (0, 1), None),
    ("RY", 1, np.pi / 7),
    ("CZ", (1, 2), None),
    ("RY", 2, np.pi / 6),
    ("CZ", (2, 3), None),
    ("RY", 3, np.pi / 8),
]
_NQUBITS = 4


def _qibo_circuit():
    from qibo import models, gates

    circuit = models.Circuit(_NQUBITS)
    for name, target, param in _GATES:
        if name == "RY":
            circuit.add(gates.RY(target, param))
        elif name == "CZ":
            circuit.add(gates.CZ(*target))
    return circuit


def _pennylane_tape():
    import pennylane as qml

    with qml.tape.QuantumTape() as tape:
        for name, target, param in _GATES:
            if name == "RY":
                qml.RY(param, wires=target)
            elif name == "CZ":
                qml.CZ(wires=list(target))
    return tape


@unittest.skipUnless(PENNYLANE_AVAILABLE, "pennylane not installed")
class PennylaneConversionTest(BaseTimedTestCase):
    def test_circuit_to_dag_roundtrip(self):
        """A PennyLane tape converts to a DAG and back to the same gate sequence."""
        from Qdislib.utils.graph_pennylane import (
            circuit_pennylane_to_dag,
            dag_to_circuit_pennylane,
            PennylaneCircuit,
        )

        dag = circuit_pennylane_to_dag(_pennylane_tape())
        # Canonical (normalised) gate names in the DAG.
        gate_names = sorted(d["gate"] for _, d in dag.nodes(data=True))
        self.assertEqual(gate_names, ["CZ", "CZ", "CZ", "RY", "RY", "RY", "RY"])

        built = dag_to_circuit_pennylane(dag, _NQUBITS)
        circuit, observables = built
        self.assertIsInstance(circuit, PennylaneCircuit)
        self.assertEqual(circuit.num_qubits, _NQUBITS)
        self.assertEqual(observables, "Z" * _NQUBITS)
        # Same multiset of operations survives the round-trip.
        self.assertEqual(sorted(op[0] for op in circuit.ops), gate_names)

    def test_iter_gates_and_num_qubits(self):
        """find_cut's backend-neutral helpers understand a PennyLane circuit."""
        from Qdislib.utils.backends import iter_gates, num_qubits

        tape = _pennylane_tape()
        self.assertEqual(num_qubits(tape), _NQUBITS)
        names = [name for name, _ in iter_gates(tape)]
        self.assertEqual(names.count("RY"), 4)
        self.assertEqual(names.count("CZ"), 3)

    def test_from_qasm_returns_pennylane_tape(self):
        import pennylane as qml
        import Qdislib as qd

        qasm = 'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n'
        tape = qd.from_qasm(qasm, software="pennylane")
        self.assertIsInstance(tape, qml.tape.QuantumScript)
        names = [op.name for op in tape.operations]
        self.assertEqual(names, ["Hadamard", "CNOT"])

    @unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
    def test_semantic_key_matches_qibo(self):
        """The same circuit hashes identically for PennyLane and Qibo (cache dedup)."""
        from Qdislib.utils.circuit_cache import semantic_key
        from Qdislib.utils.graph_pennylane import PennylaneCircuit
        from qibo import models, gates

        pennylane_circuit = PennylaneCircuit(
            [("H", (0,), []), ("CNOT", (0, 1), [])], 2
        )
        qibo_circuit = models.Circuit(2)
        qibo_circuit.add(gates.H(0))
        qibo_circuit.add(gates.CNOT(0, 1))

        key_pennylane = semantic_key(pennylane_circuit, "pennylane")
        key_qibo = semantic_key(qibo_circuit, "qibo")
        if key_pennylane is not None and key_qibo is not None:  # skip if PyZX absent
            self.assertEqual(key_pennylane, key_qibo)


@unittest.skipUnless(
    PENNYLANE_AVAILABLE and QIBO_AVAILABLE, "pennylane and qibo required"
)
class PennylaneCuttingTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        from Qdislib.utils.circuit import analytical_solution

        self.analytical = analytical_solution(_qibo_circuit(), "Z" * _NQUBITS)

    def test_gate_cutting(self):
        import Qdislib as qd

        circuit = _qibo_circuit()
        cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)
        reconstruction = qd.gate_cutting(circuit, cut, software="pennylane", shots=20000)
        self.assertLess(abs(reconstruction - self.analytical), 0.2)

    def test_wire_cutting(self):
        import Qdislib as qd

        circuit = _qibo_circuit()
        cut = qd.find_cut(circuit, gate_cut=False)
        reconstruction = qd.wire_cutting(circuit, cut, software="pennylane", shots=20000)
        self.assertLess(abs(reconstruction - self.analytical), 0.2)

    def test_gate_cutting_three_step(self):
        import Qdislib as qd

        circuit = _qibo_circuit()
        cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)
        subcircuits = qd.gate_cutting_subcircuits(circuit, cut, software="pennylane")
        values = [qd.expectation_value(s, shots=20000) for s in subcircuits]
        reconstruction = subcircuits.reconstruct(values)
        self.assertLess(abs(reconstruction - self.analytical), 0.2)

    def test_gate_cutting_pennylane_input_circuit(self):
        """A PennyLane circuit as the cut input (exercises circuit_pennylane_to_dag)."""
        import Qdislib as qd

        tape = _pennylane_tape()
        cut = qd.find_cut(tape, gate_cut=True, wire_cut=False)
        self.assertTrue(cut)  # find_cut works on a PennyLane circuit
        reconstruction = qd.gate_cutting(tape, cut, software="pennylane", shots=20000)
        self.assertLess(abs(reconstruction - self.analytical), 0.2)


if __name__ == "__main__":
    unittest.main()
