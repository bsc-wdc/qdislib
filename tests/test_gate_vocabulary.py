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

"""One gate vocabulary (OpenQASM 3) for every DAG, and cut names in any SDK's."""

import importlib.util
import math
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None
QIBO = importlib.util.find_spec("qibo") is not None


class CanonicalTest(BaseTimedTestCase):
    def test_names_and_parameters(self):
        from Qdislib.utils.gates import canonical

        self.assertEqual(canonical("CNOT"), ("cx", ()))
        self.assertEqual(canonical("cx"), ("cx", ()))
        self.assertEqual(canonical("U1", (0.3,)), ("p", (0.3,)))
        self.assertEqual(canonical("CU1", (0.3,)), ("cp", (0.3,)))
        self.assertEqual(canonical("U3", (0.1, 0.2, 0.3)), ("u", (0.1, 0.2, 0.3)))
        self.assertEqual(canonical("U2", (0.2, 0.3)), ("u", (math.pi / 2, 0.2, 0.3)))
        # qibo's CU3 controls the SU(2) U3; OpenQASM's cu carries the phase
        self.assertEqual(
            canonical("CU3", (0.1, 0.2, 0.3)), ("cu", (0.1, 0.2, 0.3, -0.25))
        )
        self.assertEqual(
            canonical("cu3", (0.1, 0.2, 0.3), sdk="qiskit"),
            ("cu", (0.1, 0.2, 0.3, 0.0)),
        )
        self.assertEqual(canonical("Adjoint(S)"), ("sdg", ()))
        for marker in ("measure", "Observable I", "M", "barrier"):
            self.assertEqual(canonical(marker)[0], marker)

    @unittest.skipUnless(QIBO, "qibo not installed")
    def test_lowered_cu_is_the_same_unitary(self):
        """cu through qibo (CU3 + a phase on the control) is OpenQASM's cu."""
        from Qdislib.utils.graph_qibo import gate_matrix

        theta, phi, lam, gamma = 0.7, 0.3, -0.4, 0.9
        u = np.array(
            [
                [math.cos(theta / 2), -np.exp(1j * lam) * math.sin(theta / 2)],
                [
                    np.exp(1j * phi) * math.sin(theta / 2),
                    np.exp(1j * (phi + lam)) * math.cos(theta / 2),
                ],
            ]
        )
        expected = np.eye(4, dtype=complex)
        expected[2:, 2:] = np.exp(1j * gamma) * u  # control = qubit 0 = MSB in qibo
        got = gate_matrix("cu", (theta, phi, lam, gamma), 2)
        self.assertTrue(np.allclose(got, expected))


@unittest.skipUnless(QISKIT and QIBO, "qiskit and qibo needed")
class SameDagTest(BaseTimedTestCase):
    def test_qiskit_and_qibo_give_the_same_gates(self):
        from qibo import gates, models
        from qiskit import QuantumCircuit

        import Qdislib as qd
        from Qdislib.utils.backends import convert_to_dag

        qk = QuantumCircuit(3)
        qk.h(0)
        qk.cx(0, 1)
        qk.p(0.3, 2)
        qk.cp(0.4, 1, 2)
        qk.u(0.1, 0.2, 0.3, 0)
        qb = models.Circuit(3)
        qb.add(gates.H(0))
        qb.add(gates.CNOT(0, 1))
        qb.add(gates.U1(2, 0.3))
        qb.add(gates.CU1(1, 2, 0.4))
        qb.add(gates.U3(0, 0.1, 0.2, 0.3))

        def gates_of(circuit):
            dag = convert_to_dag(circuit)
            return [
                (dag.nodes[n]["gate"], tuple(dag.nodes[n]["qubits"]))
                for n in sorted(dag.nodes, key=lambda n: int(n.rsplit("_", 1)[1]))
            ]

        self.assertEqual(gates_of(qk), gates_of(qb))
        self.assertEqual([g for g, _ in gates_of(qk)], ["h", "cx", "p", "cp", "u"])
        self.assertAlmostEqual(
            qd.gate_cutting(qk, ["CX_2"], shots=None),
            qd.gate_cutting(qb, ["CNOT_2"], shots=None),
            places=10,
        )

    def test_edges_match_the_old_rule(self):
        """blue from the last gate on each qubit, red from the last one not linked."""
        from qiskit import QuantumCircuit

        from Qdislib.utils.backends import convert_to_dag

        qc = QuantumCircuit(2)
        qc.h(0)  # H_1
        qc.cx(0, 1)  # CX_2
        qc.x(1)  # X_3
        qc.cx(0, 1)  # CX_4
        dag = convert_to_dag(qc)
        color = {(a, b): d["color"] for a, b, d in dag.edges(data=True)}
        self.assertEqual(color[("H_1", "CX_2")], "blue")
        self.assertEqual(color[("CX_2", "X_3")], "blue")
        self.assertEqual(color[("CX_2", "CX_4")], "blue")  # qubit 0's wire
        self.assertEqual(color[("X_3", "CX_4")], "blue")
        self.assertEqual(color[("H_1", "CX_4")], "red")  # qubit 0, behind CX_2


@unittest.skipUnless(QISKIT, "qiskit not installed")
class CutNamesTest(BaseTimedTestCase):
    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        for q in range(4):
            qc.ry(0.3 + q, q)
        qc.cx(0, 1)  # CX_5
        qc.cz(2, 3)  # CZ_6
        qc.cx(1, 2)  # CX_7
        for q in range(4):
            qc.ry(0.1 * q, q)
        return qc

    def test_any_sdk_name_for_the_gate(self):
        import Qdislib as qd
        from Qdislib.utils.gates import resolve_cuts

        qc = self._circuit()
        self.assertEqual(resolve_cuts(qc, ["CNOT_7"]), ["CX_7"])
        self.assertEqual(resolve_cuts(qc, ["cx_7"]), ["CX_7"])
        self.assertEqual(
            resolve_cuts(qc, {"gate": ["cnot_7"], "wire": [("RY_2", "CNOT_5")]}),
            {"gate": ["CX_7"], "wire": [("RY_2", "CX_5")]},
        )
        exact = qd.gate_cutting(qc, ["CX_7"], shots=None)
        self.assertAlmostEqual(qd.gate_cutting(qc, ["CNOT_7"], shots=None), exact, 10)
        self.assertAlmostEqual(
            qd.circuit_cutting(qc, {"gate": ["cnot_7"]}, shots=None), exact, 10
        )

    @unittest.skipUnless(importlib.util.find_spec("qibo"), "qibo not installed")
    def test_a_cu3_can_be_named(self):
        # its parameters are not needed to compare names
        from qibo import gates, models

        from Qdislib.utils.gates import resolve_cuts

        qc = models.Circuit(2)
        qc.add(gates.H(0))
        qc.add(gates.CU3(0, 1, 0.1, 0.2, 0.3))
        self.assertEqual(resolve_cuts(qc, ["cu_2"]), ["CU3_2"])
        self.assertEqual(resolve_cuts(qc, ["CU3_2"]), ["CU3_2"])
        with self.assertRaisesRegex(ValueError, "gate 2 of this circuit"):
            resolve_cuts(qc, ["cz_2"])

    def test_wrong_gate_is_named(self):
        import Qdislib as qd

        with self.assertRaisesRegex(ValueError, "gate 7 of this circuit is 'CX_7'"):
            qd.gate_cutting(self._circuit(), ["CZ_7"], shots=None)


class GateNamesTest(BaseTimedTestCase):
    def test_iter_gates_names_are_the_dag_labels(self):
        """What find_cut and draw_cut call a gate is a node of the DAG."""
        from Qdislib.utils.backends import convert_to_dag, iter_gates

        circuits = []
        if QISKIT:
            from qiskit import QuantumCircuit

            qc = QuantumCircuit(2)
            qc.h(0)
            qc.cx(0, 1)
            qc.p(0.2, 1)
            circuits.append(qc)
        if QIBO:
            from qibo import gates, models

            qb = models.Circuit(2)
            qb.add(gates.H(0))
            qb.add(gates.CNOT(0, 1))
            circuits.append(qb)
        if importlib.util.find_spec("pennylane") is not None:
            import pennylane as qml

            with qml.tape.QuantumTape() as tape:
                qml.Hadamard(0)
                qml.CNOT([0, 1])
                qml.RY(0.3, 1)
            circuits.append(tape)
        for circuit in circuits:
            names = {
                f"{name}_{i}"
                for i, (name, qubits) in enumerate(iter_gates(circuit), 1)
                if qubits
            }
            self.assertEqual(names, set(convert_to_dag(circuit).nodes), circuit)


@unittest.skipUnless(QISKIT, "qiskit not installed")
class StatevectorTest(BaseTimedTestCase):
    def test_every_canonical_gate_matches_qiskit(self):
        """Exact mode's own simulator agrees with Qiskit on every canonical gate,
        each followed by more gates (so a wrong phase would show)."""
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Statevector

        from Qdislib.utils.backends import convert_to_dag
        from Qdislib.utils.statevector import dag_probabilities

        rng = np.random.default_rng(0)
        a = lambda: float(rng.uniform(0.2, 2.8))  # noqa: E731
        qc = QuantumCircuit(4)

        def mix():
            for q in range(4):
                qc.ry(a(), q)
                qc.rz(a(), q)

        mix()
        for add in (
            lambda: (qc.h(0), qc.x(1), qc.y(2), qc.z(3)),
            lambda: (qc.s(0), qc.sdg(1), qc.t(2), qc.tdg(3)),
            lambda: (qc.sx(0), qc.id(1), qc.rx(a(), 2), qc.p(a(), 3)),
            lambda: (qc.u(a(), a(), a(), 0), qc.cx(0, 1), qc.cy(1, 2), qc.cz(2, 3)),
            lambda: (qc.ch(3, 0), qc.cp(a(), 0, 2), qc.crx(a(), 1, 3)),
            lambda: (qc.cry(a(), 2, 0), qc.crz(a(), 3, 1)),
            lambda: (qc.cu(a(), a(), a(), a(), 0, 3), qc.swap(1, 2)),
            lambda: (qc.iswap(0, 1), qc.rxx(a(), 2, 3), qc.ryy(a(), 0, 2)),
            lambda: (qc.rzz(a(), 1, 3), qc.ccx(0, 1, 2), qc.cswap(3, 0, 1)),
        ):
            add()
            mix()
        want = Statevector(qc).probabilities()
        got = dag_probabilities(convert_to_dag(qc), 4)
        self.assertTrue(np.allclose(want, got, atol=1e-12))


if __name__ == "__main__":
    unittest.main()
