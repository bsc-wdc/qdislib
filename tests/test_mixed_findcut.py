#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""find_cut(mixed=True): discover a combined gate + wire cut for gate_wire_cutting."""

import importlib.util
import unittest

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None


@unittest.skipUnless(QIBO, "qibo not installed")
class MixedFindCutTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _messenger(self):
        # Blocks A = {0, 1} and B = {2, 3}. Qubit 4 talks to A, then to B, so cutting
        # its wire once (gamma 4) beats cutting its gates; one direct CZ links A and
        # B (a gate cut, gamma 3). Best: one wire cut + one gate cut, gamma 12.
        from qibo import gates, models

        c = models.Circuit(5)
        for q in range(5):
            c.add(gates.RY(q, 0.3 + 0.1 * q))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(2, 3))
        for q in (0, 1):
            c.add(gates.CZ(4, q))
        c.add(gates.RY(4, 0.7))  # RY_10
        for q in (2, 3):
            c.add(gates.CZ(4, q))
        c.add(gates.CZ(1, 2))  # CZ_13
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(2, 3))
        return c

    def test_returns_the_circuit_cutting_dict(self):
        import Qdislib as qd

        mixed = qd.find_cut(self._messenger(), max_qubits=3, mixed=True, seed=0)
        self.assertEqual(set(mixed), {"gate", "hadamard", "wire"})
        self.assertTrue(all(isinstance(g, str) for g in mixed["gate"]))
        self.assertTrue(all(len(w) == 2 for w in mixed["wire"]))

    def test_discovers_a_genuinely_mixed_cut(self):
        import Qdislib as qd

        circuit = self._messenger()
        mixed = qd.find_cut(circuit, max_qubits=3, mixed=True, seed=0)
        self.assertEqual(mixed["gate"], ["CZ_13"])
        self.assertEqual(mixed["wire"], [("CZ_9", "RY_10")])
        self.assertEqual(qd.cutting_overhead(circuit, mixed)["gamma"], 12)
        value = qd.circuit_cutting(circuit, mixed, shots=None)
        self.assertAlmostEqual(
            value, qd.analytical_solution(circuit, "ZZZZZ"), places=10
        )

    def test_respects_forbidden(self):
        import Qdislib as qd

        mixed = qd.find_cut(
            self._messenger(), mixed=True, seed=0, forbidden_gates=["CZ_13"]
        )
        self.assertNotIn("CZ_13", mixed["gate"] + mixed["hadamard"])

    def test_reproducible_with_seed(self):
        import Qdislib as qd

        a = qd.find_cut(self._messenger(), mixed=True, seed=3)
        b = qd.find_cut(self._messenger(), mixed=True, seed=3)
        self.assertEqual(a, b)


QISKIT = importlib.util.find_spec("qiskit") is not None


def _qaoa(n, seed=0):
    import networkx as nx
    import numpy as np
    from qiskit import QuantumCircuit

    graph = nx.random_regular_graph(3, n, seed=seed)
    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    for q in range(n):
        qc.h(q)
    for a, b in graph.edges():
        qc.rzz(float(rng.uniform(0.1, 1.0)), a, b)
    for q in range(n):
        qc.rx(float(rng.uniform(0, 3)), q)
    return qc


def _brickwork(n, depth, seed=0):
    import numpy as np
    from qiskit import QuantumCircuit

    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    for d in range(depth):
        for q in range(n):
            qc.ry(float(rng.uniform(0, 3)), q)
        for q in range(d % 2, n - 1, 2):
            qc.cz(q, q + 1)
    return qc


@unittest.skipUnless(QISKIT, "qiskit not installed")
class CostAwareCutTest(BaseTimedTestCase):
    """find_cut(mixed=True): what it returns always runs, and it is priced right."""

    def _check(self, circuit, cut, max_qubits):
        """Every cut separates (the plan raises otherwise) and every fragment, its
        Hadamard ancillas included, fits ``max_qubits``."""
        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan

        plan = _CutPlan(circuit, cut, None, False, False, True)
        for fragment in plan.fragments:
            qubits = {
                q
                for node in fragment["nodes"]
                for q in plan.graph.nodes[node].get("qubits") or ()
            }
            self.assertLessEqual(len(qubits) + len(fragment["ancilla_of"]), max_qubits)

    def test_valid_within_limits_and_never_worse_than_gate_only(self):
        import Qdislib as qd

        for circuit, max_qubits in (
            (_brickwork(8, 4), 4),
            (_brickwork(12, 6), 6),
            (_qaoa(10), 5),
            (_qaoa(12), 6),
        ):
            with self.subTest(n=circuit.num_qubits, max_qubits=max_qubits):
                mixed = qd.find_cut(circuit, max_qubits=max_qubits, mixed=True, seed=0)
                self._check(circuit, mixed, max_qubits)
                gate_only = qd.find_cut(
                    circuit,
                    max_qubits=max_qubits,
                    mixed=True,
                    seed=0,
                    hadamard=False,
                    wire_cut=False,
                )
                self._check(circuit, gate_only, max_qubits)
                self.assertLessEqual(
                    qd.cutting_overhead(circuit, mixed)["gamma"],
                    qd.cutting_overhead(circuit, gate_only)["gamma"] + 1e-9,
                )

    def test_angle_aware_hadamard_cuts_on_qaoa(self):
        import Qdislib as qd

        circuit = _qaoa(10)
        mixed = qd.find_cut(circuit, max_qubits=5, mixed=True, seed=0)
        self.assertTrue(mixed["gate"] or mixed["hadamard"])
        # Both cuts follow an RZZ's angle (1 + 2|sin phi|); measure-and-prepare
        # used to cut it as 2 CZs, gamma 9 each.
        self.assertLess(qd.cutting_overhead(circuit, mixed)["gamma"], 729 / 10)

    def test_time_limit(self):
        import time

        import Qdislib as qd

        circuit = _brickwork(60, 8)
        start = time.time()
        cut = qd.find_cut(circuit, max_qubits=30, mixed=True, seed=0, time_limit=3)
        self.assertLess(time.time() - start, 3 + 4)  # set-up and final check on top
        self._check(circuit, cut, 30)


@unittest.skipUnless(QIBO, "qibo not installed")
class MixedCutExecutionTest(BaseTimedTestCase):
    """A discovered mixed cut executes and reconstructs, including the *adjacent* case
    (a wire cut sitting right after a gate cut, resolved by the wire-then-gate surgery
    order in gate_wire_cutting)."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _circuit(self):
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.6))
        c.add(gates.RY(1, 0.4))
        c.add(gates.RY(2, 0.5))
        c.add(gates.CZ(0, 1))  # CZ_4
        c.add(gates.CZ(1, 2))  # CZ_5  -> gate cut
        c.add(gates.RY(2, 0.3))  # RY_6
        c.add(gates.RX(2, 0.2))  # RX_7
        return c

    def test_adjacent_mixed_cut_reconstructs(self):
        # Wire cut (CZ_5, RY_6) sits right after the gate cut CZ_5 on qubit 2.
        import Qdislib as qd

        circuit = self._circuit()
        analytical = qd.analytical_solution(circuit, "ZZZ")
        got = qd.gate_wire_cutting(
            self._circuit(),
            gates_cut=["CZ_5"],
            wire_cut=[("CZ_5", "RY_6")],
            observable="ZZZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.2)

    def test_disjoint_mixed_cut_reconstructs(self):
        import Qdislib as qd

        circuit = self._circuit()
        analytical = qd.analytical_solution(circuit, "ZZZ")
        got = qd.gate_wire_cutting(
            self._circuit(),
            gates_cut=["CZ_4"],
            wire_cut=[("CZ_5", "RY_6")],
            observable="ZZZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.2)


if __name__ == "__main__":
    unittest.main()
