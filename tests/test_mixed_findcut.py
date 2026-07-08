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
        # Qubit 1 shuttles between 0 and 2, so a wire cut on its timeline plus a gate
        # cut elsewhere beats an all-gate or all-wire cut -> a genuinely mixed result.
        from qibo import gates, models

        c = models.Circuit(4)
        for q in range(4):
            c.add(gates.RY(q, 0.4 + 0.15 * q))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.RY(1, 0.6))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.CZ(2, 3))
        return c

    def test_returns_gate_wire_dict(self):
        import Qdislib as qd

        mixed = qd.find_cut(self._messenger(), mixed=True, seed=0)
        self.assertEqual(set(mixed), {"gate", "wire"})
        self.assertIsInstance(mixed["gate"], list)
        self.assertIsInstance(mixed["wire"], list)
        # gate cuts are gate ids; wire cuts are (source, target) pairs.
        self.assertTrue(all(isinstance(g, str) for g in mixed["gate"]))
        self.assertTrue(all(len(w) == 2 for w in mixed["wire"]))

    def test_discovers_a_genuinely_mixed_cut(self):
        import Qdislib as qd

        mixed = qd.find_cut(self._messenger(), mixed=True, seed=0)
        self.assertTrue(mixed["gate"], "expected at least one gate cut")
        self.assertTrue(mixed["wire"], "expected at least one wire cut")

    def test_respects_forbidden(self):
        import Qdislib as qd

        mixed = qd.find_cut(
            self._messenger(), mixed=True, seed=0, forbidden_gates=["CZ_6", "CZ_8"]
        )
        self.assertNotIn("CZ_6", mixed["gate"])
        self.assertNotIn("CZ_8", mixed["gate"])

    def test_reproducible_with_seed(self):
        import Qdislib as qd

        a = qd.find_cut(self._messenger(), mixed=True, seed=3)
        b = qd.find_cut(self._messenger(), mixed=True, seed=3)
        self.assertEqual(a, b)


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
