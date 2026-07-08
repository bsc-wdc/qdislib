#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Hardware-aware "SparseCut" gate-cut selection in find_cut (arXiv:2511.05492)."""

import importlib.util
import unittest

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None


@unittest.skipUnless(QIBO, "qibo not installed")
class SparseCutTest(BaseTimedTestCase):
    """``find_cut(coupling_map=...)`` cuts the longest-routing two-qubit gates."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _circuit(self):
        from qibo import gates, models

        c = models.Circuit(5)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))  # CZ_2  adjacent on a line (distance 1)
        c.add(gates.CZ(1, 2))  # CZ_3  adjacent
        c.add(gates.CZ(0, 4))  # CZ_4  distance 4 -> the routing gate
        c.add(gates.CZ(2, 3))  # CZ_5  adjacent
        return c

    _LINE = [(0, 1), (1, 2), (2, 3), (3, 4)]

    def test_picks_longest_routing_gate(self):
        import Qdislib as qd

        self.assertEqual(
            qd.find_cut(self._circuit(), coupling_map=self._LINE), ["CZ_4"]
        )
        self.assertEqual(
            qd.find_cut(self._circuit(), coupling_map=self._LINE, max_cuts=1), ["CZ_4"]
        )

    def test_adjacent_gates_are_never_cut(self):
        # On a fully-connected coupling map every gate is distance 1 -> nothing to cut.
        import Qdislib as qd

        full = [(a, b) for a in range(5) for b in range(a + 1, 5)]
        self.assertEqual(qd.find_cut(self._circuit(), coupling_map=full), [])

    def test_respects_layout(self):
        # Remap so that logical (1, 2) become the distant physical pair.
        import Qdislib as qd

        layout = {0: 2, 1: 0, 2: 4, 3: 1, 4: 3}
        self.assertEqual(
            qd.find_cut(
                self._circuit(), coupling_map=self._LINE, layout=layout, max_cuts=1
            ),
            ["CZ_3"],
        )

    def test_respects_forbidden(self):
        import Qdislib as qd

        self.assertEqual(
            qd.find_cut(
                self._circuit(),
                coupling_map=self._LINE,
                forbidden_gates=["CZ_4"],
            ),
            [],
        )

    def test_top_n_ranks_by_distance(self):
        import Qdislib as qd

        ranked = qd.find_cut(self._circuit(), coupling_map=self._LINE, top_n=3)
        self.assertEqual(ranked[0]["cut"], ["CZ_4"])
        self.assertEqual(ranked[0]["kind"], "gate")
        self.assertEqual(ranked[0]["score"], -4)  # -distance (longest first)

    def test_selected_cut_reconstructs(self):
        import Qdislib as qd

        circuit = self._circuit()
        cut = qd.find_cut(circuit, coupling_map=self._LINE, max_cuts=1)
        analytical = qd.analytical_solution(circuit, "ZZZZZ")
        got = qd.gate_cutting(
            circuit, cut, observables="ZZZZZ", software="qibo", shots=1024
        )
        self.assertLess(abs(analytical - got), 0.2)


if __name__ == "__main__":
    unittest.main()
