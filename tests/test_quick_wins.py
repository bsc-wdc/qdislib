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

"""Tests for the quick-win features: cache stats, from_qasm, estimate, draw_cut."""

import importlib.util
import os
import tempfile
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None
MPL_AVAILABLE = importlib.util.find_spec("matplotlib") is not None


def _ref_circuit():
    """H(0), CZ(0,1), RY(0, 0.8): analytical <ZZ> ~ -0.717 (clearly non-zero)."""
    from qibo import models, gates

    circ = models.Circuit(2)
    circ.add(gates.H(0))
    circ.add(gates.CZ(0, 1))
    circ.add(gates.RY(0, 0.8))
    return circ


class CacheStatsTest(BaseTimedTestCase):
    def test_stats_track_hits_and_misses(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        self.assertEqual(
            cache.stats(), {"hits": 0, "misses": 0, "hit_rate": 0.0, "size": 0}
        )
        cache.put("k", 1.0)
        cache.get("k")  # hit
        cache.get("k")  # hit
        cache.get("absent")  # miss
        stats = cache.stats()
        self.assertEqual(stats["hits"], 2)
        self.assertEqual(stats["misses"], 1)
        self.assertAlmostEqual(stats["hit_rate"], 2 / 3)
        self.assertEqual(stats["size"], 1)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class FromQasmTest(BaseTimedTestCase):
    def test_roundtrip_qibo(self):
        import Qdislib as qd

        qc = qd.from_qasm(_ref_circuit().to_qasm(), software="qibo")
        self.assertEqual(qc.nqubits, 2)

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_roundtrip_qiskit(self):
        import Qdislib as qd

        qc = qd.from_qasm(_ref_circuit().to_qasm(), software="qiskit")
        self.assertEqual(qc.num_qubits, 2)

    def test_unknown_software_raises(self):
        import Qdislib as qd

        with self.assertRaises(ValueError):
            qd.from_qasm("OPENQASM 2.0;", software="nonsense")


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class EstimateTest(BaseTimedTestCase):
    def test_matches_analytical(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_ref_circuit(), "ZZ")
        self.assertGreater(abs(analytical), 0.5)  # guard: non-zero reference
        value = qd.estimate(
            _ref_circuit(), "ZZ", gate_cut=True, wire_cut=False, shots=20000
        )
        self.assertLess(abs(value - analytical), 0.2)

    def test_return_metrics(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        value, metrics = qd.estimate(
            _ref_circuit(),
            "ZZ",
            gate_cut=True,
            wire_cut=False,
            shots=8000,
            return_metrics=True,
        )
        self.assertIsInstance(value, float)
        self.assertEqual(metrics["cut_type"], "gate")
        self.assertEqual(metrics["num_cuts"], len(metrics["cut"]))
        self.assertIn("find_cut_seconds", metrics)
        self.assertIn("run_seconds", metrics)


@unittest.skipUnless(QIBO_AVAILABLE and MPL_AVAILABLE, "qibo and matplotlib required")
class DrawCutTest(BaseTimedTestCase):
    def test_draw_to_file(self):
        import matplotlib

        matplotlib.use("Agg")  # headless: render to file, no display needed
        import Qdislib as qd

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cut.png")
            qd.draw_cut(_ref_circuit(), ["CZ_2"], filename=path)
            self.assertTrue(os.path.exists(path))
            self.assertGreater(os.path.getsize(path), 0)


if __name__ == "__main__":
    unittest.main()
