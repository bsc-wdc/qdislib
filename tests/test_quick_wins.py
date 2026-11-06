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
    def test_a_membership_test_is_not_a_lookup(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        cache.put("k", 1.0)
        self.assertIn("k", cache)
        self.assertNotIn("absent", cache)
        self.assertEqual(cache.stats()["hits"], 0)
        self.assertEqual(cache.stats()["misses"], 0)

    def test_a_cut_run_counts_each_circuit_once(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.3, 1)
        qc.cz(1, 2)
        qc.ry(0.5, 2)
        store = InMemoryCircuitCache()
        first = qd.gate_cutting(qc, ["CZ_2"], "ZZZ", shots=200, seed=1, cache=store)
        stats = store.stats()
        self.assertEqual(stats["misses"], stats["size"])  # one per circuit run
        self.assertEqual(stats["hits"], 0)  # circuits it repeats are not looked up
        second = qd.gate_cutting(qc, ["CZ_2"], "ZZZ", shots=200, seed=1, cache=store)
        self.assertEqual(first, second)
        again = store.stats()
        self.assertEqual(again["misses"], stats["misses"])
        self.assertGreaterEqual(again["hits"], stats["size"])

    def test_a_semantic_run_counts_each_circuit_once(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.3, 1)
        qc.cz(1, 2)
        qc.ry(0.5, 2)
        store = InMemoryCircuitCache()
        qd.gate_cutting(
            qc, ["CZ_2"], "ZZZ", shots=None, cache=store, cache_key="semantic"
        )
        stats = store.stats()
        self.assertEqual(stats["misses"], stats["size"])

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
        self.assertIsNone(metrics["cache"])

    def test_metrics_with_a_cache_named_as_cache_takes_it(self):
        import Qdislib as qd

        for cache in (True, "memory", qd.InMemoryCircuitCache()):
            _, metrics = qd.estimate(
                _ref_circuit(),
                "ZZ",
                gate_cut=True,
                wire_cut=False,
                shots=100,
                cache=cache,
                return_metrics=True,
            )
            self.assertGreater(metrics["cache"]["size"], 0)

    def test_warns_only_about_what_was_passed(self):
        import warnings

        import Qdislib as qd

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            qd.estimate(_ref_circuit(), "ZZ", gate_cut=True, wire_cut=False, shots=100)
        self.assertFalse([w for w in caught if "cache_mode" in str(w.message)])


@unittest.skipUnless(QIBO_AVAILABLE and MPL_AVAILABLE, "qibo and matplotlib required")
class DrawCutTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import matplotlib

        matplotlib.use("Agg")  # headless: render to file, no display needed

    def test_draw_to_file(self):
        import Qdislib as qd

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cut.png")
            qd.draw_cut(_ref_circuit(), ["CZ_2"], filename=path)
            self.assertTrue(os.path.exists(path))
            self.assertGreater(os.path.getsize(path), 0)

    def test_every_kind_of_cut(self):
        """Gate, Hadamard and wire cuts, as a list or a dict, into a given Axes."""
        import matplotlib.pyplot as plt

        import Qdislib as qd
        from tests.test_circuit_cutting import _Builder

        b = _Builder(5, seed=1)
        b.layer(0, 1, 2, 3, 4)
        b.gate("cz", 0, 1)
        had = b.gate("cz", 1, 2)
        b.layer(0, 1, 2)
        mp = b.gate("cx", 2, 4)
        b.layer(0, 1, 2, 4)
        wires = [b.wire(0, 3), b.wire(1, 3)]
        b.layer(0, 1, 3)
        for cuts in (
            {"gate": [mp], "hadamard": [had], "wire": wires},
            [mp, *wires],
            [had],
        ):
            with self.subTest(cuts=cuts):
                fig, ax = plt.subplots()
                out = qd.draw_cut(b.qc, cuts, ax=ax, hadamard=cuts == [had])
                self.assertIs(out, ax)
                self.assertIn("fragment", ax.get_title(loc="left"))
                plt.close(fig)
        fig, ax = plt.subplots()
        qd.draw_cut(b.qc, {"gate": [mp], "hadamard": [had], "wire": wires}, ax=ax)
        self.assertTrue(ax.get_title(loc="left").startswith("4 fragments"))
        plt.close(fig)

    def test_an_invalid_cut_is_drawn_not_raised(self):
        import matplotlib.pyplot as plt

        import Qdislib as qd
        from tests.test_circuit_cutting import _Builder

        b = _Builder(3, seed=8)
        b.layer(0, 1, 2)
        b.gate("cz", 0, 1)
        gate = b.gate("cz", 0, 1)
        b.layer(0, 1)
        wire = b.wire(1, 2)
        b.layer(1, 2)
        fig, ax = plt.subplots()
        qd.draw_cut(b.qc, {"gate": [gate], "wire": [wire]}, ax=ax)
        self.assertIn("Not a valid cut", ax.get_title(loc="left"))
        self.assertIn(gate, ax.get_title(loc="left"))
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
