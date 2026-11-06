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

"""Distributions wider than a dict can hold, kept by fragment."""

import importlib.util
import pickle
import unittest
from unittest import mock

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


def _circuit():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(5)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(-0.7, 2)
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    qc.ry(0.9, 4)  # a qubit on its own: a third, unconnected piece
    return qc


CUTS = {
    "gate": {"gate": ["CZ_4"]},
    "hadamard": {"hadamard": ["CZ_4"]},
    "wire": {"wire": [("CZ_4", "RY_5")]},
    "mixed": {"gate": ["CZ_4"], "wire": [("CX_3", "CZ_4")]},
}


def _exact(qc):
    from qiskit.quantum_info import Statevector

    return {k[::-1]: v for k, v in Statevector(qc).probabilities_dict().items()}


@unittest.skipUnless(QISKIT, "qiskit required")
class VeryWideTest(BaseTimedTestCase):
    def test_truth_of_a_distribution_too_wide_for_len(self):
        import numpy as np

        import Qdislib as qd

        wide = qd.FactoredDistribution(
            [(np.array([1.0, 1.0]), [("out", 0)])], [[0]], {}, 70
        )
        self.assertTrue(wide)
        self.assertEqual(wide["0" * 70], 1.0)


class FactoredTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd

        self.qd, self.qc = qd, _circuit()
        self.want = _exact(self.qc)

    def factored(self, cuts, **options):
        """The factored output of a small circuit, to compare with the exact one."""
        from Qdislib.core.cutting_algorithms import factored

        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 2):
            got = self.qd.circuit_cutting_sampling(self.qc, cuts, shots=None, **options)
        self.assertIsInstance(got, self.qd.FactoredDistribution)
        return got

    def assertClose(self, got, want, places=9):
        for key in set(got) | set(want):
            self.assertAlmostEqual(got.get(key, 0.0), want.get(key, 0.0), places, key)

    def test_every_outcome_for_every_cut_type(self):
        for name, cuts in CUTS.items():
            with self.subTest(cut=name):
                got = self.factored(cuts)
                self.assertEqual(got.width, 5)
                self.assertEqual(len(got), 32)
                for key, p in self.want.items():
                    self.assertAlmostEqual(got[key], p, 9, key)

    def test_marginal_is_the_exact_one(self):
        for name, cuts in CUTS.items():
            for qubits in ([1], [3, 0], [4, 2, 1]):
                with self.subTest(cut=name, qubits=qubits):
                    want = {}
                    for key, p in self.want.items():
                        short = "".join(key[q] for q in qubits)
                        want[short] = want.get(short, 0.0) + p
                    self.assertClose(self.factored(cuts).marginal(qubits), want)

    def test_to_dict_is_what_a_dict_output_gives(self):
        from Qdislib.core.cutting_algorithms import factored

        got = self.factored(CUTS["hadamard"])
        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 24):
            self.assertClose(got.to_dict(), self.want)
        with self.assertRaises(ValueError):  # at the real limit: 5 qubits are too many
            with mock.patch.object(factored, "MAX_DENSE_WIDTH", 2):
                got.to_dict()

    def test_output_marginal_keeps_its_key_order(self):
        qubits = [3, 1, 4]
        want = {}
        for key, p in self.want.items():
            short = "".join(key[q] for q in qubits)
            want[short] = want.get(short, 0.0) + p
        got = self.factored(CUTS["gate"], qubits=qubits)
        self.assertEqual(got.width, 3)
        for key, p in want.items():
            self.assertAlmostEqual(got[key], p, 9, key)

    def test_idle_qubit_is_never_one(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(3)  # qubit 2 has no gate
        qc.h(0)
        qc.cz(0, 1)
        qc.h(1)
        from Qdislib.core.cutting_algorithms import factored

        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 1):
            got = self.qd.circuit_cutting_sampling(qc, {"gate": ["CZ_2"]}, shots=None)
        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 24):
            self.assertClose(got.marginal([0, 2]), {"00": 0.5, "10": 0.5})
        self.assertEqual(got["001"], 0.0)
        self.assertAlmostEqual(got["000"], 0.25, 9)

    def test_most_likely_and_sample(self):
        got = self.factored(CUTS["hadamard"])
        top = got.most_likely(3)
        ranked = sorted(self.want.items(), key=lambda kv: -kv[1])[:3]
        self.assertEqual(list(top), [key for key, _ in ranked])
        for key, p in ranked:
            self.assertAlmostEqual(top[key], p, 9)
        counts = got.sample(50000, seed=4)
        self.assertEqual(sum(counts.values()), 50000)
        for key, p in self.want.items():
            self.assertAlmostEqual(counts.get(key, 0) / 50000, p, delta=0.012)
        self.assertEqual(counts, got.sample(50000, seed=4))

    def test_with_shots_it_estimates_the_same(self):
        from Qdislib.core.cutting_algorithms import factored

        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 2):
            got = self.qd.circuit_cutting_sampling(
                self.qc, CUTS["hadamard"], shots=200000, seed=2
            )
        for key, p in self.want.items():
            self.assertAlmostEqual(got[key], p, delta=0.03)

    def test_it_can_go_to_the_master(self):
        got = self.factored(CUTS["gate"])
        copy = pickle.loads(pickle.dumps(got))
        for key, p in self.want.items():
            self.assertAlmostEqual(copy[key], p, 9)

    def test_bad_reads(self):
        got = self.factored(CUTS["gate"])
        for key in ("01", "012345", "0120"):
            with self.assertRaises(KeyError):
                got[key]
        self.assertEqual(got.get("01", 7), 7)
        for qubits in ([], [0, 0], [9]):
            with self.assertRaises(ValueError):
                got.marginal(qubits)

    def test_error_bars_and_draws_say_so_before_running(self):
        from Qdislib.core.cutting_algorithms import factored

        with mock.patch.object(factored, "MAX_DENSE_WIDTH", 2):
            for options in ({"error": 5}, {"num_samples": 10}):
                with self.assertRaisesRegex(ValueError, "qubits="):
                    self.qd.circuit_cutting_sampling(
                        self.qc, CUTS["gate"], shots=10, **options
                    )

    def test_up_to_the_limit_it_is_a_dict(self):
        got = self.qd.circuit_cutting_sampling(self.qc, CUTS["gate"], shots=None)
        self.assertIs(type(got), dict)
        self.assertClose(got, self.want)


@unittest.skipUnless(QISKIT, "qiskit required")
class WideCircuitTest(BaseTimedTestCase):
    def test_a_30_qubit_circuit(self):
        # two 15-qubit GHZ chains joined by the cut gate: exactly two outcomes
        from qiskit import QuantumCircuit
        import Qdislib as qd

        n = 30
        qc = QuantumCircuit(n)
        qc.h(0)
        for q in range(n - 1):
            qc.cx(q, q + 1)
        cut = f"CX_{n // 2 + 1}"
        got = qd.circuit_cutting_sampling(qc, {"hadamard": [cut]}, shots=None)
        self.assertIsInstance(got, qd.FactoredDistribution)
        self.assertAlmostEqual(got["0" * n], 0.5, 9)
        self.assertAlmostEqual(got["1" * n], 0.5, 9)
        self.assertAlmostEqual(got["01" * (n // 2)], 0.0, 9)
        marginal = got.marginal([0, 29])
        self.assertAlmostEqual(marginal["00"], 0.5, 9)
        self.assertAlmostEqual(marginal["11"], 0.5, 9)
        self.assertEqual(set(got.most_likely(2)), {"0" * n, "1" * n})
        counts = got.sample(1000, seed=1)
        self.assertEqual(set(counts), {"0" * n, "1" * n})


if __name__ == "__main__":
    unittest.main()
