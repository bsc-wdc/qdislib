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

"""Reading counts is the work of a reconstruction: one task per fragment."""

import importlib.util
import unittest
from unittest import mock

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


def _circuit():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(-0.7, 2)
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    return qc


def _exact(qc, observable):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return Statevector(qc).expectation_value(SparsePauliOp(observable[::-1])).real


CUTS = {
    "gate": {"gate": ["CZ_4"]},
    "hadamard": {"hadamard": ["CZ_4"]},
    "wire": {"wire": [("CZ_4", "RY_5")]},
}


@unittest.skipUnless(QISKIT, "qiskit required")
class FragmentTasksTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms import mixed_cutting as mc

        self.qd, self.mc, self.qc = qd, mc, _circuit()

    def spy(self, name):
        return mock.patch.object(self.mc, name, wraps=getattr(self.mc, name))

    def test_error_bars_read_each_fragment_in_its_task(self):
        for name, cuts in CUTS.items():
            with self.subTest(cut=name), self.spy("_fragment_tensor_task") as fragment:
                with self.spy("_reconstruct_cut_task") as single:
                    got = self.qd.circuit_cutting(
                        self.qc, cuts, "ZZZZ", shots=40000, seed=1, error=20
                    )
                self.assertEqual(single.call_count, 0)
                self.assertGreaterEqual(fragment.call_count, 2)
                self.assertAlmostEqual(got, _exact(self.qc, "ZZZZ"), delta=0.1)
                self.assertGreater(got.error, 0)
                self.assertGreater(got.shots, 0)

    def test_a_few_correlators_stay_one_task(self):
        with self.spy("_fragment_tensor_task") as fragment:
            with self.spy("_reconstruct_cut_task") as single:
                self.qd.circuit_cutting(
                    self.qc, CUTS["gate"], "ZZZZ", shots=1000, seed=1
                )
        self.assertEqual(fragment.call_count, 0)
        self.assertEqual(single.call_count, 1)

    def test_the_error_is_the_repeated_runs_spread(self):
        import numpy as np

        qd, qc, cuts = self.qd, self.qc, CUTS["hadamard"]
        one = qd.circuit_cutting(qc, cuts, "ZZZZ", shots=2000, seed=11, error=60)
        runs = [
            float(qd.circuit_cutting(qc, cuts, "ZZZZ", shots=2000, seed=100 + n))
            for n in range(60)
        ]
        spread = float(np.std(runs, ddof=1))
        self.assertLess(abs(one.error - spread), 0.5 * spread)

    def test_a_group_of_terms_with_error_bars(self):
        observable = {"ZZII": 0.5, "IZZI": -1.0, "ZZZZ": 2.0}
        want = sum(c * _exact(self.qc, term) for term, c in observable.items())
        got = self.qd.circuit_cutting(
            self.qc, CUTS["gate"], observable, shots=40000, seed=2, error=10
        )
        self.assertAlmostEqual(got, want, delta=0.2)
        self.assertGreater(got.error, 0)

    def test_distribution_with_error_bars(self):
        from qiskit.quantum_info import Statevector

        want = {
            k[::-1]: v for k, v in Statevector(self.qc).probabilities_dict().items()
        }
        with self.spy("_reconstruct_cut_task") as single:
            got = self.qd.circuit_cutting_sampling(
                self.qc, CUTS["wire"], shots=40000, seed=3, error=10
            )
        self.assertEqual(single.call_count, 0)
        self.assertEqual(set(got.errors), set(got))
        self.assertGreater(max(got.errors.values()), 0)
        for key, p in want.items():
            self.assertAlmostEqual(got.get(key, 0.0), p, delta=0.05)

    def test_marginal_with_error_bars(self):
        got = self.qd.circuit_cutting_sampling(
            self.qc, CUTS["gate"], shots=20000, seed=3, error=10, qubits=[2, 0]
        )
        self.assertEqual(len(got), 4)
        self.assertEqual(set(got.errors), set(got))

    def test_exact_runs_have_no_error(self):
        got = self.qd.circuit_cutting(
            self.qc, CUTS["gate"], "ZZZZ", shots=None, error=5
        )
        self.assertAlmostEqual(got, _exact(self.qc, "ZZZZ"), places=9)
        self.assertEqual(got.error, 0.0)

    def test_drawn_configurations_stay_one_task(self):
        with self.spy("_reconstruct_cut_task") as single:
            got = self.qd.circuit_cutting(
                self.qc,
                CUTS["gate"],
                "ZZZZ",
                shots=20,
                seed=4,
                num_samples=4000,
                error=5,
            )
        self.assertEqual(single.call_count, 1)
        self.assertGreater(got.error, 0)

    def test_counts_brought_by_hand_go_fragment_by_fragment(self):
        subs = self.qd.circuit_cutting_subcircuits(self.qc, CUTS["gate"], "ZZZZ")
        results = [
            self.qd.sample_subcircuit(s, shots=40000, seed=n)
            for n, s in enumerate(subs)
        ]
        with self.spy("_fragment_tensor_task") as fragment:
            got = subs.reconstruct(results)
        self.assertGreaterEqual(fragment.call_count, 2)
        self.assertAlmostEqual(got, _exact(self.qc, "ZZZZ"), delta=0.1)


if __name__ == "__main__":
    unittest.main()
