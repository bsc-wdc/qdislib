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

"""Readable views of the three-step objects: the list, a subcircuit, a result."""

import importlib.util
import pickle
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


@unittest.skipUnless(QISKIT, "qiskit required")
class DescribeTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from tests.test_engine_cache import _circuit

        self.qd, self.qc = qd, _circuit()

    def test_list_and_subcircuit_reprs(self):
        qd = self.qd
        subs = qd.circuit_cutting_subcircuits(
            self.qc, {"gate": ["CZ_9"], "wire": [("CZ_7", "CZ_8")]}, "ZZZZZZ"
        )
        self.assertIn("subcircuits, 3 fragments", repr(subs))
        self.assertIn("gate ['CZ_9']", repr(subs))
        text = [repr(s) for s in subs]
        self.assertTrue(any("CZ_9 side 1" in t for t in text))
        self.assertTrue(any("prepare |+⟩" in t for t in text))
        self.assertTrue(any("q1'" in t for t in text))  # the wire's prepared line

    def test_describe_names_values_and_reconstructs(self):
        qd = self.qd
        subs = qd.circuit_cutting_subcircuits(self.qc, {"gate": ["CZ_9"]}, "ZZZZZZ")
        values = [qd.expectation_value(s, shots=None) for s in subs]
        table = subs.describe(values)
        self.assertEqual(len(table.rows), len(subs))
        first = table.rows[0]
        self.assertEqual(set(first), {"#", "fragment", "qubits", "terms", "values"})
        self.assertEqual(list(first["values"]), ["Z[q0,q1,q2]"])
        self.assertIn(f"{subs.reconstruct(values):+.10g}", repr(table))
        self.assertIn("<table>", table._repr_html_())
        self.assertNotIn("values", subs.describe().rows[0])  # without results

    def test_a_pauli_sum_in_its_bases(self):
        qd = self.qd
        hamiltonian = [(0.5, "ZZZZZZ"), (-0.2, "XXIIII"), (0.1, "IIIYZX")]
        subs = qd.circuit_cutting_subcircuits(self.qc, {"gate": ["CZ_9"]}, hamiltonian)
        values = [qd.expectation_value(s, shots=None) for s in subs]
        rows = subs.describe(values).rows
        labels = {label for row in rows for label in row["values"]}
        self.assertIn("X[q0,q1]", labels)
        self.assertIn("Y[q3] Z[q4] X[q5]", labels)
        self.assertEqual({row["group"] for row in rows}, {0, 1})

    def test_distribution_results_are_summarised(self):
        qd = self.qd
        subs = qd.circuit_cutting_sampling_subcircuits(self.qc, {"gate": ["CZ_9"]})
        values = [qd.sample_subcircuit(s, shots=None) for s in subs]
        self.assertIn("all 2**3 correlators", repr(values[0]))
        self.assertEqual(list(values[0].labelled()), ["Z[q0]", "Z[q1]", "Z[q2]"])
        self.assertIn("most likely", repr(subs.describe(values)))

    def test_a_result_is_still_a_pair(self):
        qd = self.qd
        subs = qd.circuit_cutting_subcircuits(self.qc, {"gate": ["CZ_9"]}, "ZZZZZZ")
        result = qd.expectation_value(subs[0], shots=1000, seed=1)
        corr, qubits = result
        self.assertIs(corr, result.corr)
        self.assertEqual(qubits, [0, 1, 2])
        again = pickle.loads(pickle.dumps(result))
        self.assertEqual(type(again).__name__, "FragmentResult")
        self.assertEqual(again.labelled(), result.labelled())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            store = qd.InMemoryCircuitCache()
            qd.expectation_value(subs[0], shots=1000, cache=store)
            cached = qd.expectation_value(subs[0], shots=1000, cache=store)
        self.assertIn("shots, any correlator", repr(cached))


if __name__ == "__main__":
    unittest.main()
