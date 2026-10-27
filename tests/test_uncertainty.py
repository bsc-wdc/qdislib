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

"""Error bars (bootstrap over the shots), total_shots= and shots_for."""

import importlib.util
import math
import pickle
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

_CUTS = {"gate": ["CZ_9"], "wire": [("CZ_7", "CZ_8")]}
_H = [(0.5, "ZZZZZZ"), (0.3, "ZIZIZI"), (-0.2, "XXIIZZ"), (0.7, "IIIIII")]


class EstimateTest(BaseTimedTestCase):
    def test_a_float_with_an_error(self):
        from Qdislib import Estimate
        from Qdislib.core.cutting_algorithms.uncertainty import combine

        value = Estimate(0.25, 0.01, 1000)
        self.assertEqual(value + 1, 1.25)
        self.assertEqual((value.error, value.shots), (0.01, 1000))
        self.assertEqual(repr(value), "0.25 ± 0.01")
        again = pickle.loads(pickle.dumps(value))
        self.assertEqual((float(again), again.error, again.shots), (0.25, 0.01, 1000))
        total = combine([2.0, 1.0, 1.0], [value, Estimate(1.0, 0.03, 10), 0.5])
        self.assertAlmostEqual(float(total), 2.0)
        self.assertAlmostEqual(total.error, math.hypot(0.02, 0.03))
        self.assertEqual(combine([1.0], [0.5]), 0.5)  # no error: a plain float


@unittest.skipUnless(QISKIT, "qiskit required")
class ErrorBarTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from qiskit.quantum_info import SparsePauliOp, Statevector
        from tests.test_engine_cache import _circuit

        self.qd, self.qc = qd, _circuit()
        state = Statevector(self.qc)
        self.ref = sum(
            c * float(state.expectation_value(SparsePauliOp(p[::-1])).real)
            for c, p in _H
        )

    def _spread(self, runs=20, **kwargs):
        values = [
            self.qd.circuit_cutting(self.qc, _CUTS, _H, seed=s, **kwargs)
            for s in range(runs)
        ]
        return float(np.std(values, ddof=1))

    def test_exact_has_no_error(self):
        value = self.qd.circuit_cutting(self.qc, _CUTS, _H, shots=None, error=True)
        self.assertIsInstance(value, self.qd.Estimate)
        self.assertEqual(value.error, 0.0)
        self.assertAlmostEqual(float(value), self.ref, places=9)

    def test_the_error_matches_repeated_runs(self):
        for kwargs in (
            {"shots": 2000},
            {"shots": 2000, "optimized": True},
            {"shots": 300, "num_samples": 100},  # every configuration fits
        ):
            with self.subTest(**kwargs):
                value = self.qd.circuit_cutting(
                    self.qc, _CUTS, _H, error=True, seed=99, **kwargs
                )
                spread = self._spread(**kwargs)
                self.assertGreater(value.error, 0.6 * spread)
                self.assertLess(value.error, 1.6 * spread)

    def test_every_cut_function_takes_it(self):
        qd = self.qd
        for call in (
            lambda **kw: qd.gate_cutting(self.qc, ["CZ_9"], "ZZZZZZ", **kw),
            lambda **kw: qd.wire_cutting(self.qc, [("CZ_7", "CZ_8")], "ZZZZZZ", **kw),
            lambda **kw: qd.gate_cutting_hadamard(self.qc, ["CZ_9"], "ZZZZZZ", **kw),
        ):
            value = call(shots=1000, seed=1, error=50)
            self.assertGreater(value.error, 0)
            self.assertGreater(value.shots, 0)


@unittest.skipUnless(QISKIT, "qiskit required")
class BudgetTest(BaseTimedTestCase):
    def test_the_budget_is_spent_and_beats_even_shots(self):
        import Qdislib as qd
        from qiskit.quantum_info import SparsePauliOp, Statevector
        from tests.test_engine_cache import _circuit

        qc = _circuit()
        ref = float(Statevector(qc).expectation_value(SparsePauliOp("ZZZZZZ")).real)
        cuts = {"hadamard": ["CZ_9"]}
        units = len(qd.circuit_cutting_subcircuits(qc, cuts, "ZZZZZZ"))
        budget = 400 * units
        value = qd.circuit_cutting(
            qc, cuts, "ZZZZZZ", total_shots=budget, seed=1, error=True
        )
        self.assertEqual(value.shots, budget)

        def rmse(**kwargs):
            values = [
                qd.circuit_cutting(qc, cuts, "ZZZZZZ", seed=s, **kwargs)
                for s in range(40)
            ]
            return math.sqrt(np.mean((np.array(values) - ref) ** 2))

        # weights differ (a Hadamard cut): the budget follows them
        self.assertLess(rmse(total_shots=budget), 1.05 * rmse(shots=budget // units))

    def test_budget_options(self):
        import Qdislib as qd
        from tests.test_engine_cache import _circuit

        qc = _circuit()
        with self.assertRaisesRegex(ValueError, "fewer than"):
            qd.circuit_cutting(qc, _CUTS, "ZZZZZZ", total_shots=5)
        with self.assertRaisesRegex(ValueError, "cache"):
            qd.circuit_cutting(qc, _CUTS, "ZZZZZZ", total_shots=10**5, cache="memory")
        # with num_samples the budget is whole shots per draw and fragment: at most
        # the budget, rounded down
        value = qd.circuit_cutting(
            qc, _CUTS, "ZZZZZZ", total_shots=20000, num_samples=200, seed=1, error=True
        )
        self.assertLessEqual(value.shots, 20000)
        self.assertGreater(value.shots, 0.95 * 20000)

    def test_shots_for(self):
        import Qdislib as qd
        from tests.test_engine_cache import _circuit

        qc = _circuit()
        gamma = qd.cutting_overhead(qc, _CUTS)["gamma"]
        self.assertEqual(
            qd.shots_for(qc, _CUTS, 0.01, 0.05),
            math.ceil(2 * gamma**2 / 0.01**2 * math.log(2 / 0.05)),
        )


if __name__ == "__main__":
    unittest.main()
