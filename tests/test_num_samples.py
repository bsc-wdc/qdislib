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
"""num_samples: configurations drawn with probability |weight| / gamma instead of all
run (Monte-Carlo quasiprobability sampling), for every kind of cut."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


@unittest.skipUnless(QISKIT, "qiskit not installed")
class NumSamplesTest(BaseTimedTestCase):
    def test_unbiased_for_every_kind_of_cut(self):
        """Exact fragments, so only the draws are random: the mean of independent
        runs sits within a few standard errors of the exact value."""
        import Qdislib as qd
        from tests.test_circuit_cutting import (
            _all_kinds,
            _chain,
            _exact,
            _hadamard_group,
        )

        for make in (_all_kinds, _chain, _hadamard_group):
            circuit, cuts = make()
            exact = _exact(circuit, "Z" * circuit.num_qubits)
            for optimized in (False, True):
                with self.subTest(case=make.__name__, optimized=optimized):
                    values = [
                        qd.circuit_cutting(
                            circuit,
                            cuts,
                            shots=None,
                            optimized=optimized,
                            num_samples=1500,
                            seed=seed,
                        )
                        for seed in range(12)
                    ]
                    error = np.std(values) / np.sqrt(len(values)) + 1e-9
                    self.assertLess(abs(np.mean(values) - exact), 4 * error)

    def test_hadamard_draws_follow_the_weights(self):
        """A Hadamard setting (s, t, g) is drawn with probability c_s c_t / phi."""
        from qiskit import QuantumCircuit

        from Qdislib.core.cutting_algorithms.hadamard_cut import _cut_pmf
        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.h(1)
        qc.rzz(0.4, 0, 1)  # RZZ_3: very unequal terms
        plan = _CutPlan(qc, {"hadamard": ["RZZ_3"]}, None, False, False, True)
        settings, probs, _ = _cut_pmf(plan.groups[0][1])
        draws = plan.draw(200_000, seed=1)
        counts = {config[0]: count for _, config, _, _, count in draws}
        for setting, p in zip(settings, probs):
            self.assertAlmostEqual(counts.get(setting, 0) / 200_000, p, delta=0.005)
        weights = {config[0]: w for _, config, _, w, _ in draws}
        gamma = plan.gamma()
        for (s, t, g), w in weights.items():
            self.assertAlmostEqual(w, gamma * (-1) ** g / 200_000)

    def test_reproducible_and_counts_add_up(self):
        import Qdislib as qd
        from tests.test_circuit_cutting import _chain

        circuit, cuts = _chain()
        a = qd.circuit_cutting(
            circuit, cuts, backend="qiskit", shots=4, num_samples=300, seed=7
        )
        b = qd.circuit_cutting(
            circuit, cuts, backend="qiskit", shots=4, num_samples=300, seed=7
        )
        self.assertEqual(a, b)
        subs = qd.circuit_cutting_subcircuits(
            circuit, cuts, num_samples=300, seed=7, optimized=True
        )
        fragments = len(subs.plan.fragments)
        # every draw reads one circuit per fragment (more with |-> / |-i> preps)
        self.assertGreaterEqual(sum(s.draws for s in subs), 300 * fragments)
        self.assertLessEqual(len(subs), subs.plan.num_units())

    def test_gate_and_wire_cutting_take_it(self):
        import Qdislib as qd
        from tests.test_standard_api import CUTS, _circuit, _exact

        circuit = _circuit()
        for family in ("gate_cutting", "wire_cutting"):
            with self.subTest(family=family):
                value = getattr(qd, family)(
                    circuit, CUTS[family], "XZZY", shots=None, num_samples=4000, seed=3
                )
                self.assertLess(abs(value - _exact(circuit, "XZZY")), 0.2)

    def test_rejects_a_bad_count(self):
        import Qdislib as qd
        from tests.test_circuit_cutting import _chain

        circuit, cuts = _chain()
        with self.assertRaises(ValueError):
            qd.circuit_cutting(circuit, cuts, num_samples=0)


class LargeJointGroupTest(BaseTimedTestCase):
    """A joint Hadamard group's 2**k product terms are never listed."""

    def test_product_terms_match_the_explicit_product(self):
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            _ProductTerms,
            local_unitary_decomposition,
            phi_factor,
        )

        decomps = [
            local_unitary_decomposition("RZZ", (0.3,)),
            local_unitary_decomposition("CZ"),
            local_unitary_decomposition("CP", (1.1,)),
        ]
        explicit = [(1.0, [])]
        for terms in decomps:
            explicit = [
                (c * c_m, factors + [(v, w)])
                for c, factors in explicit
                for c_m, v, w in terms
            ]
        lazy = _ProductTerms(decomps)
        self.assertEqual(len(lazy), len(explicit))
        for s, (c, factors) in enumerate(explicit):
            c_lazy, factors_lazy = lazy[s]
            self.assertAlmostEqual(c_lazy, c)
            for (v, w), (v2, w2) in zip(factors, factors_lazy):
                self.assertTrue(np.allclose(v, v2) and np.allclose(w, w2))
        self.assertAlmostEqual(phi_factor(lazy), phi_factor(explicit))

    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_overhead_of_a_twelve_gate_group_is_instant(self):
        import time

        import Qdislib as qd
        from tests.test_circuit_cutting import _boundary

        circuit, cuts = _boundary(12, "hadamard")
        start = time.time()
        cost = qd.cutting_overhead(circuit, cuts, joint=True)
        self.assertLess(time.time() - start, 5)
        self.assertAlmostEqual(cost["gamma"], 2**13 - 1)
        self.assertEqual(cost["configurations"], 2**12 * (2**13 - 1))


if __name__ == "__main__":
    unittest.main()
