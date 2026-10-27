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

"""A Pauli-sum observable runs one set of fragment circuits per qubit-wise
commuting group of its terms, and every term is read from it."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

_H = [
    (0.5, "ZZZZZZ"),
    (0.3, "ZIZIZI"),
    (-0.2, "XXIIZZ"),
    (0.1, "XXIIII"),
    (0.7, "IIIIII"),
    (0.25, "IZIZIZ"),
    (-0.4, "YYZZII"),
]


def _exact_sum(circuit, hamiltonian):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    state = Statevector(circuit)
    return sum(
        c * float(state.expectation_value(SparsePauliOp(p[::-1])).real)
        for c, p in hamiltonian
    )


class CommutingGroupsTest(BaseTimedTestCase):
    def test_groups_commute_qubit_wise_and_keep_every_term(self):
        from Qdislib.core.cutting_algorithms.mixed_cutting import _commuting_groups

        terms = [t for t in _H if set(t[1]) != {"I"}]
        groups = _commuting_groups(terms)
        self.assertEqual(len(groups), 3)
        self.assertCountEqual([t for _, m in groups for t in m], terms)
        for basis, members in groups:
            for _, term in members:
                for b, p in zip(basis, term):
                    self.assertTrue(p == "I" or p == b, (basis, term))


@unittest.skipUnless(QISKIT, "qiskit required")
class GroupedObservableTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        from tests.test_engine_cache import _CUT, _circuit

        self.qc, self.cut = _circuit(), _CUT

    def test_exact_on_every_path(self):
        import Qdislib as qd

        ref = _exact_sum(self.qc, _H)
        for name, value in (
            ("circuit_cutting", qd.circuit_cutting(self.qc, self.cut, _H, shots=None)),
            (
                "optimized",
                qd.circuit_cutting(self.qc, self.cut, _H, shots=None, optimized=True),
            ),
            ("gate", qd.gate_cutting(self.qc, ["CZ_9"], _H, shots=None)),
            ("hadamard", qd.gate_cutting_hadamard(self.qc, ["CZ_9"], _H, shots=None)),
            ("wire", qd.wire_cutting(self.qc, [("CZ_8", "CZ_9")], _H, shots=None)),
        ):
            with self.subTest(path=name):
                self.assertAlmostEqual(value, ref, places=9)
        subs = qd.circuit_cutting_subcircuits(self.qc, self.cut, _H)
        values = [qd.expectation_value(s, shots=None) for s in subs]
        self.assertAlmostEqual(subs.reconstruct(values), ref, places=9)

    def test_one_run_per_group(self):
        import Qdislib as qd
        from tests.test_engine_cache import _Runs

        per_term = qd.cutting_overhead(self.qc, self.cut)["subcircuits"]
        with _Runs() as runs:
            qd.circuit_cutting(self.qc, self.cut, _H, shots=None)
        # 6 non-identity terms in 3 groups; the identity term runs nothing
        self.assertEqual(runs.count, 3 * per_term)
        # three-step: the identity term is read from the first group's circuits
        self.assertEqual(
            len(qd.circuit_cutting_subcircuits(self.qc, self.cut, _H)), 3 * per_term
        )

    def test_shots_and_num_samples(self):
        import Qdislib as qd

        ref = _exact_sum(self.qc, _H)
        value = qd.circuit_cutting(self.qc, self.cut, _H, shots=20000, seed=1)
        self.assertLess(abs(value - ref), 0.05)
        mean = np.mean(
            [
                qd.circuit_cutting(
                    self.qc, self.cut, _H, shots=200, seed=s, num_samples=50
                )
                for s in range(10)
            ]
        )
        self.assertLess(abs(mean - ref), 0.05)


if __name__ == "__main__":
    unittest.main()
