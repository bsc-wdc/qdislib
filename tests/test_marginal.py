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

"""Marginal distributions: ``circuit_cutting_sampling(..., qubits=[...])``
rebuilds the distribution over only some qubits."""

import importlib.util
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

CUTS = {
    "gate": {"gate": ["CZ_4"]},
    "hadamard": {"hadamard": ["CZ_4"]},
    "wire": {"wire": [("CZ_4", "RY_5")]},
    "mixed": {"gate": ["CZ_4"], "wire": [("CX_3", "CZ_4")]},
}
KEPT = ([2], [3, 0], [1, 2, 3], [0, 1, 2, 3])


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


def _marginalise(distribution, qubits):
    out = {}
    for bits, p in distribution.items():
        key = "".join(bits[q] for q in qubits)
        out[key] = out.get(key, 0.0) + p
    return out


def _exact(qc, qubits):
    from qiskit.quantum_info import Statevector

    probs = Statevector(qc).probabilities_dict()  # qubit 0 last
    return _marginalise({k[::-1]: v for k, v in probs.items()}, qubits)


@unittest.skipUnless(QISKIT, "qiskit required")
class MarginalTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        warnings.simplefilter("ignore")
        import Qdislib as qd

        self.qd, self.qc = qd, _circuit()

    def assertClose(self, got, want, places=9):
        for key in set(got) | set(want):
            self.assertAlmostEqual(got.get(key, 0.0), want.get(key, 0.0), places, key)

    def test_exact_marginal_for_every_cut_type(self):
        for name, cuts in CUTS.items():
            for qubits in KEPT:
                with self.subTest(cut=name, qubits=qubits):
                    got = self.qd.circuit_cutting_sampling(
                        self.qc, cuts, shots=None, qubits=qubits
                    )
                    self.assertEqual(len(got), 2 ** len(qubits))
                    self.assertClose(got, _exact(self.qc, qubits))

    def test_sampled_marginal_is_the_full_one_marginalised(self):
        # the same seed: the same counts, so the kept correlators are the same
        qd = self.qd
        for name, cuts in CUTS.items():
            with self.subTest(cut=name):
                qubits = [3, 1]
                full = qd.circuit_cutting_sampling(self.qc, cuts, shots=4000, seed=7)
                got = qd.circuit_cutting_sampling(
                    self.qc, cuts, shots=4000, seed=7, qubits=qubits
                )
                want = _marginalise(full, qubits)
                # both abs-normalised: compare up to that, against the exact
                exact = _exact(self.qc, qubits)
                for key in exact:
                    self.assertAlmostEqual(got[key], exact[key], delta=0.08)
                    self.assertAlmostEqual(got[key], want[key], delta=0.05)

    def test_error_bars_and_draws(self):
        qd = self.qd
        cuts = CUTS["mixed"]
        got = qd.circuit_cutting_sampling(
            self.qc, cuts, shots=4000, seed=3, qubits=[2, 0], error=20
        )
        self.assertEqual(set(got.errors), set(got))
        drawn = qd.circuit_cutting_sampling(
            self.qc, cuts, shots=4000, seed=3, qubits=[2, 0], num_samples=200
        )
        exact = _exact(self.qc, [2, 0])
        for key in exact:
            self.assertAlmostEqual(got[key], exact[key], delta=0.08)
            self.assertAlmostEqual(drawn[key], exact[key], delta=0.15)

    def test_through_every_cut_type_function(self):
        qd = self.qd
        exact = _exact(self.qc, [1, 3])
        for function, cut in (
            (qd.gate_cutting_sampling, ["CZ_4"]),
            (qd.wire_cutting_sampling, [("CZ_4", "RY_5")]),
        ):
            got = function(self.qc, cut, shots=None, qubits=[1, 3])
            self.assertClose(got, exact)

    def test_parameter_values_keep_the_marginal(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit import Parameter

        a = Parameter("a")
        qc = QuantumCircuit(3)
        qc.ry(a, 0)
        qc.cx(0, 1)
        qc.ry(0.5, 2)
        qc.cz(1, 2)
        values = [[0.3], [1.2]]
        got = self.qd.circuit_cutting_sampling(
            qc, {"gate": ["CZ_4"]}, shots=None, qubits=[2], parameter_values=values
        )
        for one, value in zip(got, values):
            self.assertClose(one, _exact(qc.assign_parameters(value), [2]))

    def test_bad_qubits(self):
        qd = self.qd
        for qubits in ([], [0, 0], [4], [-1]):
            with self.assertRaises(ValueError):
                qd.circuit_cutting_sampling(
                    self.qc, CUTS["gate"], shots=None, qubits=qubits
                )


if __name__ == "__main__":
    unittest.main()
