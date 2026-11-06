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

"""Dense ``unitary`` gates (Qiskit's ``UnitaryGate``) in a cut circuit."""

import importlib.util
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None
QIBO = importlib.util.find_spec("qibo") is not None
CUTS = (
    {"hadamard": ["CX_2"]},
    {"gate": ["CX_2"]},
    {"wire": [("CX_2", "UNITARY_3")]},
)


def _circuit():
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import UnitaryGate
    from qiskit.quantum_info import random_unitary

    qc = QuantumCircuit(6)
    qc.append(UnitaryGate(random_unitary(8, seed=1)), [0, 2, 1])
    qc.cx(2, 3)
    qc.append(UnitaryGate(random_unitary(8, seed=2)), [5, 3, 4])
    qc.ry(0.3, 0)
    qc.append(UnitaryGate(random_unitary(4, seed=3)), [1, 0])
    return qc


def _exact(qc, qubits):
    from qiskit.quantum_info import Statevector

    out = {}
    for key, p in Statevector(qc).probabilities_dict().items():
        bits = key[::-1]
        marginal = "".join(bits[q] for q in qubits)
        out[marginal] = out.get(marginal, 0.0) + p
    return out


@unittest.skipUnless(QISKIT, "qiskit required")
class UnitaryGateTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        warnings.simplefilter("ignore")
        import Qdislib as qd

        self.qd, self.qc = qd, _circuit()

    def check(self, got, want, delta):
        for key, p in want.items():
            self.assertAlmostEqual(got.get(key, 0.0), p, delta=delta)

    def test_round_trip_keeps_the_circuit(self):
        from qiskit.quantum_info import Operator

        from Qdislib.utils.graph_qiskit import (
            circuit_qiskit_to_dag,
            dag_to_circuit_qiskit,
        )

        dag = circuit_qiskit_to_dag(self.qc)
        back = dag_to_circuit_qiskit(dag, self.qc.num_qubits)[0]
        self.assertTrue(Operator(back).equiv(Operator(self.qc)))

    def test_exact_every_cut_type(self):
        want = _exact(self.qc, [5, 0, 3])
        for cuts in CUTS:
            with self.subTest(cuts=cuts):
                got = self.qd.circuit_cutting_sampling(
                    self.qc, cuts, shots=None, qubits=[5, 0, 3]
                )
                self.check(got, want, 1e-9)

    def test_sampled_on_qiskit_and_qibo(self):
        want = _exact(self.qc, [5, 0, 3])
        backends = ["qiskit"] + (["qibo"] if QIBO else [])
        for backend in backends:
            with self.subTest(backend=backend):
                got = self.qd.circuit_cutting_sampling(
                    self.qc,
                    CUTS[0],
                    shots=50_000,
                    seed=1,
                    backend=backend,
                    qubits=[5, 0, 3],
                )
                self.check(got, want, 0.01)


@unittest.skipUnless(QISKIT, "qiskit required")
class DenseTwoQubitGateTest(BaseTimedTestCase):
    """A dense two-qubit unitary as the gate cut, and in the cut search."""

    def _circuit(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate
        from qiskit.quantum_info import random_unitary

        qc = QuantumCircuit(2)
        qc.ry(0.3, 0)
        qc.ry(0.2, 1)
        qc.append(UnitaryGate(random_unitary(4, seed=1)), [0, 1])
        qc.ry(0.3, 0)
        return qc

    def test_hadamard_cut_of_it(self):
        import Qdislib as qd

        qc = self._circuit()
        exact = qd.analytical_solution(qc, "ZZ")
        cut = ["UNITARY_3"]
        self.assertAlmostEqual(
            qd.gate_cutting_hadamard(qc, cut, "ZZ", shots=None), exact
        )
        value = qd.gate_cutting_hadamard(qc, cut, "ZZ", shots=20000, seed=1)
        self.assertAlmostEqual(value, exact, delta=0.2)

    def test_the_mixed_search_prices_it(self):
        import Qdislib as qd

        qc = self._circuit()
        exact = qd.analytical_solution(qc, "ZZ")
        cuts = qd.find_cut(qc, max_qubits=1, mixed=True, seed=0, time_limit=5)
        self.assertEqual(sum(len(c) for c in cuts.values()), 1)
        self.assertAlmostEqual(qd.circuit_cutting(qc, cuts, "ZZ", shots=None), exact)

    def test_a_cheaper_dense_gate_costs_less(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import UnitaryGate

        from Qdislib.core.find_cut.cost_search import _gate_costs
        from Qdislib.utils.backends import convert_to_dag

        cz = QuantumCircuit(2)
        cz.append(
            UnitaryGate([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]),
            [0, 1],
        )
        dense = _gate_costs(convert_to_dag(cz), "UNITARY_1", {})
        general = _gate_costs(convert_to_dag(self._circuit()), "UNITARY_3", {})
        self.assertLess(dense[1], general[1])


if __name__ == "__main__":
    unittest.main()
