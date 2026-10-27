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

"""What a cut costs: golden wire cut points (the Y terms a real state leaves out)
and the warning that points from a measure-and-prepare cut of a rotation to the
cheaper Hadamard cut."""

import importlib.util
import unittest
import warnings

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None
AER = importlib.util.find_spec("qiskit_aer") is not None


def _exact(circuit, pauli):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(pauli[::-1])).real
    )


def _chain(real=True):
    """Three blocks joined by two wire cuts in series: qubit 1 moves from block A
    to block B, qubit 3 from B to C. ``real=False`` adds an RX (a complex gate)
    in block A."""
    from qiskit import QuantumCircuit

    rng = np.random.default_rng(0)
    qc = QuantumCircuit(6)
    for q in range(6):
        qc.ry(rng.uniform(0.3, 2.8), q)
    if not real:
        qc.rx(0.7, 1)
    qc.cz(0, 1)
    qc.cz(2, 3)
    qc.cz(4, 5)
    qc.cz(1, 2)
    qc.ry(0.4, 1)
    qc.ry(0.9, 2)
    qc.cz(1, 2)
    qc.ry(0.5, 0)
    qc.cz(3, 4)
    qc.ry(1.1, 3)
    qc.ry(0.6, 4)
    qc.ry(1.3, 5)
    qc.cz(4, 5)
    qc.cz(3, 4)
    return qc


def _chain_cuts(qc):
    from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

    dag = circuit_qiskit_to_dag(qc)
    # The first CZ on each pair: CZ(0,1) -> CZ(1,2) on qubit 1, CZ(2,3) -> CZ(3,4)
    # on qubit 3.
    first = {}
    for n in dag.nodes:
        if str(n).startswith("CZ"):
            first.setdefault(tuple(dag.nodes[n]["qubits"]), n)
    return {"wire": [(first[(0, 1)], first[(1, 2)]), (first[(2, 3)], first[(3, 4)])]}


@unittest.skipUnless(QISKIT, "qiskit required")
class GoldenWireCutTest(BaseTimedTestCase):
    def _plan(self, qc, observable, optimized=False, golden=True):
        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan

        return _CutPlan(
            qc, _chain_cuts(qc), observable, False, optimized, True, golden=golden
        )

    def test_real_chain_drops_the_y_terms(self):
        qc = _chain()
        plan = self._plan(qc, "ZZZZZZ")
        self.assertEqual([len(d) for d in plan.wire_digits], [6, 6])
        self.assertEqual(plan.gamma(), 9.0)  # 3 * 3, not 4 * 4
        self.assertEqual(plan.num_configurations, 36)
        # optimized: A sends (X, Z) = 2; B prepares 3 states x 2 settings; C 3.
        self.assertEqual(self._plan(qc, "ZZZZZZ", optimized=True).num_units(), 11)

    def test_exact_with_and_without_golden(self):
        import Qdislib as qd
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc

        for real in (True, False):
            qc = _chain(real)
            for observable in ("ZZZZZZ", "ZXIZZX", "ZYZZZZ"):
                for optimized in (False, True):
                    value = qd.circuit_cutting(
                        qc,
                        _chain_cuts(qc),
                        observable,
                        shots=None,
                        optimized=optimized,
                    )
                    self.assertAlmostEqual(value, _exact(qc, observable), places=9)
        # and the same value with golden cut points off
        original = mc._golden_wires
        mc._golden_wires = lambda graph, fragments: set()
        try:
            qc = _chain()
            value = qd.circuit_cutting(qc, _chain_cuts(qc), "ZZZZZZ", shots=None)
        finally:
            mc._golden_wires = original
        self.assertAlmostEqual(value, _exact(qc, "ZZZZZZ"), places=9)

    def test_a_complex_gate_or_a_y_observable_keeps_the_y_terms(self):
        # RX in block A: wire 1 is not golden, so block B receives |+i> states
        # and wire 2 is not either.
        self.assertEqual(
            [len(d) for d in self._plan(_chain(False), "ZZZZZZ").wire_digits], [8, 8]
        )
        # Y on qubit 1 at the end: its basis change is a complex gate in block B,
        # so only wire 1 (sent by the real block A) stays golden.
        self.assertEqual(
            [len(d) for d in self._plan(_chain(), "ZYZZZZ").wire_digits], [6, 8]
        )

    def test_two_wires_from_one_fragment_are_not_golden(self):
        """<O Y Y> of a real state need not vanish: a fragment that sends two
        wires keeps every term."""
        from qiskit import QuantumCircuit

        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan

        qc = QuantumCircuit(4)
        for q in range(4):
            qc.ry(0.3 + 0.4 * q, q)
        qc.cz(0, 1)
        qc.cz(1, 2)
        qc.cz(0, 3)
        cuts = {"wire": [("CZ_5", "CZ_6"), ("CZ_5", "CZ_7")]}
        plan = _CutPlan(qc, cuts, "ZZZZ", False, False, True)
        self.assertEqual([len(d) for d in plan.wire_digits], [8, 8])

    def test_shots_and_num_samples_stay_unbiased(self):
        import Qdislib as qd

        qc = _chain()
        ref = _exact(qc, "ZZZZZZ")
        sampled = [
            qd.circuit_cutting(
                qc, _chain_cuts(qc), "ZZZZZZ", shots=500, seed=s, num_samples=200
            )
            for s in range(20)
        ]
        # gamma = 9 and 36 configurations, all fitting in 200 draws: the mean of 20
        # runs is well inside 0.05.
        self.assertLess(abs(np.mean(sampled) - ref), 0.05)
        value = qd.circuit_cutting(
            qc, _chain_cuts(qc), "ZZZZZZ", shots=20000, seed=3, optimized=True
        )
        self.assertLess(abs(value - ref), 0.05)


@unittest.skipUnless(QISKIT, "qiskit required")
class CheaperHadamardWarningTest(BaseTimedTestCase):
    def _circuit(self, angle):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        for q in range(4):
            qc.h(q)
        qc.cx(0, 1)
        qc.rzz(angle, 1, 2)
        qc.cx(2, 3)
        return qc

    def test_warns_for_a_weak_rotation(self):
        import Qdislib as qd

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cost = qd.cutting_overhead(self._circuit(0.1), {"gate": ["RZZ_6"]})
        self.assertEqual(cost["gamma"], 9.0)
        messages = [str(w.message) for w in caught if w.category is UserWarning]
        self.assertTrue(any("Hadamard cut is cheaper" in m for m in messages))

    def test_silent_for_cz_and_for_the_hadamard_cut(self):
        import Qdislib as qd

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            qd.cutting_overhead(self._circuit(0.1), {"hadamard": ["RZZ_6"]})
            # pi/2 is CZ-like: gamma 3 both ways, nothing to point out
            qd.cutting_overhead(self._circuit(np.pi / 2), {"gate": ["RZZ_6"]})
            qd.cutting_overhead(self._circuit(0.1), {"gate": ["CX_5"]})
        self.assertFalse(
            [w for w in caught if "Hadamard cut is cheaper" in str(w.message)]
        )


@unittest.skipUnless(QISKIT and AER, "qiskit and qiskit-aer required")
class IQMJointCountsTest(BaseTimedTestCase):
    """The IQM backend gives the engine joint final/cut-point counts (here the
    device is stood in for by Aer)."""

    def test_gate_cut_reconstructs_on_iqm(self):
        from qiskit import QuantumCircuit
        from qiskit_aer import AerSimulator

        import Qdislib as qd
        from Qdislib.core.backend.base_backend import BackendConfig
        from Qdislib.core.backend.iqm_backend import IQMBackend

        device = object.__new__(IQMBackend)
        device.name = "IQM_BACKEND"
        device.config = BackendConfig(num_qubits=8)
        simulator = AerSimulator(seed_simulator=7)
        device._counts = lambda circuit, transpile, shots: (
            simulator.run(circuit, shots=shots).result().get_counts()
        )

        qc = QuantumCircuit(2)
        qc.ry(0.7, 0)
        qc.rx(1.1, 1)
        qc.cz(0, 1)
        qc.ry(0.9, 0)
        qc.rx(0.6, 1)
        qc.cz(0, 1)
        qc.ry(0.4, 0)
        for observable in ("ZZ", "IZ"):
            value = qd.gate_cutting(
                qc, ["CZ_3", "CZ_6"], observable, backend=device, shots=8000
            )
            self.assertLess(abs(value - _exact(qc, observable)), 0.08)


@unittest.skipUnless(QISKIT and AER, "qiskit and qiskit-aer required")
class BSCQpuJointCountsTest(BaseTimedTestCase):
    """The BSC QPU returns final counts only: a cut qubit that keeps evolving after
    its cut point has its measurement deferred, so the gate cut stays exact (the
    bridge is stood in for by an exact sampler)."""

    def test_non_terminal_cut_is_deferred(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd
        import Qdislib.core.backend.qpu_backends as qpu
        from Qdislib.core.backend.base_backend import BackendConfig

        widths = []

        def final_counts(dag, num_qubits):
            from Qdislib.core.backend.aer_backend import AerBackend

            widths.append(num_qubits)
            aer = AerBackend(BackendConfig(num_qubits=num_qubits))
            return aer.evaluate_samples(dag, num_qubits, 20000, seed=len(widths))

        device = object.__new__(qpu.BSCQpuBackend)
        device.name = "BSC_QPU_BACKEND"
        device.config = BackendConfig(num_qubits=8)
        device._final_counts = final_counts

        qc = QuantumCircuit(2)
        qc.ry(0.7, 0)
        qc.rx(1.1, 1)
        qc.cz(0, 1)
        qc.ry(0.9, 0)  # both cut qubits keep evolving after the cut
        qc.rx(0.6, 1)
        for observable in ("ZZ", "ZI"):
            value = qd.gate_cutting(
                qc, ["CZ_3"], observable, backend=device, shots=20000
            )
            self.assertLess(abs(value - _exact(qc, observable)), 0.05)
        # Each fragment is one qubit; its cut point, not terminal, adds one.
        self.assertEqual(max(widths), 2)


if __name__ == "__main__":
    unittest.main()
