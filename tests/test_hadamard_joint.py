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

"""Joint Hadamard cut: the gates crossing one boundary share one ancilla per side.

Exact (statevector) checks against the uncut value, as in
``docs/design/joint_hadamard_cut.md``: k gates across one boundary at different
times, general (KAK) gates, gates whose control switches sides, and three fragments.
"""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None


def _layer(circuit, rng):
    from qibo import gates

    for q in range(circuit.nqubits):
        circuit.add(gates.U3(q, *rng.uniform(0, 2 * np.pi, 3)))


def _last(circuit):
    return f"{circuit.queue[-1].__class__.__name__.upper()}_{len(circuit.queue)}"


def _one_boundary(k, gate="CZ", seed=7):
    """Sides A = 0..k-1 and B = k..2k-1, each linked by uncut CZs; cut gate m acts
    on (m, k+m) at its own time step, with random layers in between."""
    from qibo import gates, models

    rng = np.random.default_rng(seed)
    circuit = models.Circuit(2 * k)
    cut = []
    for m in range(k + 1):
        _layer(circuit, rng)
        for q in range(k - 1):
            circuit.add(gates.CZ(q, q + 1))
            circuit.add(gates.CZ(k + q, k + q + 1))
        if m < k:
            circuit.add(getattr(gates, gate)(m, k + m))
            cut.append(_last(circuit))
    return circuit, cut


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class JointHadamardCutTest(BaseTimedTestCase):
    def _exact(self, circuit, observable=None):
        import Qdislib as qd

        return qd.analytical_solution(circuit, observable or "Z" * circuit.nqubits)

    def _joint(self, circuit, cut, observable=None):
        from Qdislib.core.cutting_algorithms.hadamard_cut import gate_cutting_hadamard

        return gate_cutting_hadamard(
            circuit, cut, observable or "Z" * circuit.nqubits, joint=True, shots=None
        )

    def test_k_gates_across_one_boundary(self):
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuits,
            _hadamard_cut_overhead,
        )

        for gate in ("CZ", "CNOT"):
            for k in (1, 2, 3):
                circuit, cut = _one_boundary(k, gate)
                with self.subTest(gate=gate, k=k):
                    self.assertAlmostEqual(
                        self._joint(circuit, cut), self._exact(circuit), places=10
                    )
                    import Qdislib as qd

                    subs = gate_cutting_hadamard_subcircuits(
                        circuit, cut, "Z" * circuit.nqubits, joint=True
                    )
                    cost = qd.cutting_overhead(circuit, {"hadamard": cut}, joint=True)
                    self.assertEqual(cost["configurations"], 2**k * (2 ** (k + 1) - 1))
                    # one ancilla per side, however many gates cross
                    self.assertEqual({s.num_qubits for s in subs}, {k + 1})
                    phi, groups = _hadamard_cut_overhead(circuit, cut, joint=True)
                    self.assertAlmostEqual(phi, 2 ** (k + 1) - 1)
                    self.assertEqual(len(groups), 1)
                    self.assertAlmostEqual(
                        _hadamard_cut_overhead(circuit, cut)[0], 3**k
                    )

    def test_one_gate_is_todays_cut(self):
        """A group of one gate builds exactly today's fragments, setting by setting."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuits,
        )

        circuit, cut = _one_boundary(1)
        separate = gate_cutting_hadamard_subcircuits(circuit, cut, "ZZ", joint=False)
        joint = gate_cutting_hadamard_subcircuits(circuit, cut, "ZZ", joint=True)
        self.assertEqual(len(separate), len(joint))
        for a, b in zip(separate, joint):
            self.assertEqual(a.observable, b.observable)
            self.assertEqual(a.circuit.to_qasm(), b.circuit.to_qasm())

    def test_general_gates_in_one_group(self):
        from qibo import gates, models

        rng = np.random.default_rng(3)
        circuit = models.Circuit(4)
        _layer(circuit, rng)
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.iSWAP(1, 2))
        first = _last(circuit)
        _layer(circuit, rng)
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.fSim(0, 3, 0.7, 0.4))
        second = _last(circuit)
        _layer(circuit, rng)
        self.assertAlmostEqual(
            self._joint(circuit, [first, second]), self._exact(circuit), places=10
        )

    def test_control_switches_sides(self):
        from qibo import gates, models

        rng = np.random.default_rng(5)
        circuit = models.Circuit(4)
        _layer(circuit, rng)
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.CZ(1, 2))  # control on side A
        first = _last(circuit)
        _layer(circuit, rng)
        circuit.add(gates.CNOT(3, 0))  # control on side B
        second = _last(circuit)
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 3))
        _layer(circuit, rng)
        self.assertAlmostEqual(
            self._joint(circuit, [first, second]), self._exact(circuit), places=10
        )

    def test_three_fragments_two_boundaries(self):
        from qibo import gates, models

        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuits,
            _hadamard_cut_overhead,
        )

        rng = np.random.default_rng(9)
        circuit = models.Circuit(6)
        cut = []
        for step in range(2):
            _layer(circuit, rng)
            for a, b in ((0, 1), (2, 3), (4, 5)):
                circuit.add(gates.CZ(a, b))
            first, second = ((1, 2), (3, 4)) if step == 0 else ((0, 3), (2, 5))
            circuit.add(gates.CZ(*first))
            cut.append(_last(circuit))
            circuit.add(gates.CNOT(*second))
            cut.append(_last(circuit))
        _layer(circuit, rng)
        for observable in ("ZZZZZZ", "ZIIZIZ"):
            with self.subTest(observable=observable):
                self.assertAlmostEqual(
                    self._joint(circuit, cut, observable),
                    self._exact(circuit, observable),
                    places=10,
                )
        phi, groups = _hadamard_cut_overhead(circuit, cut, joint=True)
        self.assertEqual([len(gates_) for gates_, _ in groups], [2, 2])
        self.assertAlmostEqual(phi, 7 * 7)
        # the middle fragment holds one ancilla per boundary it touches
        subs = gate_cutting_hadamard_subcircuits(circuit, cut, "Z" * 6, joint=True)
        self.assertEqual(sorted({s[0].nqubits for s in subs}), [3, 4])

    def test_x_and_y_observables(self):
        """X/Y are rotated into Z before cutting, on every Hadamard path."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard,
            gate_cutting_hadamard_subcircuits,
        )

        circuit, cut = _one_boundary(2, seed=5)
        observables = ("XZYI", "YYXX", [(0.5, "XXZZ"), (-1.2, "ZIYI"), (0.3, "IIII")])
        for observable in observables:
            exact = self._exact(circuit, observable)
            for joint in (False, True):
                with self.subTest(observable=observable, joint=joint):
                    value = gate_cutting_hadamard(
                        circuit, cut, observable, joint=joint, shots=None
                    )
                    self.assertAlmostEqual(value, exact, places=10)
                    subs = gate_cutting_hadamard_subcircuits(
                        circuit, cut, observable, joint=joint
                    )
                    import Qdislib as qd

                    values = [qd.expectation_value(s, shots=None) for s in subs]
                    self.assertAlmostEqual(subs.reconstruct(values), exact, places=10)

    def test_sampling_is_unbiased(self):
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard,
        )

        circuit, cut = _one_boundary(2, seed=11)
        values = [
            gate_cutting_hadamard(
                circuit, cut, "ZZZZ", num_samples=3000, shots=None, seed=s, joint=True
            )
            for s in range(4)
        ]
        self.assertLess(abs(np.mean(values) - self._exact(circuit)), 0.1)


def _exact_distribution(circuit):
    state = np.asarray(circuit().state())
    n = circuit.nqubits
    return {
        "".join(str((i >> (n - 1 - q)) & 1) for q in range(n)): float(abs(a) ** 2)
        for i, a in enumerate(state)
    }


def _tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class HadamardSamplingTest(BaseTimedTestCase):
    """gate_cutting_hadamard_sampling: the output distribution, like
    gate_cutting_sampling / wire_cutting_sampling."""

    def test_exact_distribution(self):
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
        )

        for k in (1, 2):
            circuit, cut = _one_boundary(k, gate="CNOT", seed=3)
            exact = _exact_distribution(circuit)
            for joint in (False, True):
                dist = gate_cutting_hadamard_sampling(
                    circuit, cut, joint=joint, shots=None
                )
                self.assertLess(_tvd(dist, exact), 1e-9, (k, joint))

    def test_sampled_and_three_step(self):
        import Qdislib as qd

        circuit, cut = _one_boundary(2, seed=8)
        exact = _exact_distribution(circuit)
        dist = qd.gate_cutting_hadamard_sampling(
            circuit, cut, backend="qiskit", shots=20000, joint=True, seed=1
        )
        self.assertLess(_tvd(dist, exact), 0.03)
        subs = qd.gate_cutting_hadamard_sampling_subcircuits(circuit, cut, joint=True)
        sampled = [qd.sample_subcircuit(s, shots=None) for s in subs]
        self.assertLess(_tvd(subs.reconstruct(sampled), exact), 1e-9)

    def test_old_expectation_call_is_refused(self):
        """Up to 1.0.1 this name estimated an expectation value."""
        import Qdislib as qd

        circuit, cut = _one_boundary(1)
        with self.assertRaises(TypeError):
            qd.gate_cutting_hadamard_sampling(
                circuit, cut, "ZZ", num_samples=100, shots=None
            )
        with self.assertRaises(TypeError):
            qd.gate_cutting_hadamard_sampling(circuit, cut, observable="ZZ", shots=None)


@unittest.skipUnless(
    QIBO_AVAILABLE and importlib.util.find_spec("qiskit") is not None,
    "qibo and qiskit needed",
)
class QiskitCircuitTest(BaseTimedTestCase):
    """A Qiskit circuit keeps Qiskit's gate names (u, p, cx, cp) in its DAG; the
    fragments must still build."""

    def test_qiskit_gate_names(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import SparsePauliOp, Statevector

        rng = np.random.default_rng(5)
        qc = QuantumCircuit(4)
        for q in range(4):
            qc.u(*rng.uniform(0, 2 * np.pi, 3), q)
        qc.cx(0, 1)
        qc.cp(0.4, 2, 3)
        qc.cx(1, 2)
        qc.rzz(0.7, 0, 3)
        for q in range(4):
            qc.p(rng.uniform(0, 2 * np.pi), q)
            qc.u(*rng.uniform(0, 2 * np.pi, 3), q)
        qc.cz(0, 1)
        qc.cx(2, 3)
        exact = Statevector(qc).expectation_value(SparsePauliOp("ZZZZ")).real
        for joint in (True, False):
            value = qd.gate_cutting_hadamard(
                qc, ["CX_7", "RZZ_8"], joint=joint, shots=None
            )
            self.assertAlmostEqual(value, exact, places=9)


if __name__ == "__main__":
    unittest.main()
