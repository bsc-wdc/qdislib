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

"""Joint measurement of the cut wires a fragment sends (wire_cutting(joint=True)).

k wires leaving one fragment are measured with 2**k + 1 circuits instead of 3**k,
each reading a commuting group of Paulis, every Pauli exactly once.
"""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


class JointSettingsTest(BaseTimedTestCase):
    def test_each_pauli_is_read_once(self):
        import itertools

        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _joint_measurement_settings,
        )

        for k in (1, 2, 3, 4):
            settings = _joint_measurement_settings(k)
            self.assertEqual(len(settings), 2**k + 1)
            read = [
                p for _, readout in settings for p, _, _ in readout if set(p) != {"I"}
            ]
            expected = {
                p for p in itertools.product("IXYZ", repeat=k) if set(p) != {"I"}
            }
            self.assertEqual(sorted(read), sorted(expected))  # each once, none missing


class _Builder:
    """A Qiskit circuit with random single-qubit layers, tracking gate names."""

    def __init__(self, num_qubits, seed):
        from qiskit import QuantumCircuit

        self.qc = QuantumCircuit(num_qubits)
        self.rng = np.random.default_rng(seed)
        self.last = {}

    def layer(self, *qubits):
        for q in qubits:
            self.qc.u(*self.rng.uniform(0, 2 * np.pi, 3), q)
            self.last[q] = f"U_{len(self.qc.data)}"

    def cz(self, a, b):
        self.qc.cz(a, b)
        name = f"CZ_{len(self.qc.data)}"
        self.last[a] = self.last[b] = name
        return name

    def cut_into(self, qubit, other):
        """A wire cut on ``qubit`` right before a new CZ(qubit, other)."""
        source = self.last[qubit]
        return (source, self.cz(qubit, other))


def _exact(circuit, observable):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(observable[::-1])).real
    )


@unittest.skipUnless(QISKIT, "qiskit not installed")
class JointWireMeasurementTest(BaseTimedTestCase):
    def _check(self, circuit, cuts, optimized_count, joint_count):
        import Qdislib as qd

        obs = "Z" * circuit.num_qubits
        for joint, count in ((False, optimized_count), (True, joint_count)):
            subs = qd.wire_cutting_subcircuits(
                circuit,
                cuts,
                observables=obs,
                software="qiskit",
                optimized=True,
                joint=joint,
            )
            self.assertEqual(len(subs), count)
        value = qd.wire_cutting(
            circuit,
            cuts,
            observables=obs,
            software="qiskit",
            shots=40000,
            optimized=True,
        )
        self.assertLess(abs(value - _exact(circuit, obs)), 0.05)

    def test_two_wires_to_one_fragment(self):
        b = _Builder(4, seed=4)
        b.layer(0, 1, 2, 3)
        b.cz(0, 1), b.cz(1, 2), b.cz(0, 2)
        b.layer(0, 1, 2)
        cuts = [b.cut_into(1, 3), b.cut_into(2, 3)]
        b.layer(1, 2, 3)
        self._check(b.qc, cuts, optimized_count=9 + 16, joint_count=5 + 16)

    def test_two_wires_to_two_fragments(self):
        b = _Builder(5, seed=5)
        b.layer(0, 1, 2, 3, 4)
        b.cz(0, 1), b.cz(1, 2), b.cz(0, 2)
        b.layer(0, 1, 2)
        cuts = [b.cut_into(1, 3), b.cut_into(2, 4)]
        b.layer(1, 2, 3, 4)
        self._check(b.qc, cuts, optimized_count=9 + 4 + 4, joint_count=5 + 4 + 4)

    def test_mixed_fragment_sends_two_wires(self):
        b = _Builder(5, seed=6)
        b.layer(0, 1, 2, 3, 4)
        b.cz(0, 1)
        b.layer(0, 1)
        cuts = [b.cut_into(1, 2)]
        b.cz(2, 3), b.cz(1, 3)
        b.layer(1, 2, 3)
        cuts += [b.cut_into(2, 4), b.cut_into(3, 4)]
        b.layer(2, 3, 4)
        # A: 3 measure circuits; middle (receives 1, sends 2): 9*4 -> 5*4; C: 16.
        self._check(b.qc, cuts, optimized_count=3 + 36 + 16, joint_count=3 + 20 + 16)

    def test_sampling_distribution(self):
        """The output-distribution path measures jointly too (the default)."""
        import Qdislib as qd
        from qiskit.quantum_info import Statevector

        b = _Builder(4, seed=4)
        b.layer(0, 1, 2, 3)
        b.cz(0, 1), b.cz(1, 2), b.cz(0, 2)
        b.layer(0, 1, 2)
        cuts = [b.cut_into(1, 3), b.cut_into(2, 3)]
        b.layer(1, 2, 3)
        exact = {
            key[::-1]: p for key, p in Statevector(b.qc).probabilities_dict().items()
        }
        for joint in (False, True):
            dist = qd.wire_cutting_sampling(
                b.qc,
                cuts,
                shots=30000,
                software="qiskit",
                optimized=True,
                joint=joint,
                seed=3,
            )
            tvd = 0.5 * sum(
                abs(dist.get(k, 0) - exact.get(k, 0)) for k in set(dist) | set(exact)
            )
            self.assertLess(tvd, 0.03, joint)


if __name__ == "__main__":
    unittest.main()
