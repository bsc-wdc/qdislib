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


# -*- coding: utf-8 -*-

"""Circuits that need no cut at all.

Every strategy :func:`~Qdislib.find_cut` runs searches for the cheapest *cut*.
None asked whether one was needed, so a circuit whose qubit interaction graph is
already disconnected got a gate cut anyway -- ``gamma = 3``, nine times the
shots, for a split that was free.

The case that found it: Bernstein-Vazirani only entangles a qubit where its
secret bit is 1, so a 24-qubit instance leaves twelve qubits untouched by any
two-qubit gate. It splits 12/12 at zero cost, and ``find_cut`` was returning a
gate cut with lopsided 1q and 11q fragments.
"""

import unittest

import numpy as np

from tests import BaseTimedTestCase


def two_blocks():
    """Two independent two-qubit blocks sharing no gate."""
    from qiskit import QuantumCircuit

    circuit = QuantumCircuit(4)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.h(2)
    circuit.cx(2, 3)
    circuit.rx(0.3, 1)
    circuit.rx(0.4, 3)
    return circuit


def a_ring(width=6):
    """Genuinely connected: every qubit reachable from every other."""
    from qiskit import QuantumCircuit

    circuit = QuantumCircuit(width)
    for qubit in range(width - 1):
        circuit.cx(qubit, qubit + 1)
    circuit.cx(0, width - 1)
    return circuit


class SeparableComponentsTest(BaseTimedTestCase):

    def test_two_blocks_separate_for_free(self):
        from Qdislib.core.find_cut.find_cut import separable_components

        groups = separable_components(two_blocks(), max_qubits=2)
        self.assertEqual(sorted(len(g) for g in groups), [2, 2])

    def test_a_connected_circuit_does_not(self):
        from Qdislib.core.find_cut.find_cut import separable_components

        self.assertIsNone(separable_components(a_ring(6), max_qubits=3))

    def test_a_component_wider_than_the_budget_still_needs_cutting(self):
        """Disconnected is not enough; each piece must also fit."""
        from Qdislib.core.find_cut.find_cut import separable_components

        self.assertIsNone(separable_components(two_blocks(), max_qubits=1))

    def test_components_are_packed_not_merely_listed(self):
        """Four independent pairs into a 4-qubit budget is two fragments, not four."""
        from qiskit import QuantumCircuit

        from Qdislib.core.find_cut.find_cut import separable_components

        circuit = QuantumCircuit(8)
        for pair in range(4):
            circuit.h(2 * pair)
            circuit.cx(2 * pair, 2 * pair + 1)
        groups = separable_components(circuit, max_qubits=4)
        self.assertEqual(sorted(len(g) for g in groups), [4, 4])

    def test_no_budget_means_any_split_will_do(self):
        from Qdislib.core.find_cut.find_cut import separable_components

        groups = separable_components(two_blocks(), max_qubits=None)
        self.assertEqual(len(groups), 2)

    def test_a_circuit_with_no_two_qubit_gates_is_trivially_separable(self):
        from qiskit import QuantumCircuit

        from Qdislib.core.find_cut.find_cut import separable_components

        circuit = QuantumCircuit(4)
        for qubit in range(4):
            circuit.h(qubit)
        self.assertIsNotNone(separable_components(circuit, max_qubits=2))


class FindCutReturnsNoneNeededTest(BaseTimedTestCase):

    def test_it_returns_no_cuts_when_none_are_needed(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        self.assertEqual(find_cut(two_blocks(), max_qubits=2), [])

    def test_it_still_cuts_a_connected_circuit(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        self.assertTrue(find_cut(a_ring(6), max_qubits=3))

    def test_the_ranked_mode_reports_the_free_split_too(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        ranked = find_cut(two_blocks(), max_qubits=2, top_n=3)
        self.assertEqual(ranked, [{"cut": [], "kind": "none", "score": 0.0}])

    def test_the_mixed_mode_does_as_well(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        self.assertEqual(
            find_cut(two_blocks(), max_qubits=2, mixed=True),
            {"gate": [], "hadamard": [], "wire": []},
        )


class ZeroCutReconstructionTest(BaseTimedTestCase):
    """Cutting nothing has to give the right answer, not merely run."""

    def _exact(self, circuit, observable):
        from qiskit.quantum_info import DensityMatrix, Pauli

        state = DensityMatrix.from_instruction(circuit).data
        # Qdislib reads observables with qubit 0 leftmost; qiskit.Pauli is the
        # other way round, so the reference string is reversed.
        return float(np.real(np.trace(Pauli(observable[::-1]).to_matrix() @ state)))

    def test_it_reconstructs_the_uncut_value(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        circuit = two_blocks()
        for observable in ("ZZZZ", "ZZII", "IIZZ"):
            with self.subTest(observable=observable):
                self.assertAlmostEqual(
                    gate_cutting(circuit, [], observable=observable, shots=200_000),
                    self._exact(circuit, observable),
                    delta=0.03,
                )

    def test_isolated_qubits_keep_their_share_of_the_observable(self):
        """A component with only single-qubit gates is still a fragment.

        Dropping it would silently ignore whatever the observable puts there.
        """
        from qiskit import QuantumCircuit

        from Qdislib.core.cutting_algorithms.gate_cutting import (
            gate_cutting,
            gate_cutting_subcircuits,
        )

        circuit = QuantumCircuit(4)
        circuit.h(0)
        circuit.cx(0, 1)
        circuit.rx(0.9, 2)
        circuit.ry(1.2, 3)
        pieces = gate_cutting_subcircuits(
            circuit, [], observable="ZZZZ", software="qiskit"
        )
        self.assertEqual(sorted(p[0].num_qubits for p in pieces), [1, 1, 2])
        self.assertAlmostEqual(
            gate_cutting(circuit, [], observable="ZZZZ", shots=200_000),
            self._exact(circuit, "ZZZZ"),
            delta=0.03,
        )

    def test_the_overhead_is_one(self):
        """No cut, no quasi-probability, no sampling penalty."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            gate_cutting_subcircuits,
        )

        pieces = gate_cutting_subcircuits(two_blocks(), [], software="qiskit")
        self.assertEqual(pieces.number_cuts, 0)


if __name__ == "__main__":
    unittest.main()
