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

"""Utils Circuit Unit tests."""

from tests import BaseTimedTestCase


class CircuitTest(BaseTimedTestCase):

    def test_optimal_cut_wire_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo

        qibo.set_backend("numpy")

        circuit = models.Circuit(4)
        circuit.add(gates.H(0))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1, 2))
        circuit.add(gates.RY(1, np.pi / 5))
        circuit.add(gates.RZ(2, np.pi / 3))
        circuit.add(gates.RZ(3, np.pi / 2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 3, gate_cut=False)

        reconstrucion = wire_cutting(circuit, cut)
        analytical = analytical_solution(circuit, "ZZZZ")

        if len(cut) == 1 and abs(reconstrucion - analytical) < 0.2:
            self.assertTrue(True)

        else:
            self.assertTrue(False)

    def test_optimal_cut_gate_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo

        qibo.set_backend("numpy")

        circuit = models.Circuit(4)
        circuit.add(gates.H(0))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1, 2))
        circuit.add(gates.RY(1, np.pi / 5))
        circuit.add(gates.RZ(2, np.pi / 3))
        circuit.add(gates.RZ(3, np.pi / 2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 2, wire_cut=False)

        reconstrucion = gate_cutting(circuit=circuit, cuts=cut)
        analytical = analytical_solution(circuit, "ZZZZ")

        if cut == ["CZ_5"] and abs(reconstrucion - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_zero_budgets_are_limits_not_none(self):
        """bsc-wdc/qdislib#15: max_cuts=0 allows no cut; max_qubits=0 is rejected."""
        from Qdislib.core.find_cut.find_cut import find_cut
        from qibo import models, gates

        def build():
            circuit = models.Circuit(3)
            for q in range(3):
                circuit.add(gates.RY(q, 0.3))
            circuit.add(gates.CZ(0, 1))
            circuit.add(gates.RY(1, 0.4))
            circuit.add(gates.CZ(1, 2))
            return circuit

        self.assertTrue(find_cut(build(), max_qubits=2, seed=0))
        for kwargs in ({}, {"gate_cut": False}, {"wire_cut": False}):
            self.assertEqual(
                find_cut(build(), max_qubits=2, max_cuts=0, seed=0, **kwargs), []
            )
        # Mixed mode returns an over-budget cut with a warning, but 0 is a hard stop.
        self.assertEqual(
            find_cut(build(), max_qubits=2, max_cuts=0, seed=0, mixed=True),
            {"gate": [], "wire": []},
        )
        # q1 and q2 are two hops apart here, so CZ(1, 2) is a SparseCut candidate.
        star = [(0, 1), (0, 2)]
        self.assertTrue(find_cut(build(), coupling_map=star))
        self.assertEqual(find_cut(build(), max_cuts=0, coupling_map=star), [])
        with self.assertRaises(ValueError):
            find_cut(build(), max_qubits=0)
        with self.assertRaises(ValueError):
            find_cut(build(), max_cuts=-1)

    def test_barriers_are_never_cut(self):
        """bsc-wdc/qdislib#16: a Qiskit barrier is not a cut point in any mode."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.core.find_cut.find_cut import find_cut
        from qiskit import QuantumCircuit

        def build():
            circuit = QuantumCircuit(3)
            circuit.h(0)
            circuit.barrier()
            circuit.ry(0.4, 1)
            circuit.cz(0, 1)
            circuit.cz(1, 2)
            circuit.ry(0.2, 2)
            return circuit

        wire = find_cut(build(), gate_cut=False, max_qubits=2, seed=0)
        gate = find_cut(build(), wire_cut=False, max_qubits=2, seed=0)
        mixed = find_cut(build(), max_qubits=2, seed=0, mixed=True)
        self.assertTrue(wire)
        self.assertNotIn("BARRIER", repr((wire, gate, mixed)))
        # The cut must be usable: the barrier's index is kept, so names line up.
        wire_cutting(build(), wire, observable="ZZZ", backend="qiskit", shots=1000)

    def test_wire_mode_respects_max_qubits_and_forbidden_qubits(self):
        """bsc-wdc/qdislib#10: wire cuts keep pieces within max_qubits and off
        forbidden qubits."""
        from Qdislib.core.cutting_algorithms.wire_cutting import _expand_wire_cuts
        from Qdislib.core.find_cut.find_cut import _wire_cut_fits, find_cut
        from Qdislib.utils.backends import convert_to_dag
        from qibo import models, gates

        chain = models.Circuit(6)
        for q in range(5):
            chain.add(gates.CZ(q, q + 1))
        cut = find_cut(chain, gate_cut=False, wire_cut=True, max_qubits=3, seed=0)
        self.assertTrue(cut)
        self.assertTrue(_wire_cut_fits(convert_to_dag(chain), cut, 3, None, None))

        layered = models.Circuit(4)
        for _ in range(2):
            for q in range(4):
                layered.add(gates.H(q))
            for q in range(3):
                layered.add(gates.CZ(q, q + 1))
        cut = find_cut(
            layered, gate_cut=False, wire_cut=True, forbidden_qubits=[1], seed=0
        )
        self.assertTrue(cut)
        cut_qubits = {q for _, _, q in _expand_wire_cuts(convert_to_dag(layered), cut)}
        self.assertNotIn(1, cut_qubits)
