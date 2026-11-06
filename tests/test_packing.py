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

"""Packing small circuits side by side on a large device (``Qdislib.pack``)."""

import importlib.util
import unittest
import warnings

from tests import BaseTimedTestCase

AER = importlib.util.find_spec("qiskit_aer") is not None


def _device():
    """A 57-qubit heavy-hex device with CZ gates and calibrated errors."""
    from qiskit.providers.fake_provider import GenericBackendV2
    from qiskit.transpiler import CouplingMap

    cmap = CouplingMap.from_heavy_hex(5)
    return GenericBackendV2(
        cmap.size(),
        basis_gates=["cz", "rz", "sx", "x", "id"],
        coupling_map=cmap,
        seed=11,
    )


def _run(circuits, shots, method="matrix_product_state", seed=5):
    """Noiseless counts of each (device-wide) circuit."""
    from qiskit_aer import AerSimulator

    sim = AerSimulator(method=method)
    result = sim.run(circuits, shots=shots, seed_simulator=seed).result()
    return [result.get_counts(i) for i in range(len(circuits))]


def _exact(circuit):
    """{bits: probability} keyed like Qiskit counts (qubit 0 rightmost)."""
    from qiskit.quantum_info import Statevector

    return Statevector(circuit).probabilities_dict()


def _target(couplers, num_qubits):
    """A device with a CX per coupler, ``{(a, b): error}``."""
    from qiskit.circuit import Parameter
    from qiskit.circuit.library import CXGate, Measure, RZGate, SXGate, XGate
    from qiskit.transpiler import InstructionProperties, Target

    target = Target(num_qubits=num_qubits)
    props = {}
    for (a, b), error in couplers.items():
        props[(a, b)] = props[(b, a)] = InstructionProperties(error=error)
    target.add_instruction(CXGate(), props)
    for gate in (RZGate(Parameter("t")), SXGate(), XGate(), Measure()):
        target.add_instruction(
            gate, {(q,): InstructionProperties(error=0.0) for q in range(num_qubits)}
        )
    return target


@unittest.skipUnless(AER, "qiskit-aer required")
class NoisyCouplerTest(BaseTimedTestCase):
    """A coupler over ``max_error`` carries no gate, but is still there."""

    def _bell(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        return qc

    def test_it_still_counts_for_the_distance_between_groups(self):
        import Qdislib as qd

        line = _target({(0, 1): 0.001, (1, 2): 0.5, (2, 3): 0.001, (3, 4): 0.001}, 5)
        without_buffer = qd.pack([self._bell()], line, buffer=0)
        self.assertEqual(without_buffer.regions, [(0, 1), (2, 3)])
        # (0, 1) and (2, 3) are neighbours through it: with an idle qubit between them
        self.assertEqual(
            qd.pack([self._bell()], line, buffer=1).regions, [(0, 1), (3, 4)]
        )

    def test_no_gate_runs_on_it_inside_a_group(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        triangle = _target({(0, 1): 0.001, (1, 2): 0.001, (0, 2): 0.5}, 3)
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 2)
        qc.cx(0, 1)
        qc.cx(1, 2)
        packed = qd.pack([qc], triangle, seeds=4)
        used = {
            frozenset(packed.circuits[0].find_bit(q).index for q in i.qubits)
            for i in packed.circuits[0].data
            if i.operation.num_qubits == 2
        }
        self.assertNotIn(frozenset((0, 2)), used)
        # given by hand, a group is used as it is
        given = qd.pack([self._bell()], triangle, regions=[(0, 2)], seeds=1)
        self.assertEqual(given.regions, [(0, 2)])


@unittest.skipUnless(AER, "qiskit-aer required")
class PackingTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        warnings.simplefilter("ignore")
        import Qdislib as qd

        self.qd, self.device = qd, _device()

    def test_basis_states_come_back_exactly(self):
        from qiskit import QuantumCircuit

        circuits, want = [], []
        for k in range(30):
            n = 2 + k % 3
            qc = QuantumCircuit(n)
            bits = [(k >> j) & 1 for j in range(n)]
            for j, b in enumerate(bits):
                if b:
                    qc.x(j)
            circuits.append(qc)
            want.append("".join(str(b) for b in reversed(bits)))
        packed = self.qd.pack(circuits, self.device)
        self.assertGreater(len(packed.regions), 1)
        self.assertLess(len(packed.circuits), len(circuits))
        counts = packed.unpack(_run(packed.circuits, 100, method="stabilizer"))
        for got, bits in zip(counts, want):
            self.assertEqual(got, {bits: 100})

    def test_groups_are_disjoint_and_apart(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        qc.h(0)
        for j in range(3):
            qc.cx(j, j + 1)
        packed = self.qd.pack([qc] * 12, self.device, buffer=1)
        coupled = {frozenset(e) for e in self.device.coupling_map.get_edges()}
        groups = [set(g) for g in packed.regions]
        for i, a in enumerate(groups):
            for b in groups[i + 1 :]:
                self.assertFalse(a & b)
                self.assertFalse(
                    any(frozenset((x, y)) in coupled for x in a for y in b)
                )
        for big in packed.circuits:  # only device couplers, on its own groups
            for inst in big.data:
                if inst.operation.num_qubits == 2:
                    pair = frozenset(big.find_bit(q).index for q in inst.qubits)
                    self.assertIn(pair, coupled)

    def test_random_circuits_match_their_distribution(self):
        from qiskit.circuit.random import random_circuit

        circuits = [
            random_circuit(3 + k % 2, 4, max_operands=2, seed=k) for k in range(10)
        ]
        packed = self.qd.pack(circuits, self.device)
        counts = packed.unpack(_run(packed.circuits, 20000))
        for circuit, got in zip(circuits, counts):
            want = _exact(circuit)
            for key in set(want) | set(got):
                self.assertAlmostEqual(
                    got.get(key, 0) / 20000, want.get(key, 0.0), delta=0.02
                )

    def test_cut_circuit_reconstructs_from_packed_subcircuits(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        qc.ry(0.7, 0)
        qc.rx(1.2, 1)
        qc.cx(0, 2)
        qc.cx(1, 3)
        qc.ry(0.4, 2)
        qc.rz(0.9, 3)
        qc.h(3)
        subs = self.qd.circuit_cutting_sampling_subcircuits(
            qc, {"hadamard": ["CX_3", "CX_4"]}
        )
        circuits = subs.circuits("qiskit", defer_measurements=True)
        for c in circuits:
            c.measure_all()
        packed = self.qd.pack(circuits, self.device)
        counts = packed.unpack(_run(packed.circuits, 20000))
        got = subs.reconstruct_counts(counts, bit_order="q0_last")
        exact = {k[::-1]: v for k, v in _exact(qc).items()}  # qubit-0-first keys
        for key in set(exact) | set(got):
            self.assertAlmostEqual(got.get(key, 0.0), exact.get(key, 0.0), delta=0.04)

    def test_several_registers_keep_their_keys(self):
        from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

        q = QuantumRegister(3)
        a, b = ClassicalRegister(1, "a"), ClassicalRegister(2, "b")
        qc = QuantumCircuit(q, a, b)
        qc.x(0)
        qc.x(2)
        qc.measure(0, a[0])
        qc.measure(1, b[0])
        qc.measure(2, b[1])
        alone = _run([qc], 50, method="stabilizer")[0]
        packed = self.qd.pack([qc, qc], self.device)
        got = packed.unpack(_run(packed.circuits, 50, method="stabilizer"))
        self.assertEqual(got, [alone, alone])

    def test_given_regions(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.x(1)
        first = self.qd.pack([qc] * 4, self.device)
        group = first.regions[2]
        alone = self.qd.pack([qc], self.device, regions=[group])
        self.assertEqual(alone.slots[0][2], group)
        used = {
            alone.circuits[0].find_bit(q).index
            for inst in alone.circuits[0].data
            for q in inst.qubits
        }
        self.assertTrue(used <= set(group))
        counts = alone.unpack(_run(alone.circuits, 10, method="stabilizer"))
        self.assertEqual(counts, [{"10": 10}])
        with self.assertRaises(ValueError):
            self.qd.pack([QuantumCircuit(3)], self.device, regions=[group])

    def test_bad_input(self):
        from qiskit import QuantumCircuit

        with self.assertRaises(ValueError):
            self.qd.pack([], self.device)
        with self.assertRaises(ValueError):
            self.qd.pack([QuantumCircuit(3)], self.device, size=2)
        packed = self.qd.pack([QuantumCircuit(2)], self.device)
        with self.assertRaises(ValueError):
            packed.unpack([])
        with self.assertRaises(ValueError):  # classical bits, no measurement
            self.qd.pack([QuantumCircuit(2, 2)], self.device)


if __name__ == "__main__":
    unittest.main()
