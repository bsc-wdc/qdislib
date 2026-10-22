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

"""circuit_cutting: any combination of measure-and-prepare gate, Hadamard gate and
wire cuts, as an expectation value or a distribution, all-in-one or three-step.

Exact (``shots=None``) checks against the statevector, so a bookkeeping slip in the
combined decomposition shows up as a wrong number, not as noise.
"""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


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

    def gate(self, name, a, b):
        getattr(self.qc, name)(a, b)
        node = f"{name.upper()}_{len(self.qc.data)}"
        self.last[a] = self.last[b] = node
        return node

    def wire(self, qubit, other, name="cz"):
        """A wire cut on ``qubit`` right before a new gate with ``other``."""
        source = self.last[qubit]
        return (source, self.gate(name, qubit, other))


def _all_kinds():
    """A = {0, 1} sends two wires to C = {3}; a Hadamard CZ joins A and B = {2};
    a measure-and-prepare CX joins B and D = {4}."""
    b = _Builder(5, seed=1)
    b.layer(0, 1, 2, 3, 4)
    b.gate("cz", 0, 1)
    had = b.gate("cz", 1, 2)
    b.layer(0, 1, 2)
    mp = b.gate("cx", 2, 4)
    b.layer(0, 1, 2, 4)
    wires = [b.wire(0, 3), b.wire(1, 3)]
    b.layer(0, 1, 3)
    return b.qc, {"hadamard": [had], "gate": [mp], "wire": wires}


def _chain():
    """Two wire cuts in series (the middle fragment receives one and sends one),
    plus a gate cut on the middle fragment."""
    b = _Builder(4, seed=2)
    b.layer(0, 1, 2, 3)
    b.gate("cz", 0, 1)
    b.layer(0, 1)
    wires = [b.wire(1, 2)]
    b.layer(1, 2)
    gate = b.gate("cz", 2, 3)
    b.layer(2, 3)
    wires.append(b.wire(2, 3, "cx"))
    b.layer(2, 3)
    return b.qc, {"gate": [gate], "wire": wires}


def _hadamard_group():
    """Two Hadamard-cut gates across one boundary (one joint group) and a wire."""
    b = _Builder(5, seed=3)
    b.layer(0, 1, 2, 3, 4)
    b.gate("cz", 0, 1), b.gate("cz", 2, 3)
    had = [b.gate("cz", 1, 2)]
    b.layer(0, 1, 2, 3)
    b.gate("cz", 0, 1), b.gate("cz", 2, 3)
    had.append(b.gate("cx", 0, 3))
    b.layer(0, 1, 2, 3)
    wires = [b.wire(3, 4)]
    b.layer(3, 4)
    return b.qc, {"hadamard": had, "wire": wires}


def _exact(circuit, observable):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(observable[::-1])).real
    )


def _exact_distribution(circuit):
    from qiskit.quantum_info import Statevector

    return {k[::-1]: p for k, p in Statevector(circuit).probabilities_dict().items()}


def _tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


@unittest.skipUnless(QISKIT, "qiskit not installed")
class CircuitCuttingExactTest(BaseTimedTestCase):
    MODES = (
        {"optimized": False},
        {"optimized": True},
        {"optimized": True, "joint": True},
        {"optimized": True, "joint": False},
    )

    def _check(self, circuit, cuts):
        import Qdislib as qd

        n = circuit.num_qubits
        exact_dist = _exact_distribution(circuit)
        for mode in self.MODES:
            for observable in ("Z" * n, "XY" + "Z" * (n - 2)):
                with self.subTest(mode=mode, observable=observable):
                    value = qd.circuit_cutting(
                        circuit, cuts, observable, shots=None, **mode
                    )
                    self.assertAlmostEqual(
                        value, _exact(circuit, observable), places=10
                    )
            with self.subTest(mode=mode, distribution=True):
                dist = qd.circuit_cutting_sampling(circuit, cuts, shots=None, **mode)
                self.assertLess(_tvd(dist, exact_dist), 1e-9)

    def test_every_kind_at_once(self):
        self._check(*_all_kinds())

    def test_wires_in_series_with_a_gate_cut(self):
        self._check(*_chain())

    def test_joint_hadamard_group_with_a_wire(self):
        self._check(*_hadamard_group())

    def test_three_step_and_pauli_sum(self):
        import Qdislib as qd

        circuit, cuts = _all_kinds()
        observable = [(0.5, "ZZZZZ"), (-1.5, "XIYZI")]
        exact = 0.5 * _exact(circuit, "ZZZZZ") - 1.5 * _exact(circuit, "XIYZI")
        subs = qd.circuit_cutting_subcircuits(circuit, cuts, observable, optimized=True)
        values = [qd.expectation_value(s, shots=None) for s in subs]
        self.assertAlmostEqual(subs.reconstruct(values), exact, places=10)
        self.assertAlmostEqual(
            qd.circuit_cutting(circuit, cuts, observable, shots=None, optimized=True),
            exact,
            places=10,
        )

        subs = qd.circuit_cutting_sampling_subcircuits(circuit, cuts, optimized=True)
        sampled = [qd.sample_subcircuit(s, shots=None) for s in subs]
        self.assertLess(
            _tvd(subs.reconstruct(sampled), _exact_distribution(circuit)), 1e-9
        )

    def test_with_shots(self):
        import Qdislib as qd

        circuit, cuts = _chain()
        n = circuit.num_qubits
        value = qd.circuit_cutting(
            circuit,
            cuts,
            "Z" * n,
            software="qiskit",
            shots=20000,
            optimized=True,
            seed=1,
        )
        self.assertLess(abs(value - _exact(circuit, "Z" * n)), 0.05)
        dist = qd.circuit_cutting_sampling(
            circuit, cuts, software="qiskit", shots=20000, optimized=True, seed=2
        )
        self.assertLess(_tvd(dist, _exact_distribution(circuit)), 0.05)


def _boundary(k, kind):
    """Fragments A = 0..k-1 and B = k..2k-1 joined by k cuts of one kind."""
    b = _Builder(2 * k, seed=k)
    b.layer(*range(2 * k))
    for q in range(k - 1):
        b.gate("cz", q, q + 1), b.gate("cz", k + q, k + q + 1)
    if kind == "wire":
        return b.qc, {"wire": [b.wire(m, k + m) for m in range(k)]}
    cut = [b.gate("cz", m, k + m) for m in range(k)]
    for q in range(k - 1):
        b.gate("cz", q, q + 1), b.gate("cz", k + q, k + q + 1)
    return b.qc, {kind: cut}


@unittest.skipUnless(QISKIT, "qiskit not installed")
class CircuitCuttingOverheadTest(BaseTimedTestCase):
    """How the number of circuits scales with k cuts across one boundary."""

    def _count(self, kind, k, **kw):
        import Qdislib as qd

        return qd.circuit_cutting_overhead(*_boundary(k, kind), **kw)

    def test_counts(self):
        for k in (1, 2, 3):
            with self.subTest(k=k):
                self.assertEqual(self._count("gate", k)["subcircuits"], 2 * 6**k)
                # each side of a measure-and-prepare cut has 5 distinct circuits
                self.assertEqual(
                    self._count("gate", k, optimized=True)["subcircuits"], 2 * 5**k
                )
                joint = self._count("hadamard", k, optimized=True, joint=True)
                self.assertEqual(joint["subcircuits"], 2 * 2**k * (2 ** (k + 1) - 1))
                self.assertAlmostEqual(joint["gamma"], 2 ** (k + 1) - 1)
                self.assertEqual(self._count("wire", k)["subcircuits"], 2 * 8**k)
                self.assertEqual(
                    self._count("wire", k, optimized=True, joint=False)["subcircuits"],
                    3**k + 4**k,
                )
                wire = self._count("wire", k, optimized=True)
                self.assertEqual(wire["subcircuits"], (2**k + 1 if k > 1 else 3) + 4**k)
                self.assertAlmostEqual(wire["gamma"], 4**k)

    def test_general_gate_counts_its_czs(self):
        """A SWAP is cut as its 3 CZs: 2 * 5**3 circuits, gamma 27."""
        import Qdislib as qd

        b = _Builder(4, seed=9)
        b.layer(0, 1, 2, 3)
        b.gate("cz", 0, 1), b.gate("cz", 2, 3)
        swap = b.gate("swap", 1, 2)
        b.layer(0, 1, 2, 3)
        b.gate("cz", 0, 1), b.gate("cz", 2, 3)
        count = qd.circuit_cutting_overhead(b.qc, {"gate": [swap]}, optimized=True)
        self.assertEqual((count["subcircuits"], count["gamma"]), (2 * 5**3, 27))

    def test_mixed_sums_over_fragments(self):
        import Qdislib as qd

        circuit, cuts = _chain()
        # Fragments: {0,1} sends wire 1; {1',2,...} receives it, holds one gate side
        # and sends wire 2; the last receives wire 2 and holds the other gate side.
        standard = qd.circuit_cutting_overhead(circuit, cuts)
        self.assertEqual(standard["configurations"], 6 * 8 * 8)
        self.assertEqual(standard["gamma"], 3 * 4 * 4)
        optimized = qd.circuit_cutting_overhead(circuit, cuts, optimized=True)
        self.assertEqual(optimized["subcircuits"], 3 + 5 * 4 * 3 + 5 * 4)


@unittest.skipUnless(QISKIT, "qiskit not installed")
class CutValidityTest(BaseTimedTestCase):
    """A cut is valid only if its two sides land in different fragments, judged
    with every other cut applied."""

    def test_only_the_whole_combination_separates(self):
        import itertools

        import Qdislib as qd

        # q0 and q1 linked three times: each link is cut a different way, and the
        # old q1 line only comes loose when all three are cut.
        b = _Builder(2, seed=2)
        b.layer(0, 1)
        mp = b.gate("cz", 0, 1)
        b.layer(0, 1)
        had = b.gate("cz", 0, 1)
        b.layer(0, 1)
        wire = b.wire(1, 0)
        b.layer(0, 1)
        cuts = {"gate": [mp], "hadamard": [had], "wire": [wire]}
        value = qd.circuit_cutting(b.qc, cuts, "ZZ", shots=None, optimized=True)
        self.assertAlmostEqual(value, _exact(b.qc, "ZZ"), places=10)
        for size in (1, 2):
            for kinds in itertools.combinations(cuts, size):
                with self.subTest(kinds=kinds):
                    with self.assertRaises(qd.InvalidCutError):
                        qd.circuit_cutting(
                            b.qc, {k: cuts[k] for k in kinds}, "ZZ", shots=None
                        )

    def test_a_cut_that_another_path_bypasses(self):
        """An uncut CZ keeps both sides of the cut gate together; the wire cut does
        split the circuit, but the gate cut is still not valid."""
        import Qdislib as qd

        b = _Builder(3, seed=8)
        b.layer(0, 1, 2)
        b.gate("cz", 0, 1)
        gate = b.gate("cz", 0, 1)
        b.layer(0, 1)
        wire = b.wire(1, 2)
        b.layer(1, 2)
        with self.assertRaisesRegex(qd.InvalidCutError, gate):
            qd.circuit_cutting(b.qc, {"gate": [gate], "wire": [wire]}, "ZZZ")
        with self.assertRaises(qd.InvalidCutError):
            qd.gate_wire_cutting(b.qc, [gate], [wire], "ZZZ", software="qiskit")
        with self.assertRaises(qd.InvalidCutError):
            qd.gate_cutting(b.qc, [gate], observables="ZZZ")
        qd.wire_cutting(b.qc, [wire], observables="ZZZ", software="qiskit")


class CutArgumentsTest(BaseTimedTestCase):
    def test_list_and_dict_forms(self):
        from Qdislib.core.cutting_algorithms.mixed_cutting import _normalise_cuts

        self.assertEqual(
            _normalise_cuts(["CZ_3", ("H_1", "CZ_4")], False),
            (["CZ_3"], [], [("H_1", "CZ_4")]),
        )
        self.assertEqual(_normalise_cuts(["CZ_3"], True), ([], ["CZ_3"], []))
        with self.assertRaises(ValueError):
            _normalise_cuts({"gates": ["CZ_3"]}, False)

    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_gate_wire_cutting_forwards_the_new_options(self):
        import Qdislib as qd

        circuit, cuts = _chain()
        value = qd.gate_wire_cutting(
            circuit,
            cuts["gate"],
            cuts["wire"],
            "ZZZZ",
            software="qiskit",
            shots=None,
            optimized=True,
        )
        self.assertAlmostEqual(value, _exact(circuit, "ZZZZ"), places=10)


if __name__ == "__main__":
    unittest.main()
