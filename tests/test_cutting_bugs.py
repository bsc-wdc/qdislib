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

"""Regression tests for circuit-cutting correctness bugs.

Each test is written to FAIL on the unpatched library and PASS once the
corresponding fix is applied. An independent statevector oracle is used so the
tests do not depend on Qdislib's own ``analytical_solution`` helper.
"""

import unittest

import numpy as np
import qibo
from qibo import models, gates
from qiskit import QuantumCircuit

from tests import BaseTimedTestCase

qibo.set_backend("numpy")


def exact_pauli_expectation(qibo_circuit, paulis):
    """Independent oracle for a general Pauli string (supports I/X/Y/Z)."""
    psi = qibo_circuit().state(numpy=True)
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    table = {"I": I, "X": X, "Y": Y, "Z": Z}
    op = np.array([[1]], dtype=complex)
    for p in paulis:
        op = np.kron(op, table[p])
    return float(np.real(psi.conj() @ op @ psi))


def exact_z_expectation(qibo_circuit, obs):
    """Independent oracle: <prod_q Z_q for obs[q]=='Z'> from the statevector.

    Qibo convention: qubit 0 is the most-significant bit of the state index.
    Supports only 'Z'/'I' (sufficient for these tests).
    """
    psi = qibo_circuit().state(numpy=True)
    n = qibo_circuit.nqubits
    probs = np.abs(psi) ** 2
    val = 0.0
    for i in range(2 ** n):
        sign = 1
        for q, o in enumerate(obs):
            if o == "Z":
                sign *= (-1) ** ((i >> (n - 1 - q)) & 1)
        val += sign * probs[i]
    return float(val)


class CuttingBugTest(BaseTimedTestCase):

    # ----------------------------------------------------------------- #1228
    def test_qpu_dict_none_does_not_crash(self):
        """gate_cutting(qpu=True) with default qpu_dict=None must not crash."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.8, 0)
        # Before the fix this raises TypeError ('NoneType' is not iterable).
        result = gate_cutting(qc, ["CZ_2"], qpu=True)
        self.assertIsInstance(result, float)
        self.assertTrue(np.isfinite(result))

    # ----------------------------------------------------------------- #1226
    def test_qiskit_xy_observable_does_not_crash(self):
        """X/Y observables on a Qiskit circuit must not crash and must match Qibo."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.8, 0)

        qcirc = models.Circuit(2)
        qcirc.add(gates.H(0))
        qcirc.add(gates.CZ(0, 1))
        qcirc.add(gates.RY(0, 0.8))

        # Before the fix the Qiskit call raises AttributeError (Qibo-only _change_basis).
        qiskit_val = gate_cutting(qc, ["CZ_2"], observables="XX", shots=50000)
        qibo_val = gate_cutting(qcirc, ["CZ_2"], observables="XX", shots=50000)
        self.assertAlmostEqual(qiskit_val, qibo_val, delta=0.1)

    # ----------------------------------------------------------------- #1224
    def test_multi_component_no_double_count(self):
        """Disconnected components must not be multiplied in twice."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        c = models.Circuit(4)
        c.add(gates.H(0))             # 1
        c.add(gates.H(1))             # 2
        c.add(gates.CZ(0, 1))         # 3 -> CZ_3 (the cut, component A)
        c.add(gates.H(1))             # 4
        c.add(gates.RY(1, np.pi / 5)) # 5
        c.add(gates.RY(2, np.pi / 3)) # 6 component B: <Z2> = 0.5 (!= 1)
        c.add(gates.CZ(2, 3))         # 7 edge so {2,3} is one component

        analytical = exact_z_expectation(c, "ZZZZ")
        recon = gate_cutting(c, ["CZ_3"], shots=50000)
        # Before the fix recon ~= 0.5 * analytical (component B counted twice).
        self.assertAlmostEqual(recon, analytical, delta=0.05)

    # ----------------------------------------------------------------- #1227
    @unittest.expectedFailure
    def test_expec_value_qiskit_bit_order(self):
        """_expec_value_qiskit maps observable[i] to the WRONG qubit (#1227).

        This is a confirmed bug: for a bare subcircuit, observable 'ZI' (Z on
        qubit 0) returns ~0.0 instead of ~0.5 because Qiskit count keys are
        little-endian and are zipped against the qubit-0-first observable.

        Marked expectedFailure rather than fixed: the same function is reused to
        evaluate internal reconstruction subcircuits, whose keys are
        multi-register/space-separated and rely on the current (leftmost-bit)
        pairing. A blanket ``reversed(key)`` fixes this case but breaks
        reconstruction. The correct fix is a register-aware qubit->bit mapping
        that separates the two callers; see issue #1227.

        NOTE: the same ``zip(key, observables)`` pattern appears at FIVE Qiskit
        sites in gate_cutting.py (CPU ``_expec_value_qiskit`` ~781, the QPU freq
        path ~705, ``_expec_value_qiskit_qpu`` ~838, and the two GPU variants
        ~907/~973). Only the CPU path is covered here because the QPU/GPU paths
        need hardware/CUDA. A real fix must address all five together, or the
        GPU/QPU users silently keep the wrong-qubit mapping.
        """
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit

        # qubit0: <Z>=cos(pi/3)=0.5 ; qubit1: <Z>=cos(pi/2)=0.0  (no cut, no reindex)
        sub = QuantumCircuit(2)
        sub.ry(np.pi / 3, 0)
        sub.ry(np.pi / 2, 1)

        # obs_I=[1] -> "ZI" -> Z on qubit 0 -> expect 0.5
        zi = _expec_value_qiskit((sub.copy(), [1]), shots=200000)
        # obs_I=[0] -> "IZ" -> Z on qubit 1 -> expect 0.0
        iz = _expec_value_qiskit((sub.copy(), [0]), shots=200000)
        self.assertAlmostEqual(zi, 0.5, delta=0.02)
        self.assertAlmostEqual(iz, 0.0, delta=0.02)

    # ----------------------------------------------------------------- #1225
    def test_wire_cutting_forwards_shots(self):
        """User-supplied shots must reach the subcircuit evaluator."""
        from Qdislib.core.cutting_algorithms import wire_cutting as wcmod

        seen = []
        original = wcmod._expec_value_qibo

        def recorder(subcirc, shots, backend, *a, **k):
            seen.append(shots)
            return 0.0

        wcmod._expec_value_qibo = recorder
        try:
            circuit = models.Circuit(3)
            circuit.add(gates.H(0))
            circuit.add(gates.H(1))
            circuit.add(gates.CZ(0, 1))
            circuit.add(gates.H(1))
            circuit.add(gates.CZ(1, 2))
            circuit.add(gates.RY(1, np.pi / 5))
            wcmod.wire_cutting(circuit, [("CZ_3", "H_4")], shots=4242)
        finally:
            wcmod._expec_value_qibo = original

        # Before the fix every call sees the hardcoded default 10000, never 4242.
        self.assertTrue(seen, "evaluator was never called")
        self.assertTrue(all(s == 4242 for s in seen),
                        f"expected all shots==4242, saw {set(seen)}")

    # ----------------------------------------------------------------- #1229
    def test_find_cut_deterministic_and_nonempty(self):
        """find_cut must be reproducible and find a valid cut when one exists."""
        from Qdislib.core.find_cut.find_cut import find_cut

        def build():
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
            return circuit

        import random

        random.seed(1)
        np.random.seed(1)
        cut_a = find_cut(build(), 3, gate_cut=False)
        # Perturb the global RNG state: an internally-seeded search must be immune.
        random.seed(987654321)
        np.random.seed(987654321)
        cut_b = find_cut(build(), 3, gate_cut=False)
        self.assertTrue(cut_a, "find_cut returned an empty cut for a cuttable circuit")
        self.assertEqual(cut_a, cut_b,
                         "find_cut depends on global RNG state (not internally seeded)")

    # ------------------------------------------------------- light sweep
    def test_qiskit_qibo_consistency(self):
        """Identical circuit + cut must reconstruct equally via Qiskit and Qibo."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.8, 0)

        qcirc = models.Circuit(2)
        qcirc.add(gates.H(0))
        qcirc.add(gates.CZ(0, 1))
        qcirc.add(gates.RY(0, 0.8))

        v_qiskit = gate_cutting(qc, ["CZ_2"], shots=50000)
        v_qibo = gate_cutting(qcirc, ["CZ_2"], shots=50000)
        self.assertAlmostEqual(v_qiskit, v_qibo, delta=0.1)

    def test_single_cut_reconstruction_tight(self):
        """Tight-tolerance check that the gate-cut coefficient math is correct."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        qcirc = models.Circuit(2)
        qcirc.add(gates.H(0))
        qcirc.add(gates.CZ(0, 1))
        qcirc.add(gates.RY(0, 0.8))

        analytical = exact_z_expectation(qcirc, "ZZ")
        recon = gate_cutting(qcirc, ["CZ_2"], shots=200000)
        # High shots -> tight tolerance exercises the math, not just a smoke test.
        self.assertAlmostEqual(recon, analytical, delta=0.03)

    def test_reconstruction_pure_math(self):
        """gate_cutting_subcircuit_reconstruction with number_cuts=0 is a plain product."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            gate_cutting_subcircuit_reconstruction,
        )

        self.assertAlmostEqual(
            gate_cutting_subcircuit_reconstruction([0.5, 0.4], 0), 0.2, delta=1e-9
        )

    def test_wire_cutting_reconstruction_tight(self):
        """Tight-tolerance check of the wire-cut coefficient math."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.H(1))
        c.add(gates.CZ(0, 1))
        c.add(gates.H(1))
        c.add(gates.CZ(1, 2))
        c.add(gates.RY(1, np.pi / 5))
        analytical = exact_z_expectation(c, "ZZZ")
        recon = wire_cutting(c, [("CZ_3", "H_4")], shots=200000)
        self.assertAlmostEqual(recon, analytical, delta=0.03)

    def test_no_cut_returns_full_expectation(self):
        """An empty cut list must reconstruct the full-circuit expectation."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RY(0, 0.8))
        analytical = exact_z_expectation(c, "ZZ")
        self.assertAlmostEqual(gate_cutting(c, [], shots=200000), analytical, delta=0.03)

    def test_single_qubit_circuit(self):
        """A trivial single-qubit circuit with no cut is handled correctly."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        c = models.Circuit(1)
        c.add(gates.RY(0, 0.8))
        analytical = exact_z_expectation(c, "Z")
        self.assertAlmostEqual(gate_cutting(c, [], shots=200000), analytical, delta=0.03)

    def test_change_basis_does_not_mutate_input(self):
        """gate_cutting with observables must not mutate the caller's circuit."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.8, 0)
        depth_before = qc.depth()
        v1 = gate_cutting(qc, ["CZ_2"], observables="XX", shots=50000)
        self.assertEqual(qc.depth(), depth_before, "input circuit was mutated")
        # Reusing the same circuit must give the same answer (no gate accumulation).
        v2 = gate_cutting(qc, ["CZ_2"], observables="XX", shots=50000)
        self.assertAlmostEqual(v1, v2, delta=0.1)

    def test_y_observable_matches_analytic(self):
        """Y-basis change must yield the correct <Y..Y>, not just agree across backends."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        c = models.Circuit(2)
        c.add(gates.RX(0, np.pi / 2))
        c.add(gates.CZ(0, 1))
        c.add(gates.RX(1, np.pi / 2))
        analytical = exact_pauli_expectation(c, "YY")
        recon = gate_cutting(c, ["CZ_2"], observables="YY", shots=200000)
        self.assertAlmostEqual(recon, analytical, delta=0.04)

    def test_two_cut_reconstruction(self):
        """Gate cutting with k=2 cuts (6**2 subcircuit groups) reconstructs correctly."""
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        c = models.Circuit(3)
        c.add(gates.H(0))            # 1
        c.add(gates.H(1))            # 2
        c.add(gates.CZ(0, 1))        # 3 -> CZ_3
        c.add(gates.RY(1, 0.5))      # 4
        c.add(gates.CZ(1, 2))        # 5 -> CZ_5
        c.add(gates.RY(2, 0.3))      # 6
        analytical = exact_z_expectation(c, "ZZZ")
        recon = gate_cutting(c, ["CZ_3", "CZ_5"], shots=200000)
        self.assertAlmostEqual(recon, analytical, delta=0.05)

    def test_wire_cutting_multi_component_forwards_shots(self):
        """The multi-component wire-cutting path must also forward shots."""
        from Qdislib.core.cutting_algorithms import wire_cutting as wcmod

        seen = []
        original = wcmod._expec_value_qibo

        def recorder(subcirc, shots, backend, *a, **k):
            seen.append(shots)
            return 0.0

        wcmod._expec_value_qibo = recorder
        try:
            c = models.Circuit(4)
            c.add(gates.H(0))           # component A (cuttable)
            c.add(gates.H(1))
            c.add(gates.CZ(0, 1))
            c.add(gates.H(1))
            c.add(gates.CZ(1, 2))
            c.add(gates.RY(1, np.pi / 5))
            c.add(gates.RY(3, 0.4))     # component B (isolated qubit 3)
            wcmod.wire_cutting(c, [("CZ_3", "H_4")], shots=7777)
        finally:
            wcmod._expec_value_qibo = original

        self.assertTrue(seen, "evaluator was never called")
        self.assertTrue(all(s == 7777 for s in seen),
                        f"multi-component path dropped shots, saw {set(seen)}")

    def test_multiparam_gate_roundtrip(self):
        """dag_to_circuit_qibo must unpack multi-parameter gates (U3/U2/U1)."""
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag, dag_to_circuit_qibo

        for g in (gates.U3(0, 0.3, 0.4, 0.5), gates.U2(0, 0.3, 0.4), gates.U1(0, 0.3)):
            c = models.Circuit(1)
            c.add(g)
            rebuilt = dag_to_circuit_qibo(circuit_qibo_to_dag(c), 1)[0]
            fidelity = abs(np.vdot(c().state(numpy=True), rebuilt().state(numpy=True)))
            self.assertAlmostEqual(fidelity, 1.0, delta=1e-6,
                                   msg=f"{g.__class__.__name__} did not round-trip")

    def test_find_cut_accepts_qiskit_input(self):
        """find_cut must not crash on a Qiskit circuit (num_qubits vs nqubits)."""
        from Qdislib.core.find_cut.find_cut import find_cut

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)
        qc.cz(1, 2)
        cut = find_cut(qc, 2, gate_cut=False)
        self.assertTrue(cut, "find_cut returned no cut for a cuttable Qiskit circuit")

    def test_find_cut_small_graph_no_crash(self):
        """find_cut must return gracefully (not crash) on tiny/sparse circuits."""
        from Qdislib.core.find_cut.find_cut import find_cut

        # Single 2-qubit gate -> the wire (gate) graph has one node; the KL and
        # spectral partitioners must not raise on it.
        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        cut = find_cut(c, 2, gate_cut=False)
        self.assertIsInstance(cut, list)
