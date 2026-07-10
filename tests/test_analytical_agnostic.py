#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""``analytical_solution`` / ``exact_statevector`` are backend-agnostic.

The exact-expectation oracle must accept a circuit from any SDK (Qibo, Qiskit,
CUDA-Q, PennyLane) and return the SAME value -- each SDK's native endianness is
canonicalised to qubit-0-first so an ``observables`` string lines up regardless of
input. (CUDA-Q needs a GPU/Linux runtime and is covered separately in the container;
here we check the three SDKs installable on the test env.)
"""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None
QISKIT = importlib.util.find_spec("qiskit") is not None
PENNYLANE = importlib.util.find_spec("pennylane") is not None

_PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _numpy_oracle(psi, obs):
    """Independent <psi|P|psi> via an explicit Kronecker product (qubit-0-first)."""
    op = np.array([[1]], dtype=complex)
    for ch in obs:
        op = np.kron(op, _PAULI[ch])
    return float(np.real(psi.conj() @ op @ psi))


# Reference state for RY(0, 0.7); CNOT(0, 1); RX(2, 0.4) on 3 qubits, built by hand
# so the oracle does not depend on any SDK. Qubit 0 is the most-significant bit.
def _reference_state():
    a, b = np.cos(0.35), np.sin(0.35)  # RY(0.7) on qubit 0
    two = np.zeros(4, dtype=complex)
    two[0b00] = a
    two[0b11] = b  # CNOT(0,1): |0.>->|00>, |1.>->|11>
    c, d = np.cos(0.2), -1j * np.sin(0.2)  # RX(0.4) on qubit 2
    q2 = np.array([c, d], dtype=complex)
    return np.kron(two, q2)


_OBS = ["ZZZ", "ZIZ", "IZI", "IIZ", "ZII", "XYZ", "ZXI", "YIY"]


class AgnosticAnalyticalSolutionTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        if QIBO:
            import qibo

            qibo.set_backend("numpy")

    def _qibo(self):
        from qibo import models, gates

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.7))
        c.add(gates.CNOT(0, 1))
        c.add(gates.RX(2, 0.4))
        return c

    def _qiskit(self):
        from qiskit import QuantumCircuit

        c = QuantumCircuit(3)
        c.ry(0.7, 0)
        c.cx(0, 1)
        c.rx(0.4, 2)
        return c

    def _pennylane(self):
        import pennylane as qml

        with qml.tape.QuantumTape() as tape:
            qml.RY(0.7, wires=0)
            qml.CNOT(wires=[0, 1])
            qml.RX(0.4, wires=2)
        return tape

    @unittest.skipUnless(QIBO, "qibo not installed")
    def test_qibo_matches_oracle(self):
        import Qdislib as qd

        psi = _reference_state()
        for obs in _OBS:
            self.assertAlmostEqual(
                qd.analytical_solution(self._qibo(), obs),
                _numpy_oracle(psi, obs),
                9,
                obs,
            )

    @unittest.skipUnless(QIBO and QISKIT, "needs qibo + qiskit")
    def test_qiskit_matches_oracle(self):
        import Qdislib as qd

        psi = _reference_state()
        for obs in _OBS:
            self.assertAlmostEqual(
                qd.analytical_solution(self._qiskit(), obs),
                _numpy_oracle(psi, obs),
                9,
                obs,
            )

    @unittest.skipUnless(QIBO and PENNYLANE, "needs qibo + pennylane")
    def test_pennylane_matches_oracle(self):
        import Qdislib as qd

        psi = _reference_state()
        for obs in _OBS:
            self.assertAlmostEqual(
                qd.analytical_solution(self._pennylane(), obs),
                _numpy_oracle(psi, obs),
                9,
                obs,
            )

    @unittest.skipUnless(QIBO and QISKIT and PENNYLANE, "needs all three SDKs")
    def test_all_backends_agree(self):
        import Qdislib as qd

        for obs in _OBS + [[(0.5, "ZZZ"), (0.3, "XIX")]]:  # incl. a weighted sum
            a = qd.analytical_solution(self._qibo(), obs)
            b = qd.analytical_solution(self._qiskit(), obs)
            c = qd.analytical_solution(self._pennylane(), obs)
            self.assertAlmostEqual(a, b, 9, f"qibo vs qiskit {obs}")
            self.assertAlmostEqual(a, c, 9, f"qibo vs pennylane {obs}")

    @unittest.skipUnless(QIBO and QISKIT, "needs qibo + qiskit")
    def test_exact_statevector_qubit0_first(self):
        from Qdislib.utils.circuit import exact_statevector

        psi = _reference_state()
        # up to a global phase; RY/RX/CNOT keep it real here, so compare directly.
        self.assertTrue(np.allclose(exact_statevector(self._qibo()), psi, atol=1e-9))
        self.assertTrue(np.allclose(exact_statevector(self._qiskit()), psi, atol=1e-9))

    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_scales_beyond_dense_matrix(self):
        # A 16-qubit dense 2**16 x 2**16 operator would be ~68 GB; the per-qubit
        # oracle must handle it from the 2**16 statevector without OOM.
        import Qdislib as qd
        from qiskit import QuantumCircuit

        n = 16
        qc = QuantumCircuit(n)
        for i in range(n):
            qc.ry(0.1 * (i + 1), i)
        for i in range(n - 1):
            qc.cz(i, i + 1)
        value = qd.analytical_solution(qc, "Z" * n)
        self.assertIsInstance(value, float)
        self.assertLessEqual(abs(value), 1.0 + 1e-9)

    def test_unsupported_type_raises(self):
        from Qdislib.utils.circuit import exact_statevector

        with self.assertRaises(TypeError):
            exact_statevector(object())


if __name__ == "__main__":
    unittest.main()
