#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Regression tests for two cutting-correctness fixes.

1. ``wire_cutting`` forwards the user's ``shots``/``backend`` into subcircuit
   evaluation on the standard (non-optimized) path (they were dropped, so every
   subcircuit ran at the hardcoded default).
2. A gate/wire cut that isolates a qubit into a component whose only marker is an
   Observable-I node (stripped before conversion) left an EMPTY component, which
   crashed ``_update_qubits_serie`` with ``ValueError: max() arg is empty``. An
   empty component contributes <I> = 1 and must reconstruct without crashing.
"""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None
QISKIT = importlib.util.find_spec("qiskit") is not None


def _z_oracle(qibo_circuit, obs):
    """Independent statevector oracle for a Z/I Pauli string (qubit-0-first)."""
    psi = qibo_circuit().state(numpy=True)
    z = np.diag([1.0, -1.0]).astype(complex)
    eye = np.eye(2)
    op = np.array([[1]], dtype=complex)
    for ch in obs:
        op = np.kron(op, z if ch == "Z" else eye)
    return float(np.real(psi.conj() @ op @ psi))


@unittest.skipUnless(QIBO, "qibo not installed")
class WireCuttingForwardsShotsTest(BaseTimedTestCase):
    """The standard wire-cutting path must honour the user's ``shots`` argument."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _spy_shots(self, circuit, cut, user_shots):
        import Qdislib.core.cutting_algorithms.wire_cutting as wc

        seen = []
        original = wc._build_and_run_subcircuit

        def spy(s, software, shots, backend, qpu):
            seen.append(shots)
            return original(s, software, shots, backend, qpu)

        wc._build_and_run_subcircuit = spy
        try:
            wc.wire_cutting(circuit, cut, shots=user_shots, software="qibo")
        finally:
            wc._build_and_run_subcircuit = original
        return seen

    def test_single_component_forwards_shots(self):
        from qibo import models, gates

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))  # CZ_2
        c.add(gates.RY(1, 0.5))  # RY_3
        c.add(gates.CZ(1, 2))
        seen = self._spy_shots(c, [("CZ_2", "RY_3")], user_shots=7777)
        self.assertTrue(seen, "no subcircuits were evaluated")
        # Every subcircuit must run at the user's shot count, not the old hardcoded
        # _generate_wire_cutting default (10000).
        self.assertEqual(set(seen), {7777})

    def test_multi_component_forwards_shots(self):
        from qibo import models, gates

        c = models.Circuit(4)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))  # CZ_2  (component A: q0, q1)
        c.add(gates.RY(1, 0.4))  # RY_3
        c.add(gates.H(2))
        c.add(gates.CZ(2, 3))  # component B: q2, q3
        seen = self._spy_shots(c, [("CZ_2", "RY_3")], user_shots=7777)
        self.assertTrue(seen, "no subcircuits were evaluated")
        self.assertEqual(set(seen), {7777})


@unittest.skipUnless(QIBO, "qibo not installed")
class IsolatedQubitAfterCutTest(BaseTimedTestCase):
    """A cut that isolates a qubit must not crash when that side is observed with a
    single-qubit Z (previously an empty component hit ``max() arg is empty``)."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _circuit(self):
        from qibo import models, gates

        # Cutting CZ(0, 2) isolates qubit 1 (carries RY_2, no other edge) as its own
        # component; observing "ZII"/"IZI" leaves the other component's observable
        # all-identity -> an empty subcircuit after Observable-I stripping.
        c = models.Circuit(3)
        c.add(gates.RY(0, 0.9))  # RY_1
        c.add(gates.RY(1, 1.1))  # RY_2 -> isolated component after the cut
        c.add(gates.CZ(0, 2))  # CZ_3 -> gate cut
        return c

    def test_single_z_on_isolated_side_reconstructs(self):
        import Qdislib as qd

        for obs in ("ZII", "IZI", "IIZ", "ZIZ", "ZZZ"):
            oracle = _z_oracle(self._circuit(), obs)
            got = qd.gate_cutting(
                self._circuit(), ["CZ_3"], observables=obs, software="qibo", shots=8192
            )
            self.assertLess(abs(got - oracle), 0.15, f"qibo {obs}")

    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_single_z_on_isolated_side_reconstructs_qiskit(self):
        import Qdislib as qd

        for obs in ("ZII", "IZI"):  # the two that used to crash
            oracle = _z_oracle(self._circuit(), obs)
            got = qd.gate_cutting(
                self._circuit(),
                ["CZ_3"],
                observables=obs,
                software="qiskit",
                shots=8192,
            )
            self.assertLess(abs(got - oracle), 0.15, f"qiskit {obs}")

    def test_update_qubits_serie_empty_component(self):
        # Unit check: an empty component reports one idle |0> qubit (contributes 1.0)
        # instead of raising on max()/min() of an empty set.
        import networkx
        from Qdislib.utils.graph_qibo import _update_qubits_serie

        s, highest, smallest = _update_qubits_serie(networkx.DiGraph())
        self.assertEqual((highest, smallest), (1, 0))


if __name__ == "__main__":
    unittest.main()
