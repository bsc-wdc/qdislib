#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""General (non-CZ/CX) two-qubit gate cutting in the measure-and-prepare method.

A general gate is decomposed into CZ + single-qubit gates (via qiskit's CZ-basis
decomposer) and its CZs are cut -- so ``gate_cutting`` handles any two-qubit gate, the
same way the Hadamard cut reuses the KAK decomposition."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None
QISKIT = importlib.util.find_spec("qiskit") is not None


@unittest.skipUnless(QIBO and QISKIT, "needs qibo + qiskit")
class DecomposeToCZTest(BaseTimedTestCase):
    """The general-gate -> CZ + single-qubit decomposition (no cutting)."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _reconstruct(self, ops):
        from qibo import gates, models

        c = models.Circuit(2)
        for name, qubits, params in ops:
            cls = getattr(gates, name)
            c.add(cls(*qubits, *params) if params else cls(*qubits))
        import qibo

        return np.asarray(c.unitary(qibo.get_backend()))

    def test_decomposes_common_gates(self):
        import qibo
        from qibo import gates
        from Qdislib.core.cutting_algorithms.gate_cutting import _decompose_gate_to_cz

        be = qibo.get_backend()
        for label, gate in [
            ("iSWAP", gates.iSWAP(0, 1)),
            ("SWAP", gates.SWAP(0, 1)),
            ("fSim", gates.fSim(0, 1, 0.4, 0.7)),
            ("CRX", gates.CRX(0, 1, 0.8)),
        ]:
            matrix = np.asarray(gate.matrix(be))
            ops = _decompose_gate_to_cz(matrix)
            rebuilt = self._reconstruct(ops)
            fid = abs(np.trace(rebuilt.conj().T @ matrix)) / 4
            self.assertGreater(fid, 1 - 1e-6, label)
            self.assertTrue(any(o[0] == "CZ" for o in ops), label)


@unittest.skipUnless(QIBO and QISKIT, "needs qibo + qiskit")
class GeneralGateCuttingTest(BaseTimedTestCase):
    """Cutting a general gate reconstructs the analytical expectation (all observables,
    including identity on a cut qubit -- the 3-CZ decomposition exercises the
    multiple-cut-per-qubit sign recovery)."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _iswap_circuit(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.9))
        c.add(gates.RY(1, 0.4))
        c.add(gates.iSWAP(0, 1))
        c.add(gates.RX(0, 0.3))
        return c

    def _fsim_circuit(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.6))
        c.add(gates.RY(1, 0.9))
        c.add(gates.fSim(0, 1, 0.4, 0.7))
        c.add(gates.RZ(1, 0.5))
        return c

    def test_iswap_all_in_one(self):
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        for obs in ("ZZ", "ZI", "IZ"):
            analytical = qd.analytical_solution(self._iswap_circuit(), obs)
            got = gate_cutting(
                self._iswap_circuit(),
                ["ISWAP_3"],
                observables=obs,
                software="qibo",
                shots=1024,
            )
            self.assertLess(abs(analytical - got), 0.2, f"iSWAP {obs}")

    def test_fsim_three_cz_all_observables(self):
        # fSim needs 3 CZs -> 3 cuts on the same qubit pair; identity on a cut qubit
        # would silently mis-reconstruct without the per-cut sign recovery.
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting

        for obs in ("ZZ", "ZI", "IZ"):
            analytical = qd.analytical_solution(self._fsim_circuit(), obs)
            got = gate_cutting(
                self._fsim_circuit(),
                ["FSIM_3"],
                observables=obs,
                software="qibo",
                shots=1024,
            )
            self.assertLess(abs(analytical - got), 0.2, f"fSim {obs}")

    def test_general_gate_three_step(self):
        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(
            self._iswap_circuit(), ["ISWAP_3"], software="qibo"
        )
        values = [qd.expectation_value(s, shots=1024) for s in subs]
        analytical = qd.analytical_solution(self._iswap_circuit(), "ZZ")
        self.assertLess(abs(analytical - subs.reconstruct(values)), 0.2)

    def test_general_gate_in_mixed_cutting(self):
        # A general gate cut combined with a wire cut on the same circuit.
        import Qdislib as qd
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.iSWAP(0, 1))  # ISWAP_2 -> general gate cut
        c.add(gates.RY(1, 0.5))  # RY_3
        c.add(gates.CZ(1, 2))  # CZ_4
        c.add(gates.RY(2, 0.3))
        analytical = qd.analytical_solution(c, "ZZZ")
        got = qd.gate_wire_cutting(
            c,
            gates_cut=["ISWAP_2"],
            wire_cut=[("RY_3", "CZ_4")],
            observable="ZZZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.25)


if __name__ == "__main__":
    unittest.main()
