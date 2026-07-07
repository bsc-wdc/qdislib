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

"""Utils Graph Unit tests -- the qibo circuit <-> neutral DAG round trip."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class GateResolverTest(BaseTimedTestCase):
    """The case-insensitive gate-class resolver (mixed-case qibo classes)."""

    def test_exact_and_case_insensitive(self):
        from qibo import gates
        from Qdislib.utils.graph_qibo import _resolve_qibo_gate

        # Exact match (the common path).
        self.assertIs(_resolve_qibo_gate("H"), gates.H)
        # Mixed-case classes qibo exposes as iSWAP / fSim -- an upper- or lower-cased
        # name (as the DAG pipeline produces) must still resolve to the real class.
        self.assertIs(_resolve_qibo_gate("ISWAP"), gates.iSWAP)
        self.assertIs(_resolve_qibo_gate("iswap"), gates.iSWAP)
        self.assertIs(_resolve_qibo_gate("FSIM"), gates.fSim)

    def test_unknown_gate_raises(self):
        from Qdislib.utils.graph_qibo import _resolve_qibo_gate

        with self.assertRaises(AttributeError):
            _resolve_qibo_gate("NotAGate")


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class RoundTripTest(BaseTimedTestCase):
    """circuit -> DAG -> circuit must preserve the state for every gate type,
    including multi-parameter (U2/U3/fSim) and mixed-case (iSWAP/fSim) gates that a
    naive ``getattr(gates, name.upper())(*qubits, params)`` silently mangled."""

    def _fidelity(self, circuit):
        import qibo
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag, dag_to_circuit_qibo

        qibo.set_backend("numpy")
        dag = circuit_qibo_to_dag(circuit)
        rebuilt = dag_to_circuit_qibo(dag, circuit.nqubits)
        rebuilt = rebuilt[0] if isinstance(rebuilt, (list, tuple)) else rebuilt
        before = np.asarray(circuit().state())
        after = np.asarray(rebuilt().state())
        return abs(np.vdot(before, after))  # global-phase-insensitive

    def _assert_roundtrips(self, circuit, label):
        self.assertGreater(self._fidelity(circuit), 1 - 1e-9, label)

    def _circuit(self, *added):
        from qibo import models

        n = 1 + max((max(g.qubits) for g in added), default=0)
        c = models.Circuit(n)
        for g in added:
            c.add(g)
        return c

    def test_parameterless_and_single_param(self):
        from qibo import gates

        c = self._circuit(
            gates.H(0),
            gates.X(1),
            gates.S(0),
            gates.T(1),
            gates.RX(0, 0.3),
            gates.RY(1, 0.7),
            gates.RZ(0, 1.1),
            gates.CNOT(0, 1),
            gates.CZ(0, 1),
        )
        self._assert_roundtrips(c, "parameterless + single-param")

    def test_multi_parameter_gates(self):
        # These were broken: parameters were passed as one packed tuple.
        from qibo import gates

        for label, g in [
            ("U2", gates.U2(1, 0.4, 0.9)),
            ("U3", gates.U3(1, 0.4, 0.9, 1.2)),
            ("CU3", gates.CU3(0, 1, 0.4, 0.9, 1.2)),
            ("fSim", gates.fSim(0, 1, 0.3, 0.6)),
        ]:
            self._assert_roundtrips(self._circuit(gates.H(0), g), label)

    def test_mixed_case_gate_names(self):
        # These were broken: name.upper() -> "ISWAP"/"FSIM" miss the real classes.
        from qibo import gates

        self._assert_roundtrips(self._circuit(gates.H(0), gates.iSWAP(0, 1)), "iSWAP")
        self._assert_roundtrips(
            self._circuit(gates.H(0), gates.fSim(0, 1, 0.3, 0.6)), "fSim"
        )

    def test_assorted_two_qubit_gates(self):
        from qibo import gates

        for label, g in [
            ("CU1", gates.CU1(0, 1, 0.7)),
            ("RXX", gates.RXX(0, 1, 0.5)),
            ("RZZ", gates.RZZ(0, 1, 0.5)),
            ("SWAP", gates.SWAP(0, 1)),
            ("FSWAP", gates.FSWAP(0, 1)),
            ("GIVENS", gates.GIVENS(0, 1, 0.4)),
            ("CRX", gates.CRX(0, 1, 0.8)),
        ]:
            self._assert_roundtrips(self._circuit(gates.H(0), g), label)

    def test_three_qubit_gate(self):
        from qibo import gates

        c = self._circuit(gates.H(0), gates.H(1), gates.TOFFOLI(0, 1, 2))
        self._assert_roundtrips(c, "TOFFOLI")


if __name__ == "__main__":
    unittest.main()
