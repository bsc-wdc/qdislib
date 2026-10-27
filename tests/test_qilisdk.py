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

"""QiliSDK circuits in, QiliSDK subcircuits out, and backend="qilisdk"."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QILISDK = importlib.util.find_spec("qilisdk") is not None
QISKIT = importlib.util.find_spec("qiskit") is not None


def _pair(seed=1):
    """The same 4-qubit circuit in QiliSDK and in Qiskit (gate 11 joins the halves)."""
    from qilisdk.digital import CNOT, CZ, RY, Circuit, M
    from qiskit import QuantumCircuit

    rng = np.random.default_rng(seed)
    qili, qk = Circuit(4), QuantumCircuit(4)

    def layer():
        for q in range(4):
            t = float(rng.uniform(0.2, 2.8))
            qili.add(RY(q, theta=t))
            qk.ry(t, q)

    layer()
    qili.add(CZ(0, 1))
    qk.cz(0, 1)
    qili.add(CNOT(2, 3))
    qk.cx(2, 3)
    layer()
    qili.add(CNOT(1, 2))  # CNOT_11
    qk.cx(1, 2)
    layer()
    qili.add(M(0, 1, 2, 3))
    return qili, qk


def _qili_state(circuit):
    from qilisdk.backends import QiliSim
    from qilisdk.functionals import DigitalPropagation
    from qilisdk.readout import Readout

    state = (
        QiliSim()
        .execute(DigitalPropagation(circuit), Readout().with_state_tomography())
        .get_state()
    )
    dense = state.dense() if callable(getattr(state, "dense", None)) else state.dense
    return np.asarray(dense).ravel()


@unittest.skipUnless(QILISDK and QISKIT, "qilisdk and qiskit needed")
class QiliSDKTest(BaseTimedTestCase):
    def test_every_gate_converts_exactly(self):
        """Qiskit -> DAG -> QiliSDK keeps the state (decomposed gates included)."""
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Statevector

        from Qdislib.utils.backends import convert_to_dag
        from Qdislib.utils.graph_qilisdk import dag_to_circuit_qilisdk

        qc = QuantumCircuit(3)
        rng = np.random.default_rng(0)
        a = lambda: float(rng.uniform(0.2, 2.8))  # noqa: E731
        for q in range(3):
            qc.ry(a(), q)
        qc.h(0), qc.sdg(1), qc.tdg(2), qc.sx(0), qc.p(a(), 1), qc.u(a(), a(), a(), 2)
        qc.cy(0, 1), qc.ch(1, 2), qc.cp(a(), 2, 0), qc.crz(a(), 0, 2)
        qc.cu(a(), a(), a(), a(), 1, 0), qc.rzz(a(), 0, 1), qc.rxx(a(), 1, 2)
        qc.ryy(a(), 0, 2), qc.iswap(0, 1), qc.ccx(0, 1, 2), qc.swap(1, 2)
        want = Statevector(qc).reverse_qargs().data  # QiliSDK: qubit 0 first
        circuit, _ = dag_to_circuit_qilisdk(convert_to_dag(qc), 3)
        self.assertAlmostEqual(abs(np.vdot(want, _qili_state(circuit))), 1.0, places=10)

    def test_qilisdk_circuit_to_dag(self):
        from Qdislib.utils.backends import convert_to_dag, iter_gates, num_qubits

        qili, qk = _pair()
        dag = convert_to_dag(qili)
        self.assertIn("CNOT_11", dag)
        self.assertEqual(dag.nodes["CNOT_11"]["gate"], "cx")
        self.assertEqual(num_qubits(qili), 4)
        names = {f"{n}_{i}" for i, (n, q) in enumerate(iter_gates(qili), 1) if q}
        self.assertEqual(names, set(dag.nodes))  # the final M is no node
        same = convert_to_dag(qk)
        self.assertEqual(
            sorted(d["gate"] for _, d in dag.nodes(data=True)),
            sorted(d["gate"] for _, d in same.nodes(data=True)),
        )

    def test_cutting_a_qilisdk_circuit(self):
        import Qdislib as qd
        from qiskit.quantum_info import SparsePauliOp, Statevector

        qili, qk = _pair()
        exact = Statevector(qk).expectation_value(SparsePauliOp("ZZZZ")).real
        self.assertAlmostEqual(
            qd.gate_cutting(qili, ["CNOT_11"], shots=None), exact, 10
        )
        for value in (
            qd.gate_cutting(qili, ["CNOT_11"], shots=20000, seed=1),
            qd.circuit_cutting(qili, {"gate": ["CNOT_11"]}, shots=20000, seed=1),
            qd.gate_cutting_hadamard(qili, ["CX_11"], shots=20000, seed=1),
            qd.wire_cutting(qili, [("RY_8", "CNOT_11")], shots=20000, seed=1),
            qd.gate_cutting(qk, ["CX_11"], backend="qilisdk", shots=20000, seed=1),
        ):
            self.assertLess(abs(value - exact), 0.05)
        self.assertEqual(qd.find_cut(qili, max_qubits=2), ["CNOT_11"])

    def test_qilisdk_subcircuits(self):
        import Qdislib as qd
        from qiskit.quantum_info import SparsePauliOp, Statevector

        _, qk = _pair()
        exact = Statevector(qk).expectation_value(SparsePauliOp("ZZZZ")).real
        subs = qd.circuit_cutting_subcircuits(qk, {"gate": ["CX_11"]})
        circuits = subs.circuits("qilisdk")  # a Qiskit circuit's cut, in QiliSDK
        self.assertTrue(all(type(c).__module__.startswith("qilisdk") for c in circuits))
        values = [
            qd.expectation_value(s, shots=20000, seed=i, backend="qilisdk")
            for i, s in enumerate(subs)
        ]
        self.assertLess(abs(subs.reconstruct(values) - exact), 0.05)


@unittest.skipUnless(QILISDK and QISKIT, "qilisdk and qiskit needed")
class QiliSDKExactTest(BaseTimedTestCase):
    def test_exact_on_qilisim(self):
        """shots=None on QiliSDK subcircuits simulates them on QiliSim itself."""
        import Qdislib as qd
        from qiskit.quantum_info import SparsePauliOp, Statevector

        qili, qk = _pair()
        exact = Statevector(qk).expectation_value(SparsePauliOp("ZZZZ")).real
        for cuts in ({"gate": ["CNOT_11"]}, {"hadamard": ["CNOT_11"]}):
            value = qd.circuit_cutting(qili, cuts, backend="qilisdk", shots=None)
            self.assertAlmostEqual(value, exact, places=9)
        self.assertAlmostEqual(qd.analytical_solution(qili, "ZZZZ"), exact, places=9)
        again = qd.from_qasm(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n',
            software="qilisdk",
        )
        self.assertEqual(type(again).__module__.split(".")[0], "qilisdk")
        self.assertAlmostEqual(qd.analytical_solution(again, "ZZ"), 1.0, places=9)


class RegistryTest(BaseTimedTestCase):
    def test_every_sdk_is_registered(self):
        from Qdislib.utils.sdks import sdk, sdk_names

        self.assertEqual(
            sdk_names(), ("qibo", "qiskit", "cudaq", "pennylane", "qilisdk")
        )
        self.assertEqual(sdk("qilisdk").backend, "QILISDK_BACKEND")


if __name__ == "__main__":
    unittest.main()
