#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Backend-agnostic DAG conversion: the neutral DAG must reproduce the same state on
every SDK (Qiskit / PennyLane / CUDA-Q), including the Hadamard-cut gadget gate set
(H/S/U3/CU3/U1/CNOT) whose ``U3``/``CU3`` carry SDK-specific phase conventions."""

import importlib.util
import unittest

import numpy as np

from unittest import mock

from tests import BaseTimedTestCase

_HAS = lambda m: importlib.util.find_spec(m) is not None  # noqa: E731
QIBO = _HAS("qibo")


def _gate_cases():
    """(label, builder) circuits exercising the gadget + common gates. Built lazily so
    the module imports without qibo."""
    from qibo import gates, models

    def circ(*ops):
        n = 1 + max(max(o.qubits) for o in ops)
        c = models.Circuit(n)
        for o in ops:
            c.add(o)
        return c

    return [
        ("single", circ(gates.H(0), gates.RX(0, 0.5), gates.RZ(0, 0.3))),
        ("U3", circ(gates.H(0), gates.U3(0, 0.4, 0.9, 1.2))),
        ("U1", circ(gates.H(0), gates.U1(0, 0.7))),
        ("CU3", circ(gates.H(0), gates.CU3(0, 1, 0.3, 0.7, 1.1))),
        ("CU1", circ(gates.H(0), gates.CU1(0, 1, 0.6))),
        ("iSWAP", circ(gates.H(0), gates.RX(1, 0.4), gates.iSWAP(0, 1))),
        ("RZZ", circ(gates.H(0), gates.RX(1, 0.4), gates.RZZ(0, 1, 0.5))),
        ("SWAP", circ(gates.H(0), gates.RX(1, 0.5), gates.SWAP(0, 1))),
        (
            "gadget",  # the exact gate set the Hadamard-cut fragment emits
            circ(
                gates.H(0),
                gates.S(1),
                gates.U3(0, 0.4, 0.9, 1.2),
                gates.CU3(0, 1, 0.3, 0.7, 1.1),
                gates.U1(1, 0.6),
                gates.CNOT(0, 1),
            ),
        ),
    ]


def _fidelity(ref, sv):
    ref, sv = np.asarray(ref), np.asarray(sv)
    return abs(np.vdot(ref, sv)) / (np.linalg.norm(ref) * np.linalg.norm(sv))


@unittest.skipUnless(QIBO, "qibo not installed")
class AgnosticConverterTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _ref(self, circuit):
        return np.asarray(circuit().state())

    @unittest.skipUnless(_HAS("qiskit"), "qiskit not installed")
    def test_qiskit_matches_qibo(self):
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Statevector
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        for label, c in _gate_cases():
            n = c.nqubits
            qc = dag_to_circuit_qiskit(circuit_qibo_to_dag(c), n)
            qc = qc[0] if isinstance(qc, (list, tuple)) else qc
            plain = QuantumCircuit(n)  # drop the classical measure register
            for instr in qc.data:
                if instr.operation.name in ("measure", "barrier"):
                    continue
                plain.append(
                    instr.operation, [qc.find_bit(q).index for q in instr.qubits]
                )
            sv = np.asarray(Statevector(plain))  # qiskit is little-endian
            sv = sv.reshape([2] * n).transpose(list(range(n - 1, -1, -1))).flatten()
            self.assertGreater(_fidelity(self._ref(c), sv), 1 - 1e-6, label)

    @unittest.skipUnless(_HAS("pennylane"), "pennylane not installed")
    def test_pennylane_matches_qibo(self):
        import networkx as nx
        import pennylane as qml
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        from Qdislib.utils.graph_pennylane import _apply_pennylane_gate

        for label, c in _gate_cases():
            n = c.nqubits
            dag = circuit_qibo_to_dag(c)
            ops = [
                (
                    dag.nodes[nd]["gate"].upper(),
                    tuple(dag.nodes[nd]["qubits"]),
                    list(dag.nodes[nd]["parameters"] or []),
                )
                for nd in nx.topological_sort(dag)
            ]
            dev = qml.device("default.qubit", wires=n)

            @qml.qnode(dev)
            def qf():
                for gname, qb, pr in ops:
                    _apply_pennylane_gate(gname, qb, pr)
                return qml.state()

            self.assertGreater(
                _fidelity(self._ref(c), np.asarray(qf())), 1 - 1e-6, label
            )

    @unittest.skipUnless(_HAS("cudaq"), "cudaq not installed")
    def test_cudaq_matches_qibo(self):
        import cudaq
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel

        for label, c in _gate_cases():
            n = c.nqubits
            r = dag_to_circuit_cudaq_kernel(circuit_qibo_to_dag(c), n)
            kernel = r[0] if isinstance(r, (list, tuple)) else r
            sv = np.asarray(cudaq.get_state(kernel))  # cudaq is little-endian here
            if n > 1:
                sv = sv.reshape([2] * n).transpose(list(range(n - 1, -1, -1))).flatten()
            self.assertGreater(_fidelity(self._ref(c), sv), 1 - 1e-6, label)


@unittest.skipUnless(_HAS("cudaq"), "cudaq not installed")
class CudaqKernelSourceTest(BaseTimedTestCase):
    """A kernel is read from its source: what cannot be read must not be left out,
    which would cut, and answer for, another circuit."""

    def test_what_the_source_does_not_spell_out_is_refused(self):
        import cudaq

        import Qdislib as qd

        @cudaq.kernel
        def loop():
            q = cudaq.qvector(2)
            for i in range(2):
                h(q[i])  # noqa: F821
            x.ctrl(q[0], q[1])  # noqa: F821

        @cudaq.kernel
        def register():
            q = cudaq.qvector(2)
            h(q)  # noqa: F821
            x.ctrl(q[0], q[1])  # noqa: F821

        @cudaq.kernel
        def variable():
            q = cudaq.qvector(2)
            theta = 0.7
            rx(theta, q[0])  # noqa: F821
            x.ctrl(q[0], q[1])  # noqa: F821

        for kernel, why in (
            (loop, "a loop"),
            (register, "not a number or q"),
            (variable, "not a number or q"),
        ):
            with self.assertRaisesRegex(ValueError, why):
                qd.gate_cutting(kernel, [], "ZZ", shots=None)

    def test_adjoint_and_measurements_are_read(self):
        import cudaq

        from Qdislib.utils.graph_cudaq import circuit_cudaq_to_dag

        @cudaq.kernel
        def kernel():
            q = cudaq.qvector(2)
            h(q[0])  # noqa: F821
            t.adj(q[0])  # noqa: F821
            s.adj(q[1])  # noqa: F821
            x.ctrl(q[0], q[1])  # noqa: F821
            mz(q)  # noqa: F821

        gates = [d["gate"] for _, d in circuit_cudaq_to_dag(kernel).nodes(data=True)]
        self.assertEqual(gates, ["h", "tdg", "sdg", "cx"])

    def test_a_qubit_no_gate_touches_is_still_a_qubit(self):
        import cudaq

        import Qdislib as qd

        @cudaq.kernel
        def kernel():
            q = cudaq.qvector(3)
            x(q[0])  # noqa: F821

        self.assertAlmostEqual(qd.analytical_solution(kernel, "ZII"), -1.0)

    def test_a_compiled_kernel_leaves_no_files(self):
        import os
        import tempfile

        from Qdislib.utils import graph_cudaq

        source = (
            "import cudaq\n"
            "@cudaq.kernel\n"
            "def cut_kernel():\n"
            "    q = cudaq.qvector(1)\n"
            "    x(q[0])\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(tempfile, "tempdir", tmp):
                graph_cudaq._CUDAQ_RUN_KERNEL_CACHE.clear()
                kernel = graph_cudaq._compile_cudaq_source(source)
                self.assertEqual(os.listdir(tmp), [])
        import cudaq

        self.assertEqual(dict(cudaq.sample(kernel, shots_count=5).items()), {"1": 5})


@unittest.skipUnless(_HAS("pennylane"), "pennylane not installed")
class PennylaneWiresTest(BaseTimedTestCase):
    def test_a_wire_no_gate_touches_counts(self):
        # with the wires counted, "Z" * 2 read Z0 Z1 of a gate on wires 0 and 2
        import pennylane as qml

        import Qdislib as qd

        with qml.tape.QuantumTape() as tape:
            qml.PauliX(0)
            qml.CNOT(wires=[0, 2])
        self.assertEqual(qd.utils.backends.num_qubits(tape), 3)
        self.assertAlmostEqual(qd.circuit_cutting(tape, [], shots=None), 1.0)
        self.assertAlmostEqual(qd.analytical_solution(tape, "ZIZ"), 1.0)
        self.assertAlmostEqual(qd.analytical_solution(tape, "ZII"), -1.0)


def _subclassed(num_qubits):
    """A circuit of a class from another package (qiskit-nature's ansatzes, a
    user's own)."""
    from qiskit import QuantumCircuit

    class MyCircuit(QuantumCircuit):
        pass

    return MyCircuit(num_qubits)


@unittest.skipUnless(_HAS("qiskit") and _HAS("qiskit_aer"), "qiskit required")
class SubclassedCircuitTest(BaseTimedTestCase):
    def _circuit(self):
        circuit = _subclassed(3)
        circuit.h(0)
        circuit.cx(0, 1)
        circuit.cz(1, 2)
        circuit.ry(0.4, 2)
        return circuit

    def test_is_read_as_a_circuit_of_its_sdk(self):
        import Qdislib as qd

        circuit = self._circuit()
        exact = qd.analytical_solution(circuit, "ZZZ")
        self.assertAlmostEqual(
            qd.circuit_cutting(circuit, ["CZ_3"], "ZZZ", shots=None), exact, places=9
        )
        self.assertEqual(qd.find_cut(circuit, max_cuts=0), [])

    def test_native_form_of_expectation_value(self):
        import Qdislib as qd

        circuit = _subclassed(1)
        circuit.x(0)
        self.assertEqual(qd.expectation_value(circuit, shots=100, seed=1), -1.0)


if __name__ == "__main__":
    unittest.main()
