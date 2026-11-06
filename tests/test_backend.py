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

"""The pluggable backend layer: registry, converter, and QdBackEnd routing."""

import importlib.util
import types
import unittest
from unittest import mock

from tests import BaseTimedTestCase

from Qdislib.core.backend import registry
from Qdislib.core.backend.base_backend import BaseBackend, SampleJobResult

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None
AER_AVAILABLE = importlib.util.find_spec("qiskit_aer") is not None
CUDAQ_AVAILABLE = importlib.util.find_spec("cudaq") is not None

_TOL = 0.05


def _qibo_circuit():
    """RY(0,0.8), CZ(0,1), RY(1,0.5) -- a 2-qubit circuit with <ZZ> ~ 0.611."""
    from qibo import models, gates

    c = models.Circuit(2)
    c.add(gates.RY(0, 0.8))
    c.add(gates.CZ(0, 1))
    c.add(gates.RY(1, 0.5))
    return c


class RegistryTest(BaseTimedTestCase):
    def test_default_backends_registered(self):
        from Qdislib.core.backend import available_backends

        names = available_backends()
        for expected in ("AER_BACKEND", "QIBO_BACKEND", "CUDAQ_BACKEND", "IQM_BACKEND"):
            self.assertIn(expected, names)

    def test_register_and_build_unknown(self):
        from Qdislib.core.backend import (
            register_backend,
            available_backends,
            build_backend,
        )
        from Qdislib.core.backend.base_backend import BackendConfig

        register_backend("DUMMY_BACKEND", "some.module", "SomeClass")
        self.assertIn("DUMMY_BACKEND", available_backends())
        with self.assertRaises(ValueError):
            build_backend("NOPE_BACKEND", BackendConfig(num_qubits=1))


@unittest.skipUnless(QIBO_AVAILABLE and QISKIT_AVAILABLE, "qibo+qiskit required")
class _Works(BaseBackend):
    """Samples anything."""

    def __init__(self, config):
        super().__init__(config.extra.get("name", "WORKS"), config)

    def _execute_samples(self, circuit, transpile, shots, job):
        return SampleJobResult({"0": shots}, self.name, shots, job.job_id)


class _Fails(_Works):
    def _execute_samples(self, circuit, transpile, shots, job):
        raise RuntimeError("the device is down")

    def _execute_estimator(self, circuit, transpile, shots, observable, job):
        raise RuntimeError("the device is down")


class QdBackEndRoutingTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        patcher = mock.patch.dict(registry._BACKEND_REGISTRY)
        patcher.start()
        self.addCleanup(patcher.stop)
        import Qdislib as qd

        self.qd = qd
        qd.register_backend("WORKS", __name__, "_Works")
        qd.register_backend("FAILS", __name__, "_Fails")

    def _circuit(self, width):
        return types.SimpleNamespace(num_qubits=width)

    def test_a_failure_goes_on_to_the_next_backend_and_empties_its_queue(self):
        unified = self.qd.QdBackEnd(
            {"FAILS": {"num_qubits": 5}, "WORKS": {"num_qubits": 5}}
        )
        result = unified.sample_circuit(self._circuit(2), shots=7)
        self.assertEqual(result.counts, {"0": 7})
        self.assertEqual([b.queue_length for b in unified.backends], [0, 0])

    def test_the_last_failure_is_the_cause_when_all_fail(self):
        unified = self.qd.QdBackEnd({"FAILS": {"num_qubits": 5}})
        with self.assertRaisesRegex(RuntimeError, "No backend could sample") as caught:
            unified.sample_circuit(self._circuit(2))
        self.assertIn("down", str(caught.exception.__cause__))
        with self.assertRaisesRegex(
            RuntimeError, "No backend could estimate"
        ) as caught:
            unified.get_expectation_value(self._circuit(2), "ZZ")
        self.assertIn("down", str(caught.exception.__cause__))

    def test_a_backend_without_estimator_is_skipped(self):
        unified = self.qd.QdBackEnd({"WORKS": {"num_qubits": 5}})
        with self.assertRaisesRegex(RuntimeError, "Last error: None"):
            unified.get_expectation_value(self._circuit(2), "ZZ")

    def test_too_wide_and_unknown_names_are_refused(self):
        unified = self.qd.QdBackEnd({"WORKS": {"num_qubits": 2}})
        with self.assertRaisesRegex(ValueError, "large enough"):
            unified.sample_circuit(self._circuit(3))
        with self.assertRaisesRegex(ValueError, "not found"):
            unified.sample_circuit(self._circuit(1), backend_name="NOPE")

    def test_the_shortest_queue_goes_first(self):
        unified = self.qd.QdBackEnd(
            {
                "WORKS": {"num_qubits": 5, "name": "A"},
                "FAILS": {"num_qubits": 5, "name": "B"},
            }
        )
        unified.backends[0]._enqueue("samples", 1, 1)
        order = unified._select_backends(self._circuit(1), None)
        self.assertEqual([b.name for b in order], ["B", "A"])

    def test_a_width_not_stated_is_the_default(self):
        backend = self.qd.build_backend("WORKS", self.qd.BackendConfig())
        self.assertEqual(backend.num_qubits, registry.DEFAULT_NUM_QUBITS)


class ConverterTest(BaseTimedTestCase):
    def test_roundtrip_and_num_qubits(self):
        from Qdislib.core.backend import CircuitConverter
        from qiskit import QuantumCircuit

        qibo_c = _qibo_circuit()
        qiskit_c = CircuitConverter.to_qiskit(qibo_c)
        self.assertIsInstance(qiskit_c, QuantumCircuit)
        self.assertEqual(CircuitConverter.num_qubits(qibo_c), 2)
        self.assertEqual(CircuitConverter.num_qubits(qiskit_c), 2)
        # (dag, num_qubits) CUDA-Q payload and unknown types.
        self.assertEqual(CircuitConverter.num_qubits((object(), 4)), 4)
        self.assertEqual(CircuitConverter.num_qubits(object()), 0)

    def test_to_qiskit_keeps_each_gate_meaning(self):
        """bsc-wdc/qdislib#13: CU3 (and gates qelib1.inc lacks) convert exactly."""
        import numpy as np
        from qibo import models, gates
        from qiskit.quantum_info import Operator
        from Qdislib.core.backend import CircuitConverter

        for gate in (
            gates.CU3(0, 1, 0.1, 0.9, 0.7),
            gates.SWAP(0, 1),
            gates.iSWAP(0, 1),
            gates.RXX(0, 1, 0.3),
            gates.RYY(0, 1, 0.3),
            gates.CRX(0, 1, 0.3),
            gates.CCZ(0, 1, 2),
        ):
            circuit = models.Circuit(3)
            circuit.add(gates.H(0))
            circuit.add(gates.H(1))
            circuit.add(gate)
            converted = CircuitConverter.to_qiskit(circuit)
            # Qiskit orders qubits little-endian, Qibo big-endian.
            self.assertTrue(
                np.allclose(Operator(converted.reverse_bits()).data, circuit.unitary()),
                gate.__class__.__name__,
            )


@unittest.skipUnless(AER_AVAILABLE and QIBO_AVAILABLE, "qiskit_aer+qibo required")
class AerBackendTest(BaseTimedTestCase):
    def test_default_is_aer(self):
        from Qdislib.core.backend import QdBackEnd

        self.assertEqual(QdBackEnd().backends[0].name, "AER_BACKEND")

    def test_sample_and_estimate(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit

        exact = qd.analytical_solution(_qibo_circuit(), "ZZ")
        qc = QuantumCircuit(2)
        qc.ry(0.8, 0)
        qc.cz(0, 1)
        qc.ry(0.5, 1)
        be = qd.QdBackEnd({"AER_BACKEND": {"num_qubits": 5}})

        est = be.get_expectation_value(qc, "ZZ", shots=20000)
        self.assertLess(abs(est.expectation_value - exact), _TOL)
        smp = be.sample_circuit(qc, shots=4000)
        self.assertEqual(sum(smp.counts.values()), 4000)

    def test_can_run_size_limit(self):
        from Qdislib.core.backend import QdBackEnd
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(3)
        qc.h(0)
        be = QdBackEnd({"AER_BACKEND": {"num_qubits": 2}})  # too small
        with self.assertRaises(ValueError):  # _select_backends: no compatible backend
            be.sample_circuit(qc, shots=100)


class SampleBitOrderTest(BaseTimedTestCase):
    """sample_circuit's keys: qubit 0 first, whichever backend runs it."""

    @unittest.skipUnless(
        QISKIT_AVAILABLE and QIBO_AVAILABLE, "qiskit and qibo required"
    )
    def test_every_simulator(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        qc = QuantumCircuit(2)
        qc.x(0)
        for name in ("AER_BACKEND", "QIBO_BACKEND"):
            with self.subTest(backend=name):
                be = qd.QdBackEnd({name: {"num_qubits": 4}})
                self.assertEqual(be.sample_circuit(qc, shots=10).counts, {"10": 10})
        if CUDAQ_AVAILABLE:
            from Qdislib.utils.backends import convert_to_dag
            from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel

            kernel = dag_to_circuit_cudaq_kernel(convert_to_dag(qc), 2)[0]
            be = qd.QdBackEnd({"CUDAQ_BACKEND": {"num_qubits": 4, "target": "qpp-cpu"}})
            self.assertEqual(
                be.backends[0].run_samples(kernel, False, 10).counts, {"10": 10}
            )

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit required")
    def test_a_circuit_with_its_own_measurements(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        qc = QuantumCircuit(3)
        qc.x(0)
        qc.measure_all()
        be = qd.QdBackEnd({"AER_BACKEND": {"num_qubits": 4}})
        self.assertEqual(be.sample_circuit(qc, shots=10).counts, {"100": 10})

    def test_the_device_backend(self):
        from Qdislib.core.backend.base_backend import BackendConfig, Job
        from Qdislib.core.backend.iqm_backend import IQMBackend

        backend = IQMBackend.__new__(IQMBackend)
        backend.config = BackendConfig(num_qubits=2)
        backend.name = "IQM_BACKEND"
        backend._counts = lambda circuit, transpile, shots: {"01": shots}
        job = types.SimpleNamespace(job_id="j")
        result = backend._execute_samples(None, False, 10, job)
        self.assertEqual(result.counts, {"10": 10})


@unittest.skipUnless(QIBO_AVAILABLE, "qibo required")
class QiboBackendTest(BaseTimedTestCase):
    def test_sample_and_estimate(self):
        import Qdislib as qd

        exact = qd.analytical_solution(_qibo_circuit(), "ZZ")
        be = qd.QdBackEnd({"QIBO_BACKEND": {"num_qubits": 5}})
        est = be.get_expectation_value(_qibo_circuit(), "ZZ", shots=20000)
        self.assertLess(abs(est.expectation_value - exact), _TOL)
        smp = be.sample_circuit(_qibo_circuit(), shots=4000)
        self.assertEqual(sum(smp.counts.values()), 4000)

    def test_sample_leaves_caller_circuit_unchanged(self):
        """bsc-wdc/qdislib#12: sampling must not add M to the caller's circuit."""
        import Qdislib as qd

        circuit = _qibo_circuit()
        be = qd.QdBackEnd({"QIBO_BACKEND": {"num_qubits": 5}})
        be.sample_circuit(circuit, shots=100)
        self.assertEqual([g.name for g in circuit.queue], ["ry", "cz", "ry"])
        self.assertFalse(circuit.measurements)


@unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
class CudaqBackendTest(BaseTimedTestCase):
    def test_sample_and_estimate(self):
        import cudaq
        import Qdislib as qd

        exact = qd.analytical_solution(_qibo_circuit(), "ZZ")
        kernel = cudaq.make_kernel()
        q = [kernel.qalloc() for _ in range(2)]
        kernel.ry(0.8, q[0])
        kernel.cz(q[0], q[1])
        kernel.ry(0.5, q[1])
        be = qd.QdBackEnd({"CUDAQ_BACKEND": {"num_qubits": 5}})
        est = be.get_expectation_value(
            kernel, "ZZ", shots=20000, backend_name="CUDAQ_BACKEND"
        )
        self.assertLess(abs(est.expectation_value - exact), _TOL)
        smp = be.sample_circuit(kernel, shots=4000, backend_name="CUDAQ_BACKEND")
        self.assertEqual(sum(smp.counts.values()), 4000)
