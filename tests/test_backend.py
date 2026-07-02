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
import unittest

from tests import BaseTimedTestCase

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
