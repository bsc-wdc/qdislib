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

"""Packing inside the IBM backend and the engine, on a simulated device."""

import importlib.util
import pickle
import types
import unittest
from unittest import mock

from tests import BaseTimedTestCase

AER = importlib.util.find_spec("qiskit_aer") is not None
RUNTIME = importlib.util.find_spec("qiskit_ibm_runtime") is not None


class _Sampler:
    """SamplerV2 on a noiseless simulator, counting its jobs (plain ones: it takes
    ``mode``, and each ``run`` holds one circuit)."""

    jobs = []

    def __init__(self, mode=None, options=None):
        self.mode = mode

    def run(self, circuits, shots=None):
        from qiskit_aer.primitives import SamplerV2

        assert len(circuits) == 1, "one circuit per job"
        _Sampler.jobs.append(circuits[0].num_qubits)
        sampler = SamplerV2(
            options={"backend_options": {"method": "matrix_product_state"}}
        )
        return sampler.run(circuits, shots=shots)


def _device():
    from qiskit.providers.fake_provider import GenericBackendV2
    from qiskit.transpiler import CouplingMap

    cmap = CouplingMap.from_heavy_hex(5)
    return GenericBackendV2(
        cmap.size(),
        basis_gates=["cz", "rz", "sx", "x", "id"],
        coupling_map=cmap,
        seed=11,
    )


def _circuit():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(-0.7, 2)
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    return qc


def _exact(qc):
    from qiskit.quantum_info import Statevector

    return {k[::-1]: v for k, v in Statevector(qc).probabilities_dict().items()}


class DefaultTest(BaseTimedTestCase):
    def test_a_list_runs_one_after_another(self):
        from Qdislib.core.backend.base_backend import BackendConfig, BaseBackend

        calls = []

        class One(BaseBackend):
            def evaluate_joint_counts(self, dag, n, cuts=(), shots=1024, seed=None):
                calls.append(dag)
                return {dag: shots}

        out = One("ONE", BackendConfig(num_qubits=2)).evaluate_joint_counts_many(
            [("a", 2, (), 5, None), ("b", 2, (), 7, 1)]
        )
        self.assertEqual(out, [{"a": 5}, {"b": 7}])
        self.assertEqual(calls, ["a", "b"])
        self.assertFalse(One("ONE", BackendConfig(num_qubits=2)).packs)


@unittest.skipUnless(AER and RUNTIME, "qiskit-aer and qiskit-ibm-runtime required")
class IBMPackingTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from Qdislib.core.backend.base_backend import BackendConfig
        from Qdislib.core.backend.qpu_backends import IBMBackend
        from Qdislib.core.cutting_algorithms import gate_cutting as gc

        self.qd = qd
        device = _device()
        self.ibm = IBMBackend(
            BackendConfig(
                num_qubits=device.num_qubits,
                extra={"batch": device, "backend": device, "pack": {"seeds": 1}},
            )
        )
        _Sampler.jobs = []
        patcher = mock.patch.object(gc, "Sampler", _Sampler)
        patcher.start()
        self.addCleanup(patcher.stop)

    def assertClose(self, got, want, delta=0.06):
        for key in set(got) | set(want):
            self.assertAlmostEqual(
                got.get(key, 0.0), want.get(key, 0.0), delta=delta, msg=key
            )

    def test_distribution_in_fewer_jobs_than_subcircuits(self):
        qc = _circuit()
        for name, cuts, size in (
            ("hadamard", {"hadamard": ["CZ_4"]}, 12),
            ("wire", {"wire": [("CZ_4", "RY_5")]}, 8),
            ("gate", {"gate": ["CZ_4"]}, 12),
        ):
            with self.subTest(cut=name):
                _Sampler.jobs = []
                subs = self.qd.circuit_cutting_sampling_subcircuits(qc, cuts)
                got = self.qd.circuit_cutting_sampling(
                    qc,
                    cuts,
                    backend=self.ibm,
                    shots=40000,
                    seed=1,
                    circuits_per_task=size,
                )
                self.assertClose(got, _exact(qc))
                self.assertLess(len(_Sampler.jobs), len(subs))
                self.assertEqual(set(_Sampler.jobs), {_device().num_qubits})

    def test_expectation_value(self):
        from qiskit.quantum_info import SparsePauliOp, Statevector

        qc = _circuit()
        want = Statevector(qc).expectation_value(SparsePauliOp("ZZZZ")).real
        got = self.qd.circuit_cutting(
            qc,
            {"hadamard": ["CZ_4"]},
            "ZZZZ",
            backend=self.ibm,
            shots=40000,
            seed=2,
            circuits_per_task=12,
        )
        self.assertAlmostEqual(got, want, delta=0.1)

    def test_with_error_bars(self):
        got = self.qd.circuit_cutting(
            _circuit(),
            {"hadamard": ["CZ_4"]},
            "ZZZZ",
            backend=self.ibm,
            shots=20000,
            seed=3,
            circuits_per_task=12,
            error=5,
        )
        self.assertGreater(got.error, 0)

    def test_without_groups_it_says_so(self):
        import warnings
        from Qdislib.utils import logging_config

        logging_config._WARNED.clear()
        self.addCleanup(logging_config._WARNED.clear)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.qd.circuit_cutting_sampling(
                _circuit(), {"hadamard": ["CZ_4"]}, backend=self.ibm, shots=100, seed=1
            )
        self.assertTrue(any("circuits_per_task" in str(w.message) for w in caught))
        self.assertTrue(self.ibm.packs)


@unittest.skipUnless(RUNTIME, "qiskit-ibm-runtime required")
class IBMGoesToATaskTest(BaseTimedTestCase):
    SECRET = "not-a-real-token-1234"

    def _backend(self, **extra):
        from Qdislib.core.backend.base_backend import BackendConfig
        from Qdislib.core.backend.qpu_backends import IBMBackend

        device = types.SimpleNamespace(name="ibm_fake")
        service = types.SimpleNamespace(backend=lambda name: device)
        return IBMBackend(
            BackendConfig(
                num_qubits=5,
                extra={"service": service, "backend_name": "ibm_fake", **extra},
            )
        )

    def test_no_token_or_connection_travels(self):
        blob = pickle.dumps(self._backend(token=self.SECRET, channel="ibm_cloud"))
        self.assertNotIn(self.SECRET.encode(), blob)
        copy = pickle.loads(blob)
        self.assertNotIn("token", copy.config.extra)
        self.assertNotIn("service", copy.config.extra)
        self.assertEqual(copy.config.extra["backend_name"], "ibm_fake")
        self.assertIsNone(copy._ibm_backend)

    def test_a_copy_connects_once_when_used(self):
        from Qdislib.core.backend import qpu_backends as qb

        copy = pickle.loads(pickle.dumps(self._backend()))
        other = pickle.loads(pickle.dumps(self._backend()))
        device = types.SimpleNamespace(name="ibm_fake")
        qb._CONNECTIONS.clear()
        self.addCleanup(qb._CONNECTIONS.clear)
        with mock.patch.object(
            qb.IBMBackend, "_connect", return_value=(device, device)
        ) as connect:
            self.assertIs(copy.ibm_backend, device)
            self.assertIs(copy.session, device)
            self.assertIs(other.ibm_backend, device)
        self.assertEqual(connect.call_count, 1)

    def test_an_opened_batch_cannot_go(self):
        from Qdislib.core.backend.base_backend import BackendConfig
        from Qdislib.core.backend.qpu_backends import IBMBackend

        device = types.SimpleNamespace(name="ibm_fake")
        ibm = IBMBackend(
            BackendConfig(num_qubits=5, extra={"batch": object(), "backend": device})
        )
        with self.assertRaisesRegex(TypeError, "Batch or Session"):
            pickle.dumps(ibm)


if __name__ == "__main__":
    unittest.main()
