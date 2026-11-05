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

"""A backend on a GPU runs its circuits in tasks that ask COMPSs for a GPU."""

import importlib.util
import types
import unittest
from unittest import mock

from tests import BaseTimedTestCase

AER = importlib.util.find_spec("qiskit_aer") is not None


def _circuit():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(0.7, 2)
    qc.cx(2, 3)
    return qc


class SelectionTest(BaseTimedTestCase):
    def test_uses_gpu_of_the_backends(self):
        from Qdislib.core.backend.base_backend import BaseBackend
        from Qdislib.core.backend.cudaq_backend import CudaqBackend

        self.assertFalse(BaseBackend.uses_gpu.fget(types.SimpleNamespace()))
        for target, expected in (
            ("nvidia", True),
            ("nvidia-mgpu", True),
            ("qpp-cpu", False),
            (None, False),
        ):
            fake = types.SimpleNamespace(_target=target)
            self.assertEqual(CudaqBackend.uses_gpu.fget(fake), expected, target)

    def test_task_follows_the_backend(self):
        from Qdislib.core.cutting_algorithms import mixed_cutting as mc

        for runner in ("qiskit", "cudaq", None, types.SimpleNamespace(uses_gpu=False)):
            self.assertIs(mc._fragment_task(runner), mc._run_fragment)
            self.assertIs(mc._fragments_task(runner), mc._run_fragments)
        on_gpu = types.SimpleNamespace(uses_gpu=True)
        self.assertIs(mc._fragment_task(on_gpu), mc._run_fragment_gpu)
        self.assertIs(mc._fragments_task(on_gpu), mc._run_fragments_gpu)


@unittest.skipUnless(AER, "qiskit-aer required")
class PickleTest(BaseTimedTestCase):
    def test_a_backend_can_go_to_a_task(self):
        import pickle
        import Qdislib as qd
        from Qdislib.core.backend.aer_backend import AerBackend
        from Qdislib.core.backend.base_backend import BackendConfig

        backend = AerBackend(BackendConfig(num_qubits=4))
        copy = pickle.loads(pickle.dumps(backend))
        self.assertEqual(copy.name, backend.name)
        self.assertEqual(copy.queue_length, 0)
        got = qd.circuit_cutting_sampling(
            _circuit(), {"hadamard": ["CZ_3"]}, backend=copy, shots=500, seed=1
        )
        self.assertAlmostEqual(sum(got.values()), 1.0, places=6)


@unittest.skipUnless(AER, "qiskit-aer required")
class GpuRunTest(BaseTimedTestCase):
    """The GPU tasks give the CPU tasks' results (an Aer backend on the CPU that only
    says it uses a GPU: this runs where no GPU is)."""

    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from Qdislib.core.backend.aer_backend import AerBackend
        from Qdislib.core.backend.base_backend import BackendConfig

        class SaysGpu(AerBackend):
            uses_gpu = True

        self.qd, self.backend = qd, SaysGpu(BackendConfig(num_qubits=4))

    def _spy(self, name):
        from Qdislib.core.cutting_algorithms import mixed_cutting as mc

        return mock.patch.object(mc, name, wraps=getattr(mc, name))

    def test_distribution_goes_through_the_gpu_task(self):
        qd, qc = self.qd, _circuit()
        cuts = {"hadamard": ["CZ_3"]}
        want = qd.circuit_cutting_sampling(qc, cuts, shots=4000, seed=3)
        with self._spy("_run_fragment_gpu") as spy:
            got = qd.circuit_cutting_sampling(
                qc, cuts, backend=self.backend, shots=4000, seed=3
            )
        self.assertGreater(spy.call_count, 0)
        for key in set(want) | set(got):
            self.assertAlmostEqual(got.get(key, 0.0), want.get(key, 0.0), delta=0.08)

    def test_groups_go_through_the_gpu_task(self):
        qd, qc = self.qd, _circuit()
        with self._spy("_run_fragments_gpu") as spy:
            qd.circuit_cutting_sampling(
                qc,
                {"hadamard": ["CZ_3"]},
                backend=self.backend,
                shots=200,
                circuits_per_task=4,
            )
        self.assertGreater(spy.call_count, 0)

    def test_a_cpu_backend_stays_on_the_cpu_task(self):
        from Qdislib.core.backend.aer_backend import AerBackend
        from Qdislib.core.backend.base_backend import BackendConfig

        qd, qc = self.qd, _circuit()
        with self._spy("_run_fragment_gpu") as spy:
            qd.circuit_cutting_sampling(
                qc,
                {"hadamard": ["CZ_3"]},
                backend=AerBackend(BackendConfig(num_qubits=4)),
                shots=200,
            )
        spy.assert_not_called()


if __name__ == "__main__":
    unittest.main()
