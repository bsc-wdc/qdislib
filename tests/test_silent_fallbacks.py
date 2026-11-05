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

"""When a run does something other than what was asked, it says so."""

import importlib.util
import sys
import types
import unittest
import warnings
from unittest import mock

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None


def _fake_cudaq(fails=("nvidia",)):
    """A stand-in ``cudaq`` module recording the targets it was set to."""
    fake = types.ModuleType("cudaq")
    fake.targets = []

    def set_target(name):
        if name in fails:
            raise RuntimeError(f"{name} not available")
        fake.targets.append(name)

    fake.set_target = set_target
    return fake


class CudaqTargetTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        from Qdislib.utils import logging_config

        logging_config._WARNED.clear()
        self.addCleanup(logging_config._WARNED.clear)

    def _select(self, fake, target=None):
        from Qdislib.core.backend.cudaq_backend import select_cudaq_target

        with mock.patch.dict(sys.modules, {"cudaq": fake}):
            select_cudaq_target(target)

    def test_no_gpu_warns_once_with_the_reason(self):
        fake = _fake_cudaq()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._select(fake)
            self._select(fake)
        self.assertEqual(fake.targets, ["qpp-cpu", "qpp-cpu"])
        self.assertEqual(len(caught), 1)
        self.assertIn("nvidia not available", str(caught[0].message))

    def test_gpu_available_says_nothing(self):
        fake = _fake_cudaq(fails=())
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._select(fake)
        self.assertEqual(fake.targets, ["nvidia"])
        self.assertEqual(caught, [])

    def test_explicit_target_is_not_replaced(self):
        fake = _fake_cudaq()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._select(fake, "qpp-cpu")
        self.assertEqual(fake.targets, ["qpp-cpu"])
        self.assertEqual(caught, [])
        with self.assertRaises(
            RuntimeError
        ):  # an explicit target that fails is an error
            self._select(fake, "nvidia")


@unittest.skipUnless(QISKIT, "qiskit required")
class FallbackWarningsTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd
        from qiskit import QuantumCircuit
        from Qdislib.utils import logging_config

        logging_config._WARNED.clear()
        self.addCleanup(logging_config._WARNED.clear)
        self.qd = qd
        qc = QuantumCircuit(3)
        qc.ry(0.4, 0)
        qc.cx(0, 1)
        qc.ry(0.7, 1)
        qc.cx(1, 2)
        self.qc = qc

    def test_kak_failure_warns_and_stays_correct(self):
        import numpy as np
        from Qdislib.core.cutting_algorithms import hadamard_cut as hc

        u = np.diag([1, 1, 1, -1]).astype(complex)
        with mock.patch.object(
            hc, "_kak_local_terms", side_effect=RuntimeError("no kak")
        ):
            with self.assertWarnsRegex(UserWarning, "no kak.*larger"):
                terms = hc._general_local_terms(u)
        rebuilt = sum(c * np.kron(v, w) for c, v, w in terms)
        self.assertTrue(np.allclose(rebuilt, u))

    def test_optimized_wire_failure_warns_and_stays_correct(self):
        # the 1.x path, still behind qpu=, gpu=, method= and a custom CacheStrategy
        from qiskit.quantum_info import SparsePauliOp, Statevector
        from Qdislib.core.cutting_algorithms import wire_cutting as wc
        from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

        exact = Statevector(self.qc).expectation_value(SparsePauliOp("ZZZ")).real
        dag = circuit_qiskit_to_dag(self.qc)
        with mock.patch.object(
            wc, "_wire_cutting_optimized", side_effect=RuntimeError("boom")
        ):
            with self.assertWarnsRegex(UserWarning, "boom.*standard"):
                got = wc._wire_cutting(
                    dag,
                    [("CX_2", "RY_3")],
                    observables="ZZZ",
                    shots=20000,
                    software="qiskit",
                    optimized=True,
                )
        self.assertAlmostEqual(got, exact, delta=0.1)

    def test_exact_run_says_the_backend_does_not_run(self):
        from Qdislib.core.cutting_algorithms.mixed_cutting import _exact_runner

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.assertEqual(_exact_runner("qiskit"), "qiskit")
        self.assertEqual(caught, [])
        device = types.SimpleNamespace(name="IBM_BACKEND")
        with self.assertWarnsRegex(UserWarning, "IBM_BACKEND"):
            self.assertEqual(_exact_runner(device), "qiskit")


if __name__ == "__main__":
    unittest.main()
