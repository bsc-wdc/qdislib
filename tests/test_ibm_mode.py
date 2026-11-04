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

"""The IBM backend submits plain jobs unless a Batch or Session is asked for."""

import importlib.util
import unittest
from unittest import mock

from tests import BaseTimedTestCase

RUNTIME = importlib.util.find_spec("qiskit_ibm_runtime") is not None


class _Service:
    """A stand-in QiskitRuntimeService: no network, no credentials."""

    def __init__(self):
        self.device = object()

    def backend(self, name):
        return self.device

    def least_busy(self, **kwargs):
        return self.device


@unittest.skipUnless(RUNTIME, "qiskit-ibm-runtime required")
class IBMModeTest(BaseTimedTestCase):
    def _connect(self, **extra):
        from Qdislib.core.backend.qpu_backends import IBMBackend

        service = _Service()
        return service, IBMBackend._connect({"service": service, **extra})

    def test_default_is_plain_jobs(self):
        with (
            mock.patch("qiskit_ibm_runtime.Batch") as batch,
            mock.patch("qiskit_ibm_runtime.Session") as session,
        ):
            service, (mode, backend) = self._connect(backend_name="ibm_x")
        self.assertIs(mode, service.device)
        self.assertIs(backend, service.device)
        batch.assert_not_called()
        session.assert_not_called()

    def test_batch_and_session_on_request(self):
        for name in ("batch", "session"):
            with (
                self.subTest(mode=name),
                mock.patch("qiskit_ibm_runtime.Batch") as batch,
                mock.patch("qiskit_ibm_runtime.Session") as session,
            ):
                service, (mode, _) = self._connect(mode=name, max_time=30)
                opened = batch if name == "batch" else session
                other = session if name == "batch" else batch
                opened.assert_called_once_with(backend=service.device, max_time=30)
                other.assert_not_called()
                self.assertIs(mode, opened.return_value)

    def test_bad_mode(self):
        with self.assertRaises(ValueError):
            self._connect(mode="pipeline")

    def test_check_ibm_qc_is_plain_by_default(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import _check_ibm_qc

        service = _Service()
        with mock.patch("qiskit_ibm_runtime.Batch") as batch:
            mode, backend = _check_ibm_qc({"service": service})
        self.assertIs(mode, service.device)
        self.assertIs(backend, service.device)
        batch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
