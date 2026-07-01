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

"""Sampling-based gate cutting: reconstruct the output distribution, not just <O>."""

import importlib.util
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
AER_AVAILABLE = importlib.util.find_spec("qiskit_aer") is not None
CUDAQ_AVAILABLE = importlib.util.find_spec("cudaq") is not None

_TOL = 0.06  # total-variation distance vs the exact statevector distribution


def _circuit():
    from qibo import models, gates

    c = models.Circuit(3)
    c.add(gates.RY(0, 0.8))
    c.add(gates.CZ(0, 1))
    c.add(gates.RY(1, 0.5))
    c.add(gates.CZ(1, 2))
    c.add(gates.RY(2, 0.6))
    return c


def _exact_distribution():
    import numpy as np
    import qibo

    qibo.set_backend("numpy")
    probs = np.abs(_circuit()().state()) ** 2
    return {format(b, "03b"): float(probs[b]) for b in range(8)}


def _tvd(exact, reconstructed):
    keys = set(exact) | set(reconstructed)
    return sum(abs(exact.get(k, 0.0) - reconstructed.get(k, 0.0)) for k in keys) / 2


@unittest.skipUnless(QIBO_AVAILABLE, "qibo required")
class GateCuttingSamplingTest(BaseTimedTestCase):
    def _check(self, software, cut, shots=30000):
        import Qdislib as qd

        dist = qd.gate_cutting_sampling(_circuit(), cut, shots=shots, software=software)
        self.assertAlmostEqual(sum(dist.values()), 1.0, places=6)  # valid distribution
        self.assertLess(_tvd(_exact_distribution(), dist), _TOL, f"{software} {cut}")

    def test_qibo_one_cut(self):
        self._check("qibo", ["CZ_2"])

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_one_cut(self):
        self._check("qiskit", ["CZ_2"])

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_two_cuts(self):
        self._check("qiskit", ["CZ_2", "CZ_4"])

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
    def test_cudaq_one_cut(self):
        self._check("cudaq", ["CZ_2"])
