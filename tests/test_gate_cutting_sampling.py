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

# Total-variation distance vs the exact statevector distribution. Kept loose
# because tests run at only 1024 shots (Qibo's collapse-based sampling is slow at
# high shot counts) -- the systematic pre-fix errors were ~0.13-0.22, far above
# this, so shot noise never masks a regression.
_TOL = 0.12


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


def _entangled_circuit():
    """Two internally-entangled 2-qubit blocks joined by two inter-block CZ cuts,
    with gates AFTER the cut on cut qubits. This is the case a naive reconstruction
    (measure-and-force + no cut-point sign) gets systematically wrong: each cut
    qubit keeps evolving after the cut and shares an entangled block with the other
    cut's qubit. Guards the collapse-and-continue reconstruction.
    """
    from qibo import models, gates

    c = models.Circuit(4)
    for qubit, angle in zip(range(4), (0.8, 0.5, 0.6, 0.9)):
        c.add(gates.RY(qubit, angle))
    c.add(gates.CNOT(0, 1))  # entangle block {0, 1}
    c.add(gates.CNOT(2, 3))  # entangle block {2, 3}
    c.add(gates.CZ(0, 2))  # inter-block cut 1  (CZ_7)
    c.add(gates.CZ(1, 3))  # inter-block cut 2  (CZ_8)
    c.add(gates.RY(0, 0.7))  # gate after the cut on a cut qubit (non-terminal)
    c.add(gates.RY(3, 0.7))
    return c


def _entangled_exact():
    import numpy as np
    import qibo

    qibo.set_backend("numpy")
    probs = np.abs(_entangled_circuit()().state()) ** 2
    return {format(b, "04b"): float(probs[b]) for b in range(16)}


def _tvd(exact, reconstructed):
    keys = set(exact) | set(reconstructed)
    return sum(abs(exact.get(k, 0.0) - reconstructed.get(k, 0.0)) for k in keys) / 2


@unittest.skipUnless(QIBO_AVAILABLE, "qibo required")
class GateCuttingSamplingTest(BaseTimedTestCase):
    def _check(self, software, cut, shots=1024):
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

    def _check_three_step(self, software, cut, shots=1024):
        """The three-step API (generate -> sample -> reconstruct) must reconstruct
        the same distribution as the all-in-one gate_cutting_sampling."""
        import Qdislib as qd

        subs = qd.gate_cutting_sampling_subcircuits(_circuit(), cut, software=software)
        sampled = [qd.sample_subcircuit(s, shots=shots) for s in subs]
        dist = subs.reconstruct(sampled)
        self.assertAlmostEqual(sum(dist.values()), 1.0, places=6)
        self.assertLess(_tvd(_exact_distribution(), dist), _TOL, f"3-step {software} {cut}")

    def test_qibo_three_step(self):
        self._check_three_step("qibo", ["CZ_2"])

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_three_step(self):
        self._check_three_step("qiskit", ["CZ_2"])

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
    def test_cudaq_three_step(self):
        self._check_three_step("cudaq", ["CZ_2"])

    def _check_entangled(self, software, shots=1024):
        import Qdislib as qd

        dist = qd.gate_cutting_sampling(
            _entangled_circuit(), ["CZ_7", "CZ_8"], shots=shots, software=software
        )
        self.assertAlmostEqual(sum(dist.values()), 1.0, places=6)
        self.assertLess(
            _tvd(_entangled_exact(), dist), _TOL, f"entangled 2-cut {software}"
        )

    def test_qibo_entangled_two_cuts(self):
        # Regression: cut qubits keep evolving after the cut inside entangled
        # blocks -- was ~0.22 TVD before the collapse-and-continue fix.
        self._check_entangled("qibo")

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_entangled_two_cuts(self):
        self._check_entangled("qiskit")

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
    def test_cudaq_entangled_two_cuts(self):
        # CUDA-Q executes the gate-cut measure path via cudaq.run (which collapses
        # mid-circuit measurements); cudaq.sample does not, which silently biased
        # this case before the decorator/run refactor.
        self._check_entangled("cudaq")

    def _check_wire(self, software, shots=1024):
        """wire_cutting_sampling reconstructs the output distribution (Perlin et al.
        P(x) = 2^-N * sum_A (-1)^(A.x) <Z^A>) for a wire (space-like) cut."""
        import Qdislib as qd

        cut = qd.find_cut(_circuit(), gate_cut=False, wire_cut=True)
        dist = qd.wire_cutting_sampling(_circuit(), cut, shots=shots, software=software)
        self.assertAlmostEqual(sum(dist.values()), 1.0, places=6)
        self.assertLess(_tvd(_exact_distribution(), dist), _TOL, f"wire {software} {cut}")

    def test_qibo_wire(self):
        self._check_wire("qibo")

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_wire(self):
        self._check_wire("qiskit")

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
    def test_cudaq_wire(self):
        self._check_wire("cudaq")

    def _check_wire_three_step(self, software, shots=1024):
        """The wire three-step API (generate -> sample -> reconstruct) must
        reconstruct the same distribution as the all-in-one wire_cutting_sampling."""
        import Qdislib as qd

        cut = qd.find_cut(_circuit(), gate_cut=False, wire_cut=True)
        subs = qd.wire_cutting_sampling_subcircuits(_circuit(), cut, software=software)
        dist = subs.reconstruct([qd.sample_subcircuit(s, shots=shots) for s in subs])
        self.assertAlmostEqual(sum(dist.values()), 1.0, places=6)
        self.assertLess(
            _tvd(_exact_distribution(), dist), _TOL, f"wire 3-step {software} {cut}"
        )

    def test_qibo_wire_three_step(self):
        self._check_wire_three_step("qibo")

    @unittest.skipUnless(AER_AVAILABLE, "qiskit_aer required")
    def test_qiskit_wire_three_step(self):
        self._check_wire_three_step("qiskit")

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq required")
    def test_cudaq_wire_three_step(self):
        self._check_wire_three_step("cudaq")
