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

"""Observable handling: identity, X/Y basis rotation, and weighted Pauli sums."""

import importlib.util
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None

_TOL = 0.2


def _circuit():
    """A 3-qubit circuit used across the observable tests."""
    from qibo import models, gates

    c = models.Circuit(3)
    c.add(gates.H(0))
    c.add(gates.CZ(0, 1))
    c.add(gates.RY(0, 0.8))
    c.add(gates.CZ(1, 2))
    c.add(gates.RX(2, 0.5))
    return c


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class ObservableTest(BaseTimedTestCase):
    """End-to-end checks that reconstructed expectation values match analytics."""

    def _check(self, observables, software, tol=_TOL):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_circuit(), observables)
        gate = qd.gate_cutting(
            _circuit(),
            ["CZ_2"],
            observables=observables,
            shots=16000,
            software=software,
        )
        wire = qd.wire_cutting(
            _circuit(),
            [("CZ_2", "RY_3")],
            observables=observables,
            shots=16000,
            software=software,
        )
        self.assertLess(abs(gate - analytical), tol, f"gate {observables}/{software}")
        self.assertLess(abs(wire - analytical), tol, f"wire {observables}/{software}")

    def test_identity_observable_qibo(self):
        # Identity on any position, including the trailing qubit (regression).
        for obs in ["ZIZ", "ZZI", "IZZ", "ZII"]:
            self._check(obs, "qibo")

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_identity_observable_qiskit(self):
        # Trailing-I previously failed on Qiskit (little-endian get_counts bug).
        for obs in ["ZIZ", "ZZI", "IZZ", "ZII"]:
            self._check(obs, "qiskit")

    def test_xy_observable_qibo(self):
        for obs in ["XZI", "ZYI", "XYZ", "IXZ"]:
            self._check(obs, "qibo")

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_xy_observable_qiskit(self):
        for obs in ["XZI", "ZYI", "XYZ"]:
            self._check(obs, "qiskit")

    def test_weighted_sum_list_qibo(self):
        self._check([(0.5, "ZZI"), (1.0, "XZI"), (-0.3, "IZZ")], "qibo")

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_weighted_sum_dict_qiskit(self):
        self._check({"ZZI": 0.5, "XZI": 1.0, "IZZ": -0.3}, "qiskit")

    def test_weighted_sum_matches_term_by_term(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        terms = [(0.5, "ZZI"), (1.0, "XZI"), (-0.3, "IZZ")]
        combined = qd.analytical_solution(_circuit(), terms)
        manual = sum(c * qd.analytical_solution(_circuit(), p) for c, p in terms)
        self.assertAlmostEqual(combined, manual, places=6)


class NormalizeObservablesTest(BaseTimedTestCase):
    """Unit tests for the observable-specification normaliser."""

    def test_string(self):
        from Qdislib.utils.backends import normalize_observables

        self.assertEqual(normalize_observables("ZZI"), [(1.0, "ZZI")])

    def test_dict_and_list(self):
        from Qdislib.utils.backends import normalize_observables

        self.assertEqual(
            normalize_observables({"ZZ": 0.5, "XI": -1.0}),
            [(0.5, "ZZ"), (-1.0, "XI")],
        )
        self.assertEqual(
            normalize_observables([(0.5, "ZZ"), (-1.0, "XI")]),
            [(0.5, "ZZ"), (-1.0, "XI")],
        )

    def test_invalid_pauli_character(self):
        from Qdislib.utils.backends import normalize_observables

        with self.assertRaises(ValueError):
            normalize_observables("ZQZ")

    def test_inconsistent_lengths(self):
        from Qdislib.utils.backends import normalize_observables

        with self.assertRaises(ValueError):
            normalize_observables([(1.0, "ZZ"), (1.0, "Z")])

    def test_empty(self):
        from Qdislib.utils.backends import normalize_observables

        with self.assertRaises(ValueError):
            normalize_observables([])
