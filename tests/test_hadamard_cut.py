#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Tests for the ancilla / double-Hadamard-test gate cut (arXiv:2403.01018)."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None

from Qdislib.core.cutting_algorithms.hadamard_cut import (  # noqa: E402
    gamma_factor,
    local_unitary_decomposition,
)

_TOL = 1e-9


class DecompositionTest(BaseTimedTestCase):
    """The local-unitary decomposition core (no qibo needed)."""

    def _reconstruct(self, name, params=()):
        import cmath

        terms = local_unitary_decomposition(name, params)
        recon = sum(c * np.kron(V, W) for c, V, W in terms)
        for c, V, W in terms:  # each factor is unitary, c > 0
            self.assertTrue(np.allclose(V.conj().T @ V, np.eye(2)))
            self.assertTrue(np.allclose(W.conj().T @ W, np.eye(2)))
            self.assertGreater(c, 0)
        return terms, recon

    def test_cz_optimal_gamma(self):
        terms, recon = self._reconstruct("CZ")
        target = np.diag([1, 1, 1, -1]).astype(complex)
        self.assertTrue(np.allclose(recon, target))  # exact reconstruction
        self.assertAlmostEqual(gamma_factor(terms), 3.0, places=6)  # optimal

    def test_cx_reduces_to_cz(self):
        terms, _ = self._reconstruct("CX")
        self.assertAlmostEqual(gamma_factor(terms), 3.0, places=6)

    def test_weak_interaction_low_overhead(self):
        # gamma -> 1 as the interaction weakens (the clustered-Hamiltonian regime).
        terms, _ = self._reconstruct("CP", (0.2,))
        self.assertLess(gamma_factor(terms), 1.2)

    def test_unsupported_gate_raises(self):
        with self.assertRaises(NotImplementedError):
            local_unitary_decomposition("ISWAP")


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class HadamardCutTest(BaseTimedTestCase):
    """Exact reconstruction matches the analytical expectation value."""

    def _check(self, circuit, cuts, observable):
        import qibo
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import gate_cutting_hadamard

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(circuit, observable)
        got = gate_cutting_hadamard(circuit, cuts, observable)
        self.assertLess(
            abs(analytical - got), _TOL, f"{observable}: {analytical} vs {got}"
        )

    def _bell_like(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RY(0, 0.8))
        return c

    def test_cz_cut_observables(self):
        c = self._bell_like()
        for obs in ("ZZ", "ZI", "IZ"):
            self._check(c, ["CZ_2"], obs)

    def test_cx_cut(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CNOT(0, 1))
        c.add(gates.RY(1, 0.5))
        self._check(c, ["CNOT_2"], "ZZ")

    def test_controlled_phase_cut(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.9))
        c.add(gates.RY(1, 0.4))
        c.add(gates.CU1(0, 1, 0.7))
        self._check(c, ["CU1_3"], "ZZ")

    def test_rzz_cut(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.9))
        c.add(gates.RY(1, 0.4))
        c.add(gates.RZZ(0, 1, 1.3))
        self._check(c, ["RZZ_3"], "ZZ")

    def test_non_terminal_cut(self):
        # gates on the cut qubit AFTER the cut (the collapse-and-continue case).
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RX(0, 0.9))
        c.add(gates.RZ(1, 0.4))
        self._check(c, ["CZ_2"], "ZZ")

    def test_two_cuts_three_qubits(self):
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.7))
        c.add(gates.RY(1, 0.9))
        c.add(gates.RY(2, 1.1))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.RX(2, 0.3))
        for obs in ("ZZZ", "ZIZ"):
            self._check(c, ["CZ_4", "CZ_5"], obs)

    def test_mixed_cx_cp_cuts(self):
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.5))
        c.add(gates.RY(1, 0.6))
        c.add(gates.RY(2, 0.7))
        c.add(gates.CNOT(0, 1))
        c.add(gates.CU1(1, 2, 0.8))
        self._check(c, ["CNOT_4", "CU1_5"], "ZZZ")


if __name__ == "__main__":
    unittest.main()
