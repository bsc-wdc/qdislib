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

    def test_cut_pmf_normalization(self):
        # The sampling PMF must be a valid distribution whose 1-norm phi matches the
        # closed form 2||c||_1^2 - ||c||_2^2 (and equals gamma+... = phi_factor).
        from Qdislib.core.cutting_algorithms.hadamard_cut import _cut_pmf, phi_factor

        for name in ("CZ", "CP"):
            params = (0.7,) if name == "CP" else ()
            terms = local_unitary_decomposition(name, params)
            settings, probs, phi = _cut_pmf(terms)
            self.assertAlmostEqual(sum(probs), 1.0, places=12)
            self.assertTrue(all(p >= 0 for p in probs))
            self.assertAlmostEqual(phi, phi_factor(terms), places=9)
            self.assertEqual(len(settings), len(probs))

    def test_pauli_fallback_reconstructs(self):
        # The qiskit-free fallback (used when the KAK engine is unavailable) still
        # yields a correct unitary decomposition -- just not gamma-optimal.
        from Qdislib.core.cutting_algorithms.hadamard_cut import _pauli_local_terms

        iswap = np.array(
            [[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]], dtype=complex
        )
        terms = _pauli_local_terms(iswap)
        recon = sum(c * np.kron(V, W) for c, V, W in terms)
        self.assertTrue(np.allclose(recon, iswap))
        for c, V, W in terms:
            self.assertTrue(np.allclose(V.conj().T @ V, np.eye(2)))
            self.assertTrue(np.allclose(W.conj().T @ W, np.eye(2)))
            self.assertGreater(c, 0)

    @unittest.skipUnless(QIBO_AVAILABLE, "qibo builds the 2-qubit matrix")
    def test_iswap_decomposes_via_kak(self):
        # iSWAP is non-diagonal -> no ZZ fast path -> general KAK/Cartan route. It
        # still yields unitary local terms; gamma(iSWAP) = 7 is the known optimum.
        terms = local_unitary_decomposition("iSWAP")
        recon = sum(c * np.kron(V, W) for c, V, W in terms)
        expected = np.array(
            [[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]], dtype=complex
        )
        self.assertTrue(np.allclose(recon, expected))
        for c, V, W in terms:
            self.assertTrue(np.allclose(V.conj().T @ V, np.eye(2)))
            self.assertTrue(np.allclose(W.conj().T @ W, np.eye(2)))
            self.assertGreater(c, 0)
        self.assertAlmostEqual(gamma_factor(terms), 7.0, places=6)

    @unittest.skipUnless(QIBO_AVAILABLE, "qibo builds the 2-qubit matrix")
    def test_unsupported_gate_raises(self):
        # A three-qubit gate cannot be a two-qubit cut -> NotImplementedError.
        with self.assertRaises(NotImplementedError):
            local_unitary_decomposition("TOFFOLI")


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

    def test_iswap_cut(self):
        # A genuinely non-diagonal cut gate, cut via the KAK/Cartan path.
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.9))
        c.add(gates.RY(1, 0.4))
        c.add(gates.iSWAP(0, 1))
        c.add(gates.RX(0, 0.3))
        for obs in ("ZZ", "ZI", "IZ"):
            self._check(c, ["ISWAP_3"], obs)

    def test_swap_cut(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.7))
        c.add(gates.RX(1, 1.1))
        c.add(gates.SWAP(0, 1))
        self._check(c, ["SWAP_3"], "ZZ")

    def test_fsim_cut(self):
        # A parametrised non-diagonal gate; non-terminal on both sides.
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.6))
        c.add(gates.RY(1, 0.9))
        c.add(gates.fSim(0, 1, 0.4, 0.7))
        c.add(gates.RX(0, 0.5))
        c.add(gates.RZ(1, 0.8))
        for obs in ("ZZ", "ZI", "IZ"):
            self._check(c, ["FSIM_3"], obs)

    def test_non_conventional_gates_in_fragments(self):
        # The cut gate stays in the supported ZZ set, but the fragment BODIES carry
        # multi-parameter (U3) and mixed-case (iSWAP) gates -- these must round-trip
        # through the fragment builder / neutral DAG and reconstruct exactly.
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.5))  # RY_1   -> fragment {0}
        c.add(gates.U3(1, 0.4, 0.9, 1.2))  # U3_2   (3 params)
        c.add(gates.U3(2, 0.3, 0.7, 0.2))  # U3_3
        c.add(gates.CZ(0, 1))  # CZ_4   <- the cut (separates {0} | {1,2})
        c.add(gates.iSWAP(1, 2))  # ISWAP_5 (mixed-case class) inside fragment {1,2}
        c.add(gates.RX(2, 0.6))  # RX_6
        for obs in ("ZZZ", "ZIZ"):
            self._check(c, ["CZ_4"], obs)

    def test_sampling_path_matches_within_noise(self):
        # The shot-based reconstruction is unbiased -> close to analytical at high shots.
        import qibo
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import gate_cutting_hadamard

        qibo.set_backend("numpy")
        c = self._bell_like()
        analytical = qd.analytical_solution(c, "ZZ")
        sampled = gate_cutting_hadamard(c, ["CZ_2"], "ZZ", shots=8000, seed=1)
        self.assertLess(abs(analytical - sampled), 0.08)

    def test_multibackend_all_in_one(self):
        # The sampled fan-out routes each fragment through build_execution_backend, so
        # every installed backend reconstructs within shot noise.
        import importlib.util
        import qibo
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import gate_cutting_hadamard

        qibo.set_backend("numpy")
        c = self._bell_like()
        analytical = qd.analytical_solution(c, "ZZ")
        for sw in ("qibo", "qiskit", "pennylane", "cudaq"):
            if importlib.util.find_spec(sw.replace("-", "")) is None:
                continue
            got = gate_cutting_hadamard(c, ["CZ_2"], "ZZ", software=sw, shots=8000)
            self.assertLess(abs(analytical - got), 0.12, sw)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class HadamardThreeStepTest(BaseTimedTestCase):
    """The explicit three-step workflow: subcircuits -> evaluate -> reconstruct."""

    def _circuit(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RY(0, 0.8))
        return c

    def test_subcircuits_exact_reconstruction(self):
        # Evaluating each fragment exactly must reconstruct <observable> to machine eps.
        import qibo
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuits,
            _expectation_z,
        )

        qibo.set_backend("numpy")
        c = self._circuit()
        analytical = qd.analytical_solution(c, "ZZ")
        subs = gate_cutting_hadamard_subcircuits(c, ["CZ_2"], "ZZ", software="qibo")
        self.assertEqual(len(subs), len(subs.config_of))
        values = []
        for frag, obs in subs:
            z_qubits = [q for q, ch in enumerate(obs) if ch == "Z"]
            values.append(_expectation_z(frag, z_qubits))
        self.assertLess(abs(analytical - subs.reconstruct(values)), _TOL)

    def test_subcircuits_via_expectation_value(self):
        # The advertised workflow: expectation_value on each subcircuit, then reconstruct.
        import qibo
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuits,
        )

        qibo.set_backend("numpy")
        c = self._circuit()
        analytical = qd.analytical_solution(c, "ZZ")
        subs = gate_cutting_hadamard_subcircuits(c, ["CZ_2"], "ZZ", software="qibo")
        values = [qd.expectation_value(s, shots=8000) for s in subs]
        self.assertLess(abs(analytical - subs.reconstruct(values)), 0.08)

    def test_standalone_reconstruction_matches_container(self):
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_subcircuit_reconstruction,
        )

        weights = [0.5, -0.25, 0.25]
        config_of = [0, 0, 1, 2]  # two fragments in config 0, one each in 1 and 2
        values = [0.8, 0.6, 0.4, -0.2]
        # config 0: 0.5*(0.8*0.6); config 1: -0.25*0.4; config 2: 0.25*(-0.2)
        expected = 0.5 * (0.8 * 0.6) + (-0.25) * 0.4 + 0.25 * (-0.2)
        got = gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
        self.assertAlmostEqual(got, expected, places=12)

    def test_tree_reduction_matches_single(self):
        # The hierarchical reduction just regroups the same weighted-product sum, so it
        # equals the single fan-in reduction for every chunk size (configs partition
        # cleanly -- a fragment product never spans configurations).
        import numpy as np
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            _reduce_hadamard_tree,
            gate_cutting_hadamard_subcircuit_reconstruction,
        )
        from Qdislib.utils.runtime import compss_wait_on

        rng = np.random.default_rng(0)
        n_configs = 37
        weights = list(rng.normal(size=n_configs))
        values, config_of = [], []
        for cfg in range(n_configs):
            for _ in range(3):  # three fragments per configuration
                values.append(float(rng.uniform(-1, 1)))
                config_of.append(cfg)
        single = gate_cutting_hadamard_subcircuit_reconstruction(
            values, weights, config_of
        )
        for chunk in (1, 2, 5, 8, 37, 100):
            tree = compss_wait_on(
                _reduce_hadamard_tree(values, weights, config_of, chunk)
            )
            self.assertAlmostEqual(tree, single, places=9, msg=f"chunk={chunk}")


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class HadamardSamplingTest(BaseTimedTestCase):
    """The Monte-Carlo (PMF-sampling) driver -- the native Harrow-Lowe estimator."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _bell_like(self):
        from qibo import gates, models

        c = models.Circuit(2)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RY(0, 0.8))
        return c

    def test_sampling_converges_single_cut(self):
        # Config-sampling with exact fragments: only the PMF draws are stochastic, so a
        # large sample count lands within Monte-Carlo noise of the analytical value.
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
        )

        c = self._bell_like()
        analytical = qd.analytical_solution(c, "ZZ")
        got = gate_cutting_hadamard_sampling(
            c, ["CZ_2"], "ZZ", num_samples=8000, shots=None, seed=7
        )
        self.assertLess(abs(analytical - got), 0.1)

    def test_sampling_two_cuts(self):
        # Enumeration would be 36 configs; sampling stays at num_samples.
        import Qdislib as qd
        from qibo import gates, models
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
        )

        c = models.Circuit(3)
        c.add(gates.RY(0, 0.5))
        c.add(gates.RY(1, 0.7))
        c.add(gates.RY(2, 0.9))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.RX(2, 0.3))
        analytical = qd.analytical_solution(c, "ZZZ")
        got = gate_cutting_hadamard_sampling(
            c, ["CZ_4", "CZ_5"], "ZZZ", num_samples=10000, shots=None, seed=3
        )
        self.assertLess(abs(analytical - got), 0.12)

    def test_sampling_non_diagonal_cut(self):
        import Qdislib as qd
        from qibo import gates, models
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
        )

        c = models.Circuit(2)
        c.add(gates.RY(0, 0.9))
        c.add(gates.RY(1, 0.4))
        c.add(gates.iSWAP(0, 1))
        analytical = qd.analytical_solution(c, "ZZ")
        got = gate_cutting_hadamard_sampling(
            c, ["ISWAP_3"], "ZZ", num_samples=8000, shots=None, seed=5
        )
        self.assertLess(abs(analytical - got), 0.12)

    def test_sampling_subcircuits_workflow(self):
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling_subcircuits,
        )

        c = self._bell_like()
        analytical = qd.analytical_solution(c, "ZZ")
        subs = gate_cutting_hadamard_sampling_subcircuits(
            c, ["CZ_2"], "ZZ", software="qibo", num_samples=8000, seed=1
        )
        # exact per-fragment eval isolates the config-sampling noise
        from Qdislib.core.cutting_algorithms.hadamard_cut import _expectation_z

        values = [
            _expectation_z(frag, [q for q, ch in enumerate(obs) if ch == "Z"])
            for frag, obs in subs
        ]
        self.assertLess(abs(analytical - subs.reconstruct(values)), 0.1)

    def test_seed_reproducible(self):
        # The seed fixes the PMF draws (the stochastic part); the reconstructed value
        # matches to ~machine eps (only BLAS-level float noise in the simulation, not
        # the sampling, differs between runs).
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
            _sample_configs,
            _hadamard_setup,
        )

        c = self._bell_like()
        _, _, decomps = _hadamard_setup(c, ["CZ_2"])
        cfg_a, _ = _sample_configs(decomps, 2000, 42)
        cfg_b, _ = _sample_configs(decomps, 2000, 42)
        self.assertEqual(cfg_a, cfg_b)  # the draws themselves are reproducible

        a = gate_cutting_hadamard_sampling(
            c, ["CZ_2"], "ZZ", num_samples=2000, shots=None, seed=42
        )
        b = gate_cutting_hadamard_sampling(
            c, ["CZ_2"], "ZZ", num_samples=2000, shots=None, seed=42
        )
        self.assertAlmostEqual(a, b, places=10)

    def test_tree_reduction_end_to_end(self):
        # reduce_chunk (hierarchical fan-in) must not change the estimate: same seed +
        # exact fragments -> identical draws and values, only the reduction is regrouped.
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            gate_cutting_hadamard_sampling,
        )

        c = self._bell_like()
        single = gate_cutting_hadamard_sampling(
            c, ["CZ_2"], "ZZ", num_samples=3000, shots=None, seed=9
        )
        tree = gate_cutting_hadamard_sampling(
            c, ["CZ_2"], "ZZ", num_samples=3000, shots=None, seed=9, reduce_chunk=16
        )
        self.assertAlmostEqual(single, tree, places=9)


if __name__ == "__main__":
    unittest.main()
