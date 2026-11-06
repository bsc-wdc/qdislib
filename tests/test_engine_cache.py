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

"""The cutting engine's result cache: what it reuses, what it keeps apart, and
how it is turned on."""

import importlib.util
import tempfile
import unittest
import warnings
from unittest import mock

from tests import BaseTimedTestCase

from Qdislib.core.cutting_algorithms import engine_cache

QISKIT = importlib.util.find_spec("qiskit") is not None
PYZX = importlib.util.find_spec("pyzx") is not None
LMDB = importlib.util.find_spec("lmdb") is not None


def _exact(circuit, pauli):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(pauli[::-1])).real
    )


def _circuit(extra_identity=False):
    """Six qubits in a CZ line, cut at CZ_9 (between qubits 3 and 4).
    ``extra_identity`` adds an X X pair on qubit 0: a different circuit with the
    same fragments up to equivalence."""
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(6)
    for q in range(6):
        qc.ry(0.3 + 0.4 * q, q)
    for a in range(5):
        qc.cz(a, a + 1)
    for q in range(6):
        qc.rx(0.2 + 0.1 * q, q)
    if extra_identity:
        qc.x(0)
        qc.x(0)
    return qc


_CUT = {"gate": ["CZ_9"]}
_HAMILTONIAN = [(0.5, "ZZZZZZ"), (0.3, "ZIZIZI"), (-0.2, "XXIIZZ")]


class _Runs:
    """Counts the fragment circuits the engine runs (one task each)."""

    def __enter__(self):
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc

        self.mc, self.original, self.count = mc, mc._run_fragment, 0

        def counting(*args, **kwargs):
            self.count += 1
            return self.original(*args, **kwargs)

        mc._run_fragment = counting
        return self

    def __exit__(self, *exc):
        self.mc._run_fragment = self.original


@unittest.skipUnless(QISKIT, "qiskit required")
class EngineCacheTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd

        self.qd = qd
        self.store = qd.InMemoryCircuitCache()

    def _cut(self, observable=_HAMILTONIAN, **kwargs):
        with _Runs() as runs:
            value = self.qd.circuit_cutting(_circuit(), _CUT, observable, **kwargs)
        return value, runs.count

    def test_exact_results_are_reused(self):
        reference, plain = self._cut(shots=None)
        first, runs = self._cut(shots=None, cache=self.store)
        self.assertLess(runs, plain)  # the terms share the fragments they don't touch
        second, again = self._cut(shots=None, cache=self.store)
        self.assertEqual(again, 0)
        self.assertEqual(first, reference)
        self.assertEqual(second, reference)
        self.assertGreater(self.store.stats()["hits"], 0)

    def test_a_shot_result_serves_later_runs_and_other_observables(self):
        first, _ = self._cut(shots=4000, seed=1, cache=self.store)
        second, runs = self._cut(shots=4000, seed=2, cache=self.store)
        self.assertEqual(runs, 0)
        self.assertEqual(first, second)  # the same samples, reused
        other, runs = self._cut(
            observable="ZIIIIZ", shots=4000, seed=5, cache=self.store
        )
        self.assertEqual(runs, 0)  # the stored counts give any correlator
        self.assertLess(abs(other - _exact(_circuit(), "ZIIIIZ")), 0.08)

    def test_refresh_runs_again_and_replaces(self):
        first, _ = self._cut(shots=4000, seed=1, cache=self.store)
        fresh, runs = self._cut(
            shots=4000, seed=9, cache=self.store, cache_refresh=True
        )
        self.assertGreater(runs, 0)
        self.assertNotEqual(first, fresh)
        after, runs = self._cut(shots=4000, seed=1, cache=self.store)
        self.assertEqual((after, runs), (fresh, 0))

    def test_shots_and_backend_are_part_of_the_key(self):
        self._cut(shots=4000, seed=1, cache=self.store)
        _, runs = self._cut(shots=2000, seed=1, cache=self.store)
        self.assertGreater(runs, 0)
        _, runs = self._cut(shots=None, cache=self.store)
        self.assertGreater(runs, 0)
        _, runs = self._cut(shots=4000, seed=1, backend="qibo", cache=self.store)
        self.assertGreater(runs, 0)

    def test_every_cut_function_and_distributions(self):
        qd, qc = self.qd, _circuit()
        for call in (
            lambda **kw: qd.gate_cutting(qc, ["CZ_9"], "ZZZZZZ", **kw),
            lambda **kw: qd.gate_cutting_hadamard(qc, ["CZ_9"], "ZZZZZZ", **kw),
            lambda **kw: qd.wire_cutting(qc, [("CZ_8", "CZ_9")], "ZZZZZZ", **kw),
            lambda **kw: qd.circuit_cutting_sampling(qc, _CUT, **kw),
        ):
            first = call(shots=None, cache=self.store)
            with _Runs() as runs:
                second = call(shots=None, cache=self.store)
            self.assertEqual(runs.count, 0)
            self.assertEqual(first, second)

    def test_use_cache_for_a_session(self):
        qd = self.qd
        with qd.use_cache(self.store):
            self._cut(shots=None)
            _, runs = self._cut(shots=None)
            self.assertEqual(runs, 0)
            _, runs = self._cut(shots=None, cache=False)  # off for one call
            self.assertGreater(runs, 0)
        _, runs = self._cut(shots=None)  # off again after the block
        self.assertGreater(runs, 0)
        qd.use_cache(self.store)
        try:
            _, runs = self._cut(shots=None)
            self.assertEqual(runs, 0)
        finally:
            qd.use_cache(None)

    @unittest.skipUnless(PYZX, "pyzx required")
    def test_semantic_key_matches_equivalent_circuits(self):
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc

        qd = self.qd
        simulations = []
        original = mc._exact_fragment

        def counting(*args, **kwargs):
            simulations.append(1)
            return original(*args, **kwargs)

        mc._exact_fragment = counting
        try:
            for cache_key in ("exact", "semantic"):
                store = qd.InMemoryCircuitCache()
                qd.circuit_cutting(
                    _circuit(),
                    _CUT,
                    "ZZZZZZ",
                    shots=None,
                    cache=store,
                    cache_key=cache_key,
                )
                simulations.clear()
                value = qd.circuit_cutting(
                    _circuit(extra_identity=True),
                    _CUT,
                    "ZZZZZZ",
                    shots=None,
                    cache=store,
                    cache_key=cache_key,
                )
                self.assertAlmostEqual(value, _exact(_circuit(), "ZZZZZZ"), places=9)
                if cache_key == "exact":
                    missed = len(simulations)  # a different circuit: all run
                    self.assertGreater(missed, 0)
                else:
                    # equivalent: served, but ZX reduction is not a canonical form,
                    # so a circuit can still be missed (never served wrongly)
                    self.assertLess(len(simulations), missed)
        finally:
            mc._exact_fragment = original

    @unittest.skipUnless(LMDB, "lmdb required")
    def test_a_persistent_store(self):
        qd = self.qd
        with tempfile.TemporaryDirectory() as path:
            first, _ = self._cut(shots=4000, seed=1, cache=f"lmdb:{path}")
            second, runs = self._cut(shots=4000, seed=2, cache=f"lmdb:{path}")
            self.assertEqual((second, runs), (first, 0))
            exact, _ = self._cut(shots=None, cache=f"lmdb:{path}")
            again, runs = self._cut(shots=None, cache=f"lmdb:{path}")
            self.assertEqual((again, runs), (exact, 0))
        del qd

    def test_a_backends_options_are_part_of_the_key(self):
        # a noisy backend's results must not be served for a noiseless one
        from qiskit_aer.noise import NoiseModel, depolarizing_error

        noise = NoiseModel()
        noise.add_all_qubit_quantum_error(
            depolarizing_error(0.2, 1), ["ry", "rx", "rz", "z", "s", "sdg", "x"]
        )
        qd = self.qd
        noisy = qd.build_backend("AER_BACKEND", qd.BackendConfig(noise_model=noise))
        plain = qd.build_backend("AER_BACKEND", qd.BackendConfig())
        first, runs = self._cut(
            observable="ZZZZZZ", backend=noisy, shots=2000, seed=1, cache=self.store
        )
        self.assertGreater(runs, 0)
        second, runs = self._cut(
            observable="ZZZZZZ", backend=plain, shots=2000, seed=1, cache=self.store
        )
        self.assertGreater(runs, 0)
        self.assertNotEqual(first, second)
        # the same options, in another object: the same results
        again = qd.build_backend("AER_BACKEND", qd.BackendConfig())
        third, runs = self._cut(
            observable="ZZZZZZ", backend=again, shots=2000, seed=5, cache=self.store
        )
        self.assertEqual((third, runs), (second, 0))

    def test_how_a_backend_connects_is_not_part_of_the_key(self):
        from Qdislib.core.backend.base_backend import BackendConfig, BaseBackend
        from Qdislib.core.cutting_algorithms.engine_cache import backend_id

        def device(**extra):
            return BaseBackend("DEV", BackendConfig(extra=extra))

        self.assertEqual(
            backend_id(device(token="a", backend_name="x"), 10),
            backend_id(device(token="b", backend_name="x"), 10),
        )
        self.assertNotEqual(
            backend_id(device(backend_name="x"), 10),
            backend_id(device(backend_name="y"), 10),
        )
        self.assertNotEqual(
            backend_id(device(options={"a": 1}), 10),
            backend_id(device(options={"a": 2}), 10),
        )
        # the device a backend chose by itself counts as one it was told
        chosen = device()
        chosen._ibm_backend = mock.Mock()
        chosen._ibm_backend.name = "x"
        self.assertEqual(
            backend_id(chosen, 10), backend_id(device(backend_name="x"), 10)
        )
        self.assertTrue(backend_id(device(), None) == "exact")

    def test_expectation_value_with_cache_true_reuses(self):
        # it built a new store on every call: nothing was ever found
        subs = self.qd.circuit_cutting_subcircuits(_circuit(), _CUT, "ZZZZZZ")
        with mock.patch.dict(engine_cache._STORES, clear=True):
            with _Runs() as first:
                values = [
                    self.qd.expectation_value(s, shots=200, cache=True, seed=1)
                    for s in subs
                ]
            with _Runs() as second:
                again = [
                    self.qd.expectation_value(s, shots=200, cache=True, seed=2)
                    for s in subs
                ]
        self.assertGreater(first.count, 0)
        self.assertEqual(second.count, 0)
        self.assertEqual(
            subs.reconstruct(values), subs.reconstruct(again)
        )  # the same samples

    def test_expectation_value_with_cache_false_is_off(self):
        subs = self.qd.circuit_cutting_subcircuits(_circuit(), _CUT, "ZZZZZZ")
        with _Runs() as runs:
            self.qd.expectation_value(subs[0], shots=100, cache=False, seed=1)
            self.qd.expectation_value(subs[0], shots=100, cache=False, seed=1)
        self.assertEqual(runs.count, 2)

    def test_a_key_tells_complex_numbers_apart(self):
        import numpy as np

        for kind in (np.complex64, np.complex128, complex):
            self.assertNotEqual(
                engine_cache._canonical(kind(1, 1)),
                engine_cache._canonical(kind(1, -1)),
            )
        self.assertEqual(engine_cache._canonical(np.float64(0.5)), repr(0.5))

    def test_a_redis_url_keeps_its_credentials(self):
        with mock.patch.dict(engine_cache._STORES, clear=True):
            with mock.patch(
                "Qdislib.utils.circuit_cache.default_circuit_cache"
            ) as build:
                engine_cache._store("redis://me@host:7000/3")
        build.assert_called_once_with(
            "redis", host="host", port=7000, db=3, username="me"
        )

    def test_a_tls_redis_url(self):
        with mock.patch.dict(engine_cache._STORES, clear=True):
            with mock.patch(
                "Qdislib.utils.circuit_cache.default_circuit_cache"
            ) as build:
                engine_cache._store("rediss://me@host:7000/3")
        build.assert_called_once_with(
            "redis", host="host", port=7000, db=3, username="me", ssl=True
        )

    def test_a_password_in_a_url_is_refused_and_not_printed(self):
        # a store goes to PyCOMPSs tasks, whose arguments are written to files
        with mock.patch.dict(engine_cache._STORES, clear=True):
            with self.assertRaisesRegex(ValueError, "QDISLIB_REDIS_PASSWORD") as caught:
                engine_cache._store("redis://me:s3cret@host:7000/3")
        self.assertNotIn("s3cret", str(caught.exception))

    def test_a_refused_url_does_not_print_its_password(self):
        for spec in (
            "rediz://me:s3cret@host:7000/3",
            "rediz://:s3cret@host",
            "rediz://me:pa/s3cret@host",
            "rediz://me:p@s3cret@host",
        ):
            with mock.patch.dict(engine_cache._STORES, clear=True):
                with self.assertRaises(ValueError) as caught:
                    engine_cache._store(spec)
            self.assertNotIn("s3cret", str(caught.exception))
            self.assertIn("rediz://", str(caught.exception))

    def test_cache_mode_is_deprecated(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self.qd.gate_cutting(
                _circuit(),
                ["CZ_9"],
                "ZZZZZZ",
                shots=None,
                cache=self.store,
                cache_mode="counts",
            )
        self.assertTrue(
            any("cache_mode is deprecated" in str(w.message) for w in caught)
        )


@unittest.skipUnless(QISKIT, "qiskit required")
class EstimateMetricsTest(BaseTimedTestCase):
    """estimate's metrics report the cache the call used, however it was chosen."""

    def _estimate(self, **options):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.3, 1)
        qc.cz(1, 2)
        qc.ry(0.5, 2)
        _, metrics = qd.estimate(
            qc, "ZZZ", max_qubits=2, shots=None, return_metrics=True, **options
        )
        return metrics["cache"]

    def test_the_session_cache(self):
        import Qdislib as qd

        with qd.use_cache("memory") as store:
            stats = self._estimate()
            self.assertIsNotNone(stats)
            self.assertEqual(stats, store.stats())
            self.assertGreater(stats["size"], 0)

    def test_its_own_cache(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        store = InMemoryCircuitCache()
        self.assertEqual(self._estimate(cache=store), store.stats())

    def test_off(self):
        import Qdislib as qd

        self.assertIsNone(self._estimate())
        with qd.use_cache("memory"):
            self.assertIsNone(self._estimate(cache=False))


if __name__ == "__main__":
    unittest.main()
