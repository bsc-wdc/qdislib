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

"""Tests for the semantic quantum-circuit cache (Qdislib.utils.circuit_cache)."""

import importlib.util
import os
import pickle
import tempfile
import unittest
import unittest.mock

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None
PYZX_AVAILABLE = importlib.util.find_spec("pyzx") is not None
LMDB_AVAILABLE = importlib.util.find_spec("lmdb") is not None
CUDAQ_AVAILABLE = importlib.util.find_spec("cudaq") is not None


def _qibo_circuit():
    """A small fixed two-qubit Qibo circuit (fresh object each call).

    Qibo evaluators measure in place, so every evaluation needs its own circuit.
    """
    from qibo import models, gates

    circ = models.Circuit(2)
    circ.add(gates.H(0))
    circ.add(gates.CZ(0, 1))
    circ.add(gates.RX(1, 0.3))
    return circ


def _nonzero_circuit():
    """The suite's reference circuit, whose <ZZ> ~ -0.717 (clearly non-zero).

    With the RY *after* the CZ the analytical value is far from zero, so a cache
    error could not hide near zero -- a strong correctness check.
    """
    from qibo import models, gates

    circ = models.Circuit(2)
    circ.add(gates.H(0))
    circ.add(gates.CZ(0, 1))
    circ.add(gates.RY(0, 0.8))
    return circ


class CacheStoreTest(BaseTimedTestCase):
    """The store backends share one get/put/contains interface."""

    def test_in_memory_roundtrip(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        self.assertIsNone(cache.get("missing"))
        self.assertNotIn("missing", cache)
        cache.put("k", 0.5)
        self.assertEqual(cache.get("k"), 0.5)
        self.assertIn("k", cache)
        self.assertEqual(len(cache), 1)

    def test_in_memory_lru_eviction(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache(max_entries=2)
        cache.put("a", 1)
        cache.put("b", 2)
        cache.get("a")  # touch 'a' so 'b' becomes least-recently-used
        cache.put("c", 3)  # over the cap -> evict the LRU entry ('b')
        self.assertEqual(len(cache), 2)
        self.assertIsNone(cache.get("b"))
        self.assertEqual(cache.get("a"), 1)
        self.assertEqual(cache.get("c"), 3)

    def test_in_memory_ttl_expiry(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        expired = InMemoryCircuitCache(ttl=0)  # entries expire immediately
        expired.put("k", 1)
        self.assertIsNone(expired.get("k"))
        self.assertNotIn("k", expired)

        fresh = InMemoryCircuitCache(ttl=1000)  # long ttl -> still cached
        fresh.put("k", 1)
        self.assertEqual(fresh.get("k"), 1)

    def test_default_cache_backend_selection(self):
        from Qdislib.utils.circuit_cache import (
            default_circuit_cache,
            InMemoryCircuitCache,
            LMDBCircuitCache,
        )

        # The caller chooses the backend explicitly; "memory" is the default.
        self.assertIsInstance(default_circuit_cache(), InMemoryCircuitCache)
        self.assertIsInstance(default_circuit_cache("memory"), InMemoryCircuitCache)
        with self.assertRaises(ValueError):
            default_circuit_cache("lmdb")  # path is required
        with self.assertRaises(ValueError):
            default_circuit_cache("nonsense")
        if LMDB_AVAILABLE:
            with tempfile.TemporaryDirectory() as tmp:
                cache = default_circuit_cache("lmdb", path=os.path.join(tmp, "c"))
                self.assertIsInstance(cache, LMDBCircuitCache)
                cache.close()

    @unittest.skipUnless(LMDB_AVAILABLE, "lmdb not installed")
    def test_lmdb_roundtrip_and_persistence(self):
        from Qdislib.utils.circuit_cache import LMDBCircuitCache

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache.lmdb")
            cache = LMDBCircuitCache(path)
            self.assertIsNone(cache.get("k"))
            cache.put("k", {"01": 7, "10": 3})
            self.assertEqual(cache.get("k"), {"01": 7, "10": 3})
            self.assertEqual(len(cache), 1)
            cache.close()  # release this process's handle so the reopen hits disk

            # A fresh handle on the same path sees the persisted entry.
            self.assertEqual(LMDBCircuitCache(path).get("k"), {"01": 7, "10": 3})
            LMDBCircuitCache(path).close()

    @unittest.skipUnless(LMDB_AVAILABLE, "lmdb not installed")
    def test_lmdb_is_picklable(self):
        # The store must survive pickling (e.g. to a PyCOMPSs task): only the
        # config travels, the live env handle does not.
        from Qdislib.utils.circuit_cache import LMDBCircuitCache

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache.lmdb")
            cache = LMDBCircuitCache(path)
            cache.put("k", 1.0)
            data = pickle.dumps(cache)
            cache.close()  # release the handle before reviving in this process
            revived = pickle.loads(data)
            self.assertEqual(revived.get("k"), 1.0)
            revived.close()


@unittest.skipUnless(
    QIBO_AVAILABLE and PYZX_AVAILABLE, "qibo and pyzx required for the semantic key"
)
class SemanticKeyTest(BaseTimedTestCase):
    """semantic_key is stable, equivalence-aware and fails soft."""

    def test_identical_circuits_share_a_key(self):
        from Qdislib.utils.circuit_cache import semantic_key

        key = semantic_key(_qibo_circuit())
        self.assertIsInstance(key, str)
        self.assertEqual(key, semantic_key(_qibo_circuit()))

    def test_different_circuits_differ(self):
        from qibo import models, gates
        from Qdislib.utils.circuit_cache import semantic_key

        other = models.Circuit(2)
        other.add(gates.H(0))
        other.add(gates.H(1))
        self.assertNotEqual(semantic_key(_qibo_circuit()), semantic_key(other))

    def test_semantically_equivalent_circuits_share_a_key(self):
        # An appended self-inverse pair (X then X = identity) is removed by the
        # ZX full_reduce, so the key must match the circuit without it.
        from qibo import gates
        from Qdislib.utils.circuit_cache import semantic_key

        base = _qibo_circuit()
        padded = _qibo_circuit()
        padded.add(gates.X(0))
        padded.add(gates.X(0))
        self.assertEqual(semantic_key(base), semantic_key(padded))

    def test_qubit_permuted_circuits_differ(self):
        """bsc-wdc/qdislib#11: X on qubit 0 and X on qubit 1 must not share a key."""
        from qibo import models, gates
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache, semantic_key

        def x_on(qubit):
            circuit = models.Circuit(2)
            circuit.add(gates.X(qubit))
            return circuit

        self.assertNotEqual(semantic_key(x_on(0)), semantic_key(x_on(1)))

        cache = InMemoryCircuitCache()
        self.assertEqual(
            expectation_value([x_on(0), "ZI"], shots=100, cache=cache), -1.0
        )
        self.assertEqual(
            expectation_value([x_on(1), "ZI"], shots=100, cache=cache), 1.0
        )

    def test_unsupported_backend_returns_none(self):
        from Qdislib.utils.circuit_cache import semantic_key

        # An unrecognised backend (e.g. CUDA-Q) cannot be keyed -> None, not error.
        self.assertIsNone(semantic_key(_qibo_circuit(), software="cudaq"))

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_qiskit_key(self):
        from qiskit import QuantumCircuit
        from Qdislib.utils.circuit_cache import semantic_key

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        self.assertIsInstance(semantic_key(qc), str)


@unittest.skipUnless(
    QIBO_AVAILABLE and PYZX_AVAILABLE, "qibo and pyzx required for cached evaluation"
)
class CachedExpectationTest(BaseTimedTestCase):
    """expectation_value(..., cache=...) reuses results and stays correct."""

    def test_expectation_mode_reuses_and_returns_stored_value(self):
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        first = expectation_value([_qibo_circuit(), "ZZ"], shots=2048, cache=cache)
        second = expectation_value([_qibo_circuit(), "ZZ"], shots=2048, cache=cache)
        # The second identical circuit is a hit: one stored entry, exact reuse
        # (so no sampling-noise difference between the two values).
        self.assertEqual(len(cache), 1)
        self.assertEqual(first, second)

    def test_counts_mode_reuses_counts_across_observables(self):
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        # Same topology, two observables: counts mode keys on topology alone, so
        # a single counts entry is reused to derive both expectations.
        expectation_value(
            [_qibo_circuit(), "ZZ"], shots=2048, cache=cache, cache_mode="counts"
        )
        expectation_value(
            [_qibo_circuit(), "ZI"], shots=2048, cache=cache, cache_mode="counts"
        )
        self.assertEqual(len(cache), 1)

    def test_cache_true_uses_default_store(self):
        from Qdislib.utils.backends import expectation_value

        value = expectation_value([_qibo_circuit(), "ZZ"], shots=1024, cache=True)
        self.assertIsInstance(value, float)

    def test_custom_strategy_stores_user_payload(self):
        # A user-defined strategy can store its own value shape (here counts plus
        # metadata) and read it back -- the extensibility hook.
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache, CountsStrategy

        class CountsWithMeta(CountsStrategy):
            def store(self, sim):
                counts = sim.counts()
                payload = {
                    "counts": counts,
                    "shots": sim.shots,
                    "observable": sim.observables,
                }
                return payload, sim.expectation_from_counts(counts)

            def load(self, cached, sim):
                return sim.expectation_from_counts(cached["counts"])

        cache = InMemoryCircuitCache()
        first = expectation_value(
            [_qibo_circuit(), "ZZ"],
            shots=2048,
            cache=cache,
            cache_mode=CountsWithMeta(),
        )
        second = expectation_value(
            [_qibo_circuit(), "ZZ"],
            shots=2048,
            cache=cache,
            cache_mode=CountsWithMeta(),
        )
        self.assertEqual(len(cache), 1)
        self.assertEqual(first, second)
        stored = next(iter(cache._store.values()))
        self.assertIn("counts", stored)
        self.assertEqual(stored["shots"], 2048)
        self.assertEqual(stored["observable"], "ZZ")

    def test_caching_preserves_wire_cutting_result(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_nonzero_circuit(), "ZZ")
        # Guard: the reference value is genuinely non-zero, so a wrong cache
        # result could not slip through near zero.
        self.assertGreater(abs(analytical), 0.5)

        subs = qd.wire_cutting_subcircuits(_nonzero_circuit(), [("H_1", "CZ_2")])
        cache = qd.InMemoryCircuitCache()
        values = [qd.expectation_value(s, shots=20000, cache=cache) for s in subs]
        self.assertLess(abs(subs.reconstruct(values) - analytical), 0.2)

    def test_caching_preserves_gate_cutting_result(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_nonzero_circuit(), "ZZ")
        # Guard: a clearly non-zero analytical value (~ -0.717).
        self.assertGreater(abs(analytical), 0.5)

        subs = qd.gate_cutting_subcircuits(_nonzero_circuit(), ["CZ_2"])
        cache = qd.InMemoryCircuitCache()
        values = [qd.expectation_value(s, shots=20000, cache=cache) for s in subs]
        # Gate cutting with the cache must still match the (non-zero) analytical.
        self.assertLess(abs(subs.reconstruct(values) - analytical), 0.2)

    def test_all_in_one_gate_cutting_accepts_cache(self):
        # The one-call API also caches (delegates to the subcircuit workflow).
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_nonzero_circuit(), "ZZ")
        cache = qd.InMemoryCircuitCache()
        result = qd.gate_cutting(
            _nonzero_circuit(), ["CZ_2"], shots=20000, backend="qibo", cache=cache
        )
        self.assertLess(abs(result - analytical), 0.2)
        self.assertGreater(len(cache), 0)  # the all-in-one path populated the cache

    def test_all_in_one_wire_cutting_accepts_cache(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        analytical = qd.analytical_solution(_nonzero_circuit(), "ZZ")
        cache = qd.InMemoryCircuitCache()
        result = qd.wire_cutting(
            _nonzero_circuit(),
            [("H_1", "CZ_2")],
            shots=20000,
            backend="qibo",
            cache=cache,
        )
        self.assertLess(abs(result - analytical), 0.2)
        self.assertGreater(len(cache), 0)


class _FakeRedisClient:
    """What ``redis.Redis`` is used for, in a shared dict."""

    data = {}

    def __init__(self, **settings):
        self.settings = settings

    def get(self, key):
        return self.data.get(key)

    def set(self, key, value):
        self.data[key] = value

    def delete(self, key):
        self.data.pop(key, None)


class RedisStoreTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        from unittest import mock

        from Qdislib.utils import circuit_cache as cc

        _FakeRedisClient.data = {}
        fake = mock.Mock(Redis=_FakeRedisClient)
        for name, value in (("_redis", fake), ("redis_available", True)):
            patcher = mock.patch.object(cc, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.cc = cc

    def test_a_value_goes_in_and_out_as_json_under_its_prefix(self):
        store = self.cc.RedisCircuitCache(host="h", port=1, db=2, password="p")
        store.put("k", {"a": [1, 2]})
        self.assertEqual(list(_FakeRedisClient.data), ["qdislib:k"])
        self.assertEqual(store.get("k"), {"a": [1, 2]})
        self.assertIsNone(store.get("other"))
        self.assertEqual(store.stats()["hits"], 1)
        self.assertEqual(store.stats()["misses"], 1)
        self.assertIsNone(store.stats()["size"])  # a shared server: not counted
        self.assertEqual(store._client.settings["password"], "p")
        store.delete("k")
        self.assertIsNone(store.get("k"))

    def test_goes_to_a_task_by_its_settings_not_its_connection(self):
        store = self.cc.RedisCircuitCache(host="h", prefix="x:")
        store.put("k", 1.5)
        store._client.unpicklable = lambda: None
        again = pickle.loads(pickle.dumps(store))
        self.assertIsNone(again._client)
        self.assertEqual((again.host, again.prefix), ("h", "x:"))
        self.assertEqual(again.get("k"), 1.5)

    def test_redis_outside_compss_is_this_store(self):
        with unittest.mock.patch.object(self.cc, "pycompss_available", False):
            store = self.cc.default_circuit_cache("redis", host="h")
        self.assertIsInstance(store, self.cc.RedisCircuitCache)

    def test_without_the_package_it_says_so(self):
        with unittest.mock.patch.object(self.cc, "redis_available", False):
            with self.assertRaisesRegex(ImportError, "redis"):
                self.cc.RedisCircuitCache()


@unittest.skipUnless(PYZX_AVAILABLE, "pyzx required")
class PennylaneKeyTest(BaseTimedTestCase):
    def test_a_gate_without_a_mnemonic_gets_no_key(self):
        # left out of the QASM, it would not change the key: U3(0.1) and U3(1.1)
        # shared one, and a cached result
        from Qdislib.utils.circuit_cache import semantic_key
        from Qdislib.utils.graph_pennylane import PennylaneCircuit

        circuit = PennylaneCircuit([("H", (0,), []), ("U3", (0,), [0.1, 0.2, 0.3])], 1)
        self.assertIsNone(semantic_key(circuit, "pennylane"))

    def test_gates_with_a_mnemonic_and_the_markers(self):
        from Qdislib.utils.circuit_cache import semantic_key
        from Qdislib.utils.graph_pennylane import PennylaneCircuit

        def key(*ops):
            return semantic_key(PennylaneCircuit(list(ops), 2), "pennylane")

        bell = [("H", (0,), []), ("CX", (0, 1), [])]
        self.assertIsNotNone(key(*bell))
        self.assertEqual(key(*bell), key(*bell, ("MEASURE", (0,), [])))
        self.assertNotEqual(key(*bell), key(*bell, ("RZ", (1,), [0.3])))


class LoggingControlTest(BaseTimedTestCase):
    """set_log_level replaces the old per-function verbose= flags."""

    def test_set_log_level(self):
        import logging
        from Qdislib.utils.logging_config import set_log_level

        logger = set_log_level("DEBUG")
        self.assertEqual(logger.level, logging.DEBUG)
        self.assertTrue(
            any(getattr(h, "_qdislib_handler", False) for h in logger.handlers)
        )
        # Idempotent: calling again must not add a second handler.
        handler_count = len(logger.handlers)
        set_log_level("WARNING")
        self.assertEqual(len(logger.handlers), handler_count)
        self.assertEqual(logger.level, logging.WARNING)


class PlaceholderRefreshTest(BaseTimedTestCase):
    """The HPC concurrency trick: reserve a placeholder, defer, refresh later."""

    def test_is_pending(self):
        from Qdislib.utils.circuit_cache import _PENDING, is_pending

        self.assertTrue(is_pending(_PENDING))
        self.assertFalse(is_pending(0.5))
        self.assertFalse(is_pending({"01": 7}))
        self.assertFalse(is_pending(None))

    def test_pending_result_resolves_to_stored_value(self):
        from Qdislib.utils.circuit_cache import (
            InMemoryCircuitCache,
            ExpectationStrategy,
            PendingResult,
        )

        cache = InMemoryCircuitCache()
        cache.put("k", 0.42)
        pending = PendingResult(cache, "k", ExpectationStrategy(), "ZZ", "qibo", 1024)
        self.assertEqual(pending.resolve(), 0.42)

    def test_pending_result_gives_up_when_unresolved(self):
        from Qdislib.utils.circuit_cache import (
            InMemoryCircuitCache,
            ExpectationStrategy,
            PendingResult,
            _PENDING,
        )
        from Qdislib.utils.exceptions import QdislibException

        cache = InMemoryCircuitCache()
        cache.put("k", _PENDING)  # never replaced by a real value
        pending = PendingResult(
            cache, "k", ExpectationStrategy(), "ZZ", "qibo", 1024, retries=1, delay=0
        )
        # a value left out would silently change a product of fragment values
        with self.assertRaisesRegex(QdislibException, "never computed"):
            pending.resolve()
        self.assertIsNone(cache.get("k"))  # the next run computes it

    @unittest.skipUnless(QIBO_AVAILABLE and PYZX_AVAILABLE, "qibo and pyzx required")
    def test_concurrent_placeholder_defers_then_refreshes(self):
        # Simulate a concurrent task: plant the placeholder another task would
        # have written, so a lookup must defer (PendingResult) instead of
        # re-simulating; once the real value lands, the deferred result refreshes.
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import (
            InMemoryCircuitCache,
            PendingResult,
            semantic_key,
            _PENDING,
        )

        cache = InMemoryCircuitCache()
        key = semantic_key(_qibo_circuit(), "qibo")
        self.assertIsNotNone(key)
        key = f"{key}|shots=1024|ZZ"
        cache.put(key, _PENDING)  # another task is computing this key

        deferred = expectation_value([_qibo_circuit(), "ZZ"], cache=cache)
        self.assertIsInstance(deferred, PendingResult)

        cache.put(key, 0.25)  # owning task stored the real value
        self.assertEqual(deferred.resolve(), 0.25)

    @unittest.skipUnless(QIBO_AVAILABLE, "qibo required to import reconstruction")
    def test_reconstruction_resolves_pending(self):
        # The master-side refresh: a PendingResult fed to reconstruction is
        # resolved from the cache (after the COMPSs barrier) before combining.
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            wire_cutting_subcircuit_reconstruction,
        )
        from Qdislib.utils.circuit_cache import (
            InMemoryCircuitCache,
            ExpectationStrategy,
            PendingResult,
        )

        cache = InMemoryCircuitCache()
        cache.put("k|ZZ", 0.5)  # the owning task's stored value
        pending = PendingResult(
            cache, "k|ZZ", ExpectationStrategy(), "ZZ", "qibo", 1024
        )
        result = wire_cutting_subcircuit_reconstruction([pending, 1.0], number_cuts=0)
        self.assertAlmostEqual(result, 0.5)


@unittest.skipUnless(
    QISKIT_AVAILABLE and PYZX_AVAILABLE, "qiskit and pyzx required for the shot keys"
)
class ShotsInKeyTest(BaseTimedTestCase):
    """A shot run's value belongs to its shots: a later run with others runs."""

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.rz(0.3, 0)
        qc.h(0)
        qc.cz(0, 1)
        return qc

    def _assert_not_reused(self, payload, cache, **options):
        from Qdislib.utils.backends import expectation_value

        # <ZZ> = cos(0.3) = 0.955. Four shots can only give multiples of 0.5.
        expectation_value(payload, shots=4, seed=1, cache=cache, **options)
        value = expectation_value(payload, shots=100000, seed=1, cache=cache, **options)
        self.assertAlmostEqual(value, 0.9553, delta=0.02)

    def test_native_circuit(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        self._assert_not_reused(self._circuit(), InMemoryCircuitCache())

    def test_circuit_observable_pair(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        self._assert_not_reused([self._circuit(), "ZZ"], InMemoryCircuitCache())

    def test_the_shared_store_of_cache_true(self):
        # Its own circuit: the store outlives the test.
        qc = self._circuit()
        qc.rz(0.123456, 1)
        self._assert_not_reused(qc, True)

    def test_counts_mode(self):
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        cache = InMemoryCircuitCache()
        self._assert_not_reused(self._circuit(), cache, cache_mode="counts")
        self.assertEqual(len(cache), 2)

    def test_custom_strategy_keeps_its_own_key_hook(self):
        from Qdislib.utils.circuit_cache import (
            ExpectationStrategy,
            InMemoryCircuitCache,
        )

        class Keyed(ExpectationStrategy):
            def make_key(self, semantic_key, sim):
                return f"mine:{semantic_key}"

        cache = InMemoryCircuitCache()
        self._assert_not_reused(self._circuit(), cache, cache_mode=Keyed())
        self.assertEqual(len(cache), 2)

    def test_subcircuit_of_a_list(self):
        from Qdislib.utils.backends import Subcircuit, convert_to_dag
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        subcircuit = Subcircuit(convert_to_dag(self._circuit()), 2, "qiskit")
        self._assert_not_reused(subcircuit, InMemoryCircuitCache())

    def test_the_deprecated_cut_paths(self):
        import warnings

        from qiskit import QuantumCircuit

        import Qdislib as qd
        from Qdislib.utils.circuit_cache import InMemoryCircuitCache

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.rz(0.3, 0)
        qc.h(0)
        qc.cz(0, 1)
        qc.cz(1, 2)
        exact = qd.analytical_solution(qc, "ZZZ")
        for run in (
            lambda **k: qd.gate_cutting(qc, ["CZ_5"], "ZZZ", method="automatic", **k),
            lambda **k: qd.wire_cutting(
                qc, [("CZ_4", "CZ_5")], "ZZZ", method="automatic", **k
            ),
        ):
            cache = InMemoryCircuitCache()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                run(shots=4, seed=1, cache=cache)
                size = len(cache)
                value = run(shots=100000, seed=1, cache=cache)
            self.assertEqual(len(cache), 2 * size)
            self.assertAlmostEqual(value, exact, delta=0.1)


@unittest.skipUnless(
    QIBO_AVAILABLE and CUDAQ_AVAILABLE and PYZX_AVAILABLE,
    "qibo, cudaq and pyzx required for CUDA-Q caching",
)
class CudaqCacheTest(BaseTimedTestCase):
    """CUDA-Q kernels are cacheable via cudaq.translate (when it supports them)."""

    def _kernel(self):
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel

        dag = circuit_qibo_to_dag(_qibo_circuit())
        return dag_to_circuit_cudaq_kernel(dag, 2)

    def test_semantic_key_is_str_or_none(self):
        from Qdislib.utils.circuit_cache import semantic_key

        try:
            kernel, _ = self._kernel()
        except Exception:
            self.skipTest("could not build a CUDA-Q kernel here")
        key = semantic_key(kernel, "cudaq")
        # cudaq.translate supports many kernels but not all; either a key or a
        # graceful None (caching skipped) is acceptable.
        self.assertTrue(key is None or isinstance(key, str))

    def test_cached_evaluation_reuses(self):
        from Qdislib.utils.backends import expectation_value
        from Qdislib.utils.circuit_cache import semantic_key, InMemoryCircuitCache

        try:
            kernel, _ = self._kernel()
        except Exception:
            self.skipTest("could not build a CUDA-Q kernel here")
        if semantic_key(kernel, "cudaq") is None:
            self.skipTest("cudaq.translate did not produce a key for this kernel")

        cache = InMemoryCircuitCache()
        k1, o1 = self._kernel()
        k2, o2 = self._kernel()
        first = expectation_value([k1, o1], shots=512, cache=cache)
        second = expectation_value([k2, o2], shots=512, cache=cache)
        self.assertEqual(len(cache), 1)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
