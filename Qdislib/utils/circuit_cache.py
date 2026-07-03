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

"""Semantic quantum-circuit cache for Qdislib.

Content-addressed ("semantic") cache for subcircuit results. A circuit is
reduced to a canonical ZX-diagram with PyZX (``full_reduce``), converted to a
``networkx`` graph and hashed with a Weisfeiler-Lehman graph hash; that hash is
the cache key. So *functionally equivalent* subcircuits -- common when cutting
repeats fragments, or across iterations of a variational algorithm -- map to the
same key and are only simulated once.

The key (:func:`semantic_key`) is backend-agnostic (Qibo / Qiskit / CUDA-Q, each
via its OpenQASM-2 export). The caller picks a store, all sharing one
``get``/``put`` interface:

- :class:`InMemoryCircuitCache` -- a process-local ``dict`` (serial use only:
  not shared between PyCOMPSs workers).
- :class:`LMDBCircuitCache` -- a persistent on-disk LMDB store, shareable
  between processes on a common filesystem.
- :class:`RedisCircuitCache` -- a generic redis-py store for distributed sharing
  without COMPSs.
- :class:`RedisPSCOCircuitCache` -- the BSC COMPSs PSCO/Redis store, shared
  across PyCOMPSs tasks and nodes.

:func:`default_circuit_cache` builds one from an explicit ``backend`` choice
(``"memory"``/``"lmdb"``/``"redis"``); ``"redis"`` resolves to the PSCO store
under PyCOMPSs and the redis-py store otherwise. What a hit stores and returns
is governed by a :class:`CacheStrategy` (``"expectation"``, ``"counts"``, or a
custom subclass).

Caching never affects correctness: when a key cannot be computed (PyZX absent,
an unsupported gate, or a CUDA-Q kernel that ``translate`` rejects) the circuit
is simply simulated as usual.
"""

import json
import logging
import math
import time
from fractions import Fraction

import networkx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash

from Qdislib.utils.runtime import pycompss_available

logger = logging.getLogger(__name__)

# --- Optional dependencies (the cache is an optional extra). ---
try:
    import pyzx as zx

    pyzx_available = True
except ImportError:
    zx = None
    pyzx_available = False

try:
    import lmdb

    lmdb_available = True
except ImportError:
    lmdb = None
    lmdb_available = False

try:
    import redis as _redis

    redis_available = True
except ImportError:
    _redis = None
    redis_available = False

# BSC COMPSs persistent storage (Redis-backed PSCOs). Only importable inside a
# COMPSs deployment with a configured storage backend.
try:
    from storage.api import StorageObject, get_by_ID

    redis_psco_available = True
except ImportError:
    StorageObject = object
    get_by_ID = None
    redis_psco_available = False

# Circuit class sentinels for backend detection (each SDK is optional).
try:
    from qiskit.circuit import QuantumCircuit as _QiskitCircuit
except ImportError:
    _QiskitCircuit = type("QuantumCircuit", (), {})
try:
    from qibo.models import Circuit as _QiboCircuit
except ImportError:
    _QiboCircuit = type("QiboCircuit", (), {})
try:
    from cudaq.kernel.kernel_decorator import PyKernelDecorator as _CudaqKernel
except ImportError:
    _CudaqKernel = type("PyKernelDecorator", (), {})
# PennylaneCircuit is a plain Qdislib data wrapper -- always importable, no
# PennyLane install required to reference the class for type detection.
from Qdislib.utils.graph_pennylane import PennylaneCircuit as _PennylaneCircuit


__all__ = [
    "semantic_key",
    "CircuitCache",
    "InMemoryCircuitCache",
    "LMDBCircuitCache",
    "RedisCircuitCache",
    "RedisPSCOCircuitCache",
    "default_circuit_cache",
    "CacheStrategy",
    "ExpectationStrategy",
    "CountsStrategy",
    "PendingResult",
    "is_pending",
]


# --------------------------------------------------------------------------- #
# Semantic key: ZX full-reduce + Weisfeiler-Lehman graph hash
# --------------------------------------------------------------------------- #


def _detect_software(circuit):
    """Return ``"qibo"``/``"qiskit"``/``"cudaq"``/``"pennylane"`` for a circuit, or ``None``."""
    if isinstance(circuit, _QiboCircuit):
        return "qibo"
    if isinstance(circuit, _QiskitCircuit):
        return "qiskit"
    if isinstance(circuit, _CudaqKernel):
        return "cudaq"
    if isinstance(circuit, _PennylaneCircuit):
        return "pennylane"
    return None


# One-/two-qubit OpenQASM 2.0 mnemonics for Qdislib's canonical gate names, used to
# serialise a PennylaneCircuit for the semantic cache.
_QASM_ONE_QUBIT = {"H": "h", "X": "x", "Y": "y", "Z": "z", "S": "s", "T": "t",
                   "SDG": "sdg", "TDG": "tdg", "SX": "sx"}
_QASM_ROTATION = {"RX": "rx", "RY": "ry", "RZ": "rz"}
_QASM_TWO_QUBIT = {"CX": "cx", "CNOT": "cx", "CZ": "cz", "CY": "cy"}


def _pennylane_circuit_to_qasm(circuit):
    """Serialise a :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit` to OpenQASM 2.0.

    Only the unitary gates are emitted -- the gate-cut ``MEASURE`` (collapse) and
    ``Observable I`` markers are skipped, matching how the semantic hash treats the
    other backends (the collapse is folded into the key separately, via
    ``_cut_point_fingerprint``). A gate PyZX cannot parse simply makes
    :func:`semantic_key` fall back to no caching.
    """
    lines = [
        "OPENQASM 2.0;",
        'include "qelib1.inc";',
        f"qreg q[{circuit.num_qubits}];",
    ]
    for gate, qubits, params in circuit.ops:
        name = gate.upper()
        if name in _QASM_ONE_QUBIT:
            lines.append(f"{_QASM_ONE_QUBIT[name]} q[{qubits[0]}];")
        elif name in _QASM_ROTATION:
            angle = params[0] if params else 0.0
            lines.append(f"{_QASM_ROTATION[name]}({angle}) q[{qubits[0]}];")
        elif name in _QASM_TWO_QUBIT:
            lines.append(f"{_QASM_TWO_QUBIT[name]} q[{qubits[0]}],q[{qubits[1]}];")
        # MEASURE / Observable I: non-unitary markers, intentionally omitted.
    return "\n".join(lines)


def _circuit_to_qasm(circuit, software):
    """Serialise a circuit to OpenQASM 2.0 for PyZX to parse.

    All three SDKs expose an OpenQASM-2 export: Qibo ``Circuit.to_qasm()``,
    Qiskit ``qasm2.dumps()`` and CUDA-Q ``cudaq.translate(kernel,
    format="openqasm2")``. CUDA-Q's exporter only supports argument-free kernels
    and a subset of modifiers; an unsupported kernel raises here and
    :func:`semantic_key` then falls back to no caching.
    """
    if software == "qibo":
        return circuit.to_qasm()
    if software == "qiskit":
        import qiskit.qasm2

        return qiskit.qasm2.dumps(circuit)
    if software == "cudaq":
        import cudaq

        return cudaq.translate(circuit, format="openqasm2")
    if software == "pennylane":
        return _pennylane_circuit_to_qasm(circuit)
    raise ValueError(f"Cannot serialise a {software!r} circuit to QASM for caching.")


def _normalize_phase(phase_str):
    """Normalise a PyZX phase label to a stable string (no ``eval``).

    Phases come out of PyZX as fractions of pi (e.g. ``"1/2"``, ``"3/4"``,
    sometimes carrying a ``π``/``~``). Equivalent circuits must yield identical
    labels, so reduce each to a fixed-precision float; anything unparsable falls
    back to its ASCII characters (still stable).
    """
    s = str(phase_str).strip()
    if s in ("", "0"):
        return "0"
    has_pi = "π" in s or "pi" in s
    cleaned = s.replace("π", "").replace("pi", "").replace("~", "").replace("*", "")
    cleaned = cleaned.strip()
    if cleaned in ("", "+"):
        coeff = Fraction(1)
    elif cleaned == "-":
        coeff = Fraction(-1)
    else:
        try:
            coeff = Fraction(cleaned)
        except (ValueError, ZeroDivisionError):
            return "".join(c for c in s if ord(c) < 128)
    value = float(coeff) * (math.pi if has_pi else 1.0)
    return f"{value:.12f}"


def _zx_dict_to_nx(zx_dict):
    """Convert PyZX's ``graph.to_dict()`` into a labelled ``networkx`` graph."""
    graph = networkx.Graph()
    for vertex in zx_dict["vertices"]:
        graph.add_node(
            vertex["id"],
            type=str(vertex["t"]),
            phase=str(vertex.get("phase", "0")),
        )
    for source, target, edge_type in zx_dict["edges"]:
        graph.add_edge(source, target, type=str(edge_type))
    return graph


def _wl_hash_zx_graph(graph, iterations=3):
    """Weisfeiler-Lehman hash of a ZX graph, keyed on spider type + phase."""
    labelled = graph.copy()
    for node, data in labelled.nodes(data=True):
        spider = data.get("type", "")
        phase = _normalize_phase(data.get("phase", "0"))
        labelled.nodes[node]["label"] = f"{spider}:{phase}"
    for source, target, data in labelled.edges(data=True):
        labelled[source][target]["edge_label"] = str(data.get("type", ""))
    return weisfeiler_lehman_graph_hash(
        labelled, node_attr="label", edge_attr="edge_label", iterations=iterations
    )


def semantic_key(circuit, software=None):
    """Return a backend-agnostic semantic hash for ``circuit``, or ``None``.

    The hash is invariant under ZX-equivalent rewrites (``full_reduce``), so two
    functionally equivalent circuits share a key and a cached result.

    :param circuit: A Qibo ``Circuit`` or Qiskit ``QuantumCircuit``.
    :param software: ``"qibo"`` or ``"qiskit"``; auto-detected from the type when
        omitted.
    :return: A stable hash string, or ``None`` when the key cannot be computed
        (PyZX missing, an unsupported gate/backend) -- in which case the caller
        should simulate without caching. This function never raises.
    """
    if not pyzx_available:
        logger.warning(
            "PyZX is not installed; circuit caching is disabled "
            "(install the 'cache' extra: Qdislib[cache])."
        )
        return None
    if software is None:
        software = _detect_software(circuit)
    if software not in ("qibo", "qiskit", "cudaq", "pennylane"):
        logger.debug(
            "semantic_key is unsupported for software=%r; "
            "caching skipped for this circuit.",
            software,
        )
        return None
    try:
        qasm = _circuit_to_qasm(circuit, software)
        zx_circuit = zx.Circuit.from_qasm(qasm)
        graph = zx_circuit.to_graph()
        zx.simplify.full_reduce(graph)
        nx_graph = _zx_dict_to_nx(graph.to_dict(include_scalar=True))
        return _wl_hash_zx_graph(nx_graph)
    except Exception as exc:  # never let a caching helper break a computation
        logger.warning(
            "Could not compute a semantic key (%s); "
            "caching skipped for this circuit.",
            exc,
        )
        return None


# --------------------------------------------------------------------------- #
# Cache stores
# --------------------------------------------------------------------------- #


class CircuitCache:
    """Abstract semantic-cache store: a string key to a JSON-serialisable value.

    The value is a measurement-counts ``dict`` (counts mode) or an expectation
    ``float`` (expectation mode); the store treats it as opaque.
    """

    def get(self, key):
        """Return the cached value for ``key``, or ``None`` on a miss."""
        raise NotImplementedError

    def put(self, key, value):
        """Store ``value`` under ``key``."""
        raise NotImplementedError

    def delete(self, key):
        """Remove ``key`` if present (used to release a stale reservation)."""
        raise NotImplementedError

    def __contains__(self, key):
        return self.get(key) is not None

    def _record(self, value):
        """Update hit/miss counters from a ``get`` result and return the value."""
        if value is None:
            self.misses = getattr(self, "misses", 0) + 1
        else:
            self.hits = getattr(self, "hits", 0) + 1
        return value

    def stats(self):
        """Return cache statistics: ``hits``, ``misses``, ``hit_rate``, ``size``.

        Hit/miss counters are per store *instance*: under PyCOMPSs each worker
        gets its own pickled copy, so the counters reflect a single process only.
        ``size`` is still global for the shared LMDB/Redis/PSCO stores.
        """
        hits = getattr(self, "hits", 0)
        misses = getattr(self, "misses", 0)
        total = hits + misses
        try:
            size = len(self)
        except (TypeError, NotImplementedError):
            size = None
        return {
            "hits": hits,
            "misses": misses,
            "hit_rate": hits / total if total else 0.0,
            "size": size,
        }


class InMemoryCircuitCache(CircuitCache):
    """Process-local ``dict``-backed cache (serial / single-process use).

    Not shared between PyCOMPSs workers (each task receives its own copy), so use
    :class:`RedisPSCOCircuitCache` under PyCOMPSs.

    :param max_entries: Optional cap on the number of cached entries. Once the cap
        is reached, the least-recently-used entry is evicted (LRU) on the next
        insert, so a long run can't grow the cache without bound. The default
        (``None``) is unbounded, matching the previous behaviour.
    :param ttl: Optional time-to-live in seconds. An entry read after ``ttl``
        seconds is treated as a miss and dropped (lazy expiry), so a long-lived
        cache doesn't serve stale results, e.g. after an SDK upgrade. The default
        (``None``) never expires.
    """

    def __init__(self, max_entries=None, ttl=None):
        from collections import OrderedDict

        self._store = OrderedDict()  # key -> value (raw)
        self._expiry = {}  # key -> expires_at, only populated when ttl is set
        self._max_entries = max_entries
        self._ttl = ttl

    def get(self, key):
        if key not in self._store:
            return self._record(None)
        expires_at = self._expiry.get(key)
        if expires_at is not None and time.time() >= expires_at:
            self.delete(key)  # expired -> treat as a miss
            return self._record(None)
        self._store.move_to_end(key)  # mark most-recently-used (LRU)
        return self._record(self._store[key])

    def put(self, key, value):
        self._store[key] = value
        self._store.move_to_end(key)
        if self._ttl is not None:
            self._expiry[key] = time.time() + self._ttl
        if self._max_entries is not None:
            while len(self._store) > self._max_entries:
                evicted, _ = self._store.popitem(last=False)  # evict LRU
                self._expiry.pop(evicted, None)

    def delete(self, key):
        self._store.pop(key, None)
        self._expiry.pop(key, None)

    def __len__(self):
        return len(self._store)


# lmdb forbids opening the same environment path twice in one process, so share
# one open handle per path across all LMDBCircuitCache instances of this process.
_LMDB_ENVS = {}


class LMDBCircuitCache(CircuitCache):
    """Persistent on-disk LMDB cache, shareable across processes on one filesystem.

    Values are JSON-encoded. The environment is opened lazily and shared per path
    within a process; the instance holds only its configuration, so it pickles
    cleanly and a PyCOMPSs worker reopens the path on first use. Note that many
    concurrent writers contend on a shared filesystem; prefer
    :class:`RedisPSCOCircuitCache` for heavy distributed writing under PyCOMPSs.
    """

    def __init__(self, path, map_size=1 << 30, max_readers=500):
        if not lmdb_available:
            raise ImportError(
                "LMDBCircuitCache requires the 'lmdb' package "
                "(install the 'cache' extra: Qdislib[cache])."
            )
        self.path = str(path)
        self.map_size = map_size
        self.max_readers = max_readers

    def _env_handle(self):
        env = _LMDB_ENVS.get(self.path)
        if env is None:
            env = lmdb.open(
                self.path,
                map_size=self.map_size,
                subdir=True,
                max_dbs=1,
                max_readers=self.max_readers,
                readahead=False,
                meminit=False,
            )
            _LMDB_ENVS[self.path] = env
        return env

    def get(self, key):
        with self._env_handle().begin(write=False) as txn:
            raw = txn.get(key.encode("utf-8"))
        if raw is None:
            return self._record(None)
        try:
            return self._record(json.loads(raw.decode("utf-8")))
        except json.JSONDecodeError:
            return self._record(None)

    def put(self, key, value):
        raw = json.dumps(value).encode("utf-8")
        with self._env_handle().begin(write=True) as txn:
            txn.put(key.encode("utf-8"), raw)

    def delete(self, key):
        with self._env_handle().begin(write=True) as txn:
            txn.delete(key.encode("utf-8"))

    def __len__(self):
        with self._env_handle().begin() as txn:
            return txn.stat()["entries"]

    def close(self):
        """Close and release this path's shared LMDB environment."""
        env = _LMDB_ENVS.pop(self.path, None)
        if env is not None:
            env.close()


class RedisCircuitCache(CircuitCache):
    """Generic redis-py-backed cache for distributed sharing without COMPSs.

    Values are JSON-encoded under a key prefix. Use this when the chosen backend
    is Redis but the run is *not* under PyCOMPSs; under PyCOMPSs prefer
    :class:`RedisPSCOCircuitCache`, which uses the COMPSs storage layer. The
    client is created lazily and not pickled, so the store ships to a worker by
    its connection settings.
    """

    def __init__(self, host="localhost", port=6379, db=0, prefix="qdislib:", **kwargs):
        if not redis_available:
            raise ImportError(
                "RedisCircuitCache requires the 'redis' package "
                "(install the 'cache' extra: Qdislib[cache])."
            )
        self.host = host
        self.port = port
        self.db = db
        self.prefix = prefix
        self._kwargs = kwargs
        self._client = None

    def _conn(self):
        if self._client is None:
            self._client = _redis.Redis(
                host=self.host, port=self.port, db=self.db, **self._kwargs
            )
        return self._client

    def get(self, key):
        raw = self._conn().get(self.prefix + key)
        return self._record(None if raw is None else json.loads(raw))

    def put(self, key, value):
        self._conn().set(self.prefix + key, json.dumps(value))

    def delete(self, key):
        self._conn().delete(self.prefix + key)

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_client"] = None  # the live connection is not picklable
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


class _KeyValuePSCO(StorageObject):
    """@ClassFields"""

    # Mirrors the KeyValuePSCO from the published experiment's Redis variant: a
    # COMPSs persistent storage object holding one cached (key, value,
    # observables) entry, addressed by the semantic key.
    def __init__(self, key=None, value=None, observables=None):
        super(_KeyValuePSCO, self).__init__()
        self.key = key
        self.value = value
        self.observables = observables

    def get_value(self):
        """Return the stored value."""
        return self.value

    def set_value(self, value):
        """Set a new value (persistence handled externally)."""
        self.value = value

    def get_observables(self):
        """Return the stored observables."""
        return self.observables

    def set_observables(self, observables):
        """Set new observables (persistence handled externally)."""
        self.observables = observables


class RedisPSCOCircuitCache(CircuitCache):
    """BSC COMPSs PSCO/Redis cache, shared across PyCOMPSs tasks and nodes.

    Backed by the COMPSs ``storage`` API exactly as in the published experiment:
    each entry is a :class:`_KeyValuePSCO` persisted under its id (the semantic
    key), so a result written by one task is visible to every other task on any
    node. This is the right store under PyCOMPSs, where an in-memory ``dict`` is
    not shared and concurrent LMDB writers contend. Requires a COMPSs deployment
    with a Redis/PSCO storage backend configured.
    """

    def __init__(self):
        if not redis_psco_available:
            raise ImportError(
                "RedisPSCOCircuitCache requires the COMPSs 'storage' API, available "
                "in a COMPSs deployment with a Redis/PSCO storage backend."
            )

    def get(self, key):
        try:
            obj = get_by_ID(key)
        except Exception:
            # get_by_ID raises (e.g. AttributeError when it returns None
            # internally and then calls .pycompss_mark_as_unmodified()) on a
            # missing id; treat any lookup failure as a cache miss.
            return self._record(None)
        return self._record(None if obj is None else obj.get_value())

    def put(self, key, value, observables=None):
        try:
            obj = get_by_ID(key)
        except Exception:
            obj = None
        if obj is None:
            obj = _KeyValuePSCO(key, value, observables)
            obj.make_persistent(key)
        else:
            # The COMPSs Redis connector has no in-place update, so replace
            # (set_value -> delete_persistent -> make_persistent), as in the
            # experiment.
            obj.set_value(value)
            obj.delete_persistent()
            obj.make_persistent(key)

    def delete(self, key):
        try:
            obj = get_by_ID(key)
        except Exception:
            obj = None
        if obj is not None:
            obj.delete_persistent()


def default_circuit_cache(backend="memory", path=None, **kwargs):
    """Build a cache store for an explicitly chosen backend.

    The caller picks the store. ``"redis"`` then resolves by runtime: the BSC
    COMPSs PSCO store under PyCOMPSs (shared across tasks/nodes), or a generic
    redis-py store otherwise.

    :param backend: ``"memory"`` (default), ``"lmdb"`` or ``"redis"``.
    :param path: LMDB directory; required when ``backend="lmdb"``.
    :param kwargs: Forwarded to the chosen store (e.g. Redis ``host``/``port``/
        ``db``, or LMDB ``map_size``).
    :raises ValueError: For an unknown backend, or ``"lmdb"`` without a ``path``.
    :return: A :class:`CircuitCache` instance.
    """
    if backend == "memory":
        return InMemoryCircuitCache(**kwargs)
    if backend == "lmdb":
        if path is None:
            raise ValueError("default_circuit_cache(backend='lmdb') requires a 'path'.")
        return LMDBCircuitCache(path, **kwargs)
    if backend == "redis":
        # The caller chose Redis; use the COMPSs PSCO store under PyCOMPSs,
        # otherwise a generic redis-py store.
        if pycompss_available and redis_psco_available:
            return RedisPSCOCircuitCache()
        return RedisCircuitCache(**kwargs)
    raise ValueError(
        f"Unknown cache backend {backend!r}; expected 'memory', 'lmdb' or 'redis'."
    )


# --------------------------------------------------------------------------- #
# Cache strategies: what is stored per subcircuit, and how it is read back
# --------------------------------------------------------------------------- #


class CacheStrategy:
    """Policy for *what* the cache stores for a subcircuit and how to read it back.

    ``expectation_value(..., cache_mode=...)`` accepts the strings ``"expectation"``
    and ``"counts"`` (the two built-ins below) or any ``CacheStrategy`` instance,
    so you can store a custom payload -- e.g. counts plus metadata -- by
    overriding these three hooks. Each receives a *simulation context* ``sim``
    exposing the work to run on a cache miss:

    - ``sim.run_native()`` -- the backend evaluator's expectation value;
    - ``sim.counts()`` -- raw measurement counts for the circuit;
    - ``sim.expectation_from_counts(counts)`` -- expectation from those counts;
    - ``sim.observables`` / ``sim.software`` / ``sim.shots`` -- the context.

    Example -- cache counts together with metadata::

        class CountsWithMeta(CountsStrategy):
            def store(self, sim):
                counts = sim.counts()
                payload = {"counts": counts, "shots": sim.shots,
                           "observable": sim.observables}
                return payload, sim.expectation_from_counts(counts)

            def load(self, cached, sim):
                return sim.expectation_from_counts(cached["counts"])
    """

    def make_key(self, semantic_key, sim):
        """Return the cache key for this subcircuit.

        Defaults to ``semantic_key`` plus the observable, since one circuit
        topology yields different expectations under different observables.
        """
        return f"{semantic_key}|{sim.observables}"

    def store(self, sim):
        """On a cache miss: return ``(value_to_cache, expectation_to_return)``."""
        raise NotImplementedError

    def load(self, cached, sim):
        """On a cache hit: return the expectation derived from ``cached``."""
        raise NotImplementedError


class ExpectationStrategy(CacheStrategy):
    """Default: cache the scalar expectation, keyed by semantic hash + observable."""

    def store(self, sim):
        value = sim.run_native()
        return (None if value is None else float(value), value)

    def load(self, cached, sim):
        return cached


class CountsStrategy(CacheStrategy):
    """Cache raw counts keyed by the semantic hash alone; re-derive per observable.

    Matches the paper's granularity: one simulation's counts serve every
    observable measured on that circuit topology.
    """

    def make_key(self, semantic_key, sim):
        return semantic_key

    def store(self, sim):
        counts = sim.counts()
        return (counts, sim.expectation_from_counts(counts))

    def load(self, cached, sim):
        return sim.expectation_from_counts(cached)


# --------------------------------------------------------------------------- #
# Concurrency: placeholder reservation + deferred refresh (for HPC)
# --------------------------------------------------------------------------- #

# A reserved "computation in progress" marker. Under PyCOMPSs a task writes this
# under a key before simulating, so other tasks hitting the same key while it is
# in flight defer (return a PendingResult) instead of redundantly re-simulating;
# the real value overwrites it when the owning task finishes. Without it, many
# concurrent tasks on the same fragment would all miss and all simulate -- the
# hit rate collapses just when high concurrency should help most.
_PENDING = {"__qdislib_pending__": True}


def is_pending(value):
    """Return ``True`` if ``value`` is the reserved placeholder marker."""
    return isinstance(value, dict) and value.get("__qdislib_pending__") is True


class PendingResult:
    """A deferred subcircuit result whose key another task is still computing.

    Returned by the cache path when a lookup finds the placeholder. It is
    resolved at reconstruction time -- after the owning task has stored the real
    value, which COMPSs task dependencies guarantee has happened by then -- by
    re-reading the cache and applying the strategy. Resolution retries briefly to
    tolerate write lag and returns ``None`` (logged) if the value never appears
    (e.g. the owning task failed).
    """

    def __init__(
        self,
        cache,
        cache_key,
        strategy,
        observables,
        software,
        shots,
        retries=30,
        delay=0.5,
    ):
        self.cache = cache
        self.cache_key = cache_key
        self.strategy = strategy
        self.observables = observables
        self.software = software
        self.shots = shots
        self.retries = retries
        self.delay = delay

    def resolve(self):
        """Return the expectation once the owning task has stored the value."""
        import time

        from Qdislib.utils.backends import _SimContext

        sim = _SimContext(None, None, self.observables, self.software, self.shots)
        for attempt in range(self.retries):
            value = self.cache.get(self.cache_key)
            if value is not None and not is_pending(value):
                return self.strategy.load(value, sim)
            if attempt < self.retries - 1:
                time.sleep(self.delay)
        logger.warning(
            "Cache key %s never resolved (owning task may have failed); "
            "dropping this subcircuit's contribution.",
            self.cache_key,
        )
        return None


def _resolve_value(value):
    """Resolve a :class:`PendingResult` to its expectation; pass others through."""
    return value.resolve() if isinstance(value, PendingResult) else value
