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

"""The result cache of the cutting engine: fragment circuits already run are not
run again.

Turned on per call (``cache=`` on every cutting function) or for a session
(:func:`use_cache`). What it stores is what the engine reads back from a fragment
circuit, whatever the backend: the joint final / cut-point counts of a shot run,
or the exact correlators -- so one entry serves every observable later read from
that circuit.

Two keys:

- ``"exact"`` (the default): the circuit's identity -- its fragment (gates,
  qubits, cut markers) and the terms it runs (its ``spec``) -- plus the backend
  and the shots. Computed on the master without building the circuit, so a hit
  never becomes a task, and it cannot match a different circuit.
- ``"semantic"``: also matches circuits that differ but are equivalent, by
  :func:`~Qdislib.utils.circuit_cache.semantic_key` (ZX), computed in the task.

A shot run's result is reused by later runs (the key holds the backend and the
shots, not the seed), which saves the most device time; ``cache_refresh=True``
runs it again and replaces it, e.g. to repeat a run for error bars. Exact results
are the same on every simulator and are shared between them.
"""

import contextlib
import hashlib
import logging
import os
import pickle

import numpy as np

logger = logging.getLogger(__name__)

__all__ = ["use_cache"]

#: Exact results are cached for fragments of up to this many qubits (their
#: correlator vector has ``2**n`` entries).
EXACT_MAX_QUBITS = 18

_KEYS = ("exact", "semantic")
_SESSION = {"store": None, "key": "exact"}
_STORES = {}  # the stores built from True / a string, one per spec


def _in_compss_runtime():
    """True when running inside a PyCOMPSs application (not just installed)."""
    try:
        from pycompss.util.context import in_pycompss

        return bool(in_pycompss())
    except Exception:  # noqa: BLE001
        try:
            from pycompss.util.context import CONTEXT

            return bool(CONTEXT.in_pycompss())
        except Exception:  # noqa: BLE001
            return False


def _store(spec):
    """A store from what ``cache=`` takes: ``True`` (the default store: in memory,
    or LMDB on the working directory's filesystem under PyCOMPSs), ``"memory"``,
    ``"lmdb:<path>"``, ``"redis"`` / ``"redis://host:port/db"``, or a
    :class:`~Qdislib.utils.circuit_cache.CircuitCache` object."""
    from Qdislib.utils import circuit_cache as cc

    if isinstance(spec, cc.CircuitCache):
        return spec
    if spec is True:
        if _in_compss_runtime():
            path = os.environ.get(
                "QDISLIB_CACHE_DIR", os.path.join(os.getcwd(), ".qdislib_cache")
            )
            spec = f"lmdb:{path}"
        else:
            spec = "memory"
    if not isinstance(spec, str):
        raise TypeError(
            "cache= takes True, 'memory', 'lmdb:<path>', 'redis[://host:port/db]' or "
            f"a CircuitCache object, not {type(spec).__name__}"
        )
    if spec not in _STORES:
        if spec == "memory":
            if _in_compss_runtime():
                raise ValueError(
                    "cache='memory' is one process's memory: PyCOMPSs tasks cannot "
                    "write to it. Use cache=True (LMDB), 'lmdb:<path>' or 'redis'."
                )
            _STORES[spec] = cc.InMemoryCircuitCache()
        elif spec.startswith("lmdb:"):
            path = spec[len("lmdb:") :]
            os.makedirs(path, exist_ok=True)
            _STORES[spec] = cc.LMDBCircuitCache(path)
        elif spec == "redis" or spec.startswith("redis://"):
            kwargs = {}
            if spec.startswith("redis://"):
                from urllib.parse import urlparse

                url = urlparse(spec)
                kwargs = {
                    "host": url.hostname or "localhost",
                    "port": url.port or 6379,
                    "db": int((url.path or "/0").strip("/") or 0),
                }
            _STORES[spec] = cc.default_circuit_cache("redis", **kwargs)
        else:
            raise ValueError(f"Unknown cache {spec!r}")
    return _STORES[spec]


@contextlib.contextmanager
def _restoring(previous):
    try:
        yield _SESSION["store"]
    finally:
        _SESSION.update(previous)


class use_cache:  # noqa: N801 -- used as a function and as a context manager
    """Turn the engine's result cache on for every cutting call of the session::

        qd.use_cache(True)                    # in memory (LMDB under PyCOMPSs)
        qd.use_cache("lmdb:/scratch/cache")   # persistent, shared by processes
        qd.use_cache(None)                    # off again

        with qd.use_cache("memory", key="semantic"):
            ...                               # only inside the block

    A ``cache=`` given to a call wins over it (``cache=False`` turns it off for
    that call).

    :param store: ``True``, ``"memory"``, ``"lmdb:<path>"``, ``"redis"``,
        ``"redis://host:port/db"``, a ``CircuitCache`` object, or ``None`` (off).
    :param key: ``"exact"`` (default) or ``"semantic"`` (see the module).
    """

    def __init__(self, store=True, *, key="exact"):
        if key not in _KEYS:
            raise ValueError(f"key must be one of {_KEYS}, got {key!r}")
        self._previous = dict(_SESSION)
        _SESSION.update(
            store=None if store in (None, False) else _store(store), key=key
        )

    def __enter__(self):
        return _SESSION["store"]

    def __exit__(self, *exc):
        _SESSION.update(self._previous)
        return False


def resolve(cache, cache_key=None, refresh=False):
    """The cache a call uses (``None``: off): its own ``cache=``, else the
    session's."""
    if cache is False:
        return None
    if cache is None:
        store, key = _SESSION["store"], _SESSION["key"]
    else:
        store, key = _store(cache), "exact"
    if store is None:
        return None
    key = cache_key or key
    if key not in _KEYS:
        raise ValueError(f"cache_key must be one of {_KEYS}, got {key!r}")
    return EngineCache(store, key, refresh)


# --------------------------------------------------------------------------- #
# Payloads: what a fragment circuit's result is stored as (JSON)
# --------------------------------------------------------------------------- #


def to_payload(result):
    """A fragment run's ``(corr, original)`` as a JSON value, or ``None`` when it
    is not cached (an exact result too wide to store)."""
    from Qdislib.core.cutting_algorithms.gate_cutting import _CountCorrelators

    corr, original = result
    original = [int(q) for q in original]
    if isinstance(corr, _CountCorrelators):
        return {
            "kind": "counts",
            "counts": [
                ["".join(map(str, final)), "".join(map(str, cut)), int(count)]
                for (final, cut), count in corr.joint.items()
            ],
            "prep": [int(q) for q in corr.prep_sorted],
            "k": int(corr.k),
            "original": original,
        }
    if isinstance(corr, dict):  # every correlator listed (an empty fragment)
        if len(corr) != 2 ** len(original):
            return None
        corr = [corr[m] for m in range(len(corr))]
    corr = np.asarray(corr, dtype=float).reshape(-1)
    if len(corr) > 2**EXACT_MAX_QUBITS:
        return None
    return {"kind": "exact", "corr": corr.tolist(), "original": original}


def from_payload(payload):
    """``(corr, original)`` back from :func:`to_payload`."""
    from Qdislib.core.cutting_algorithms.gate_cutting import _CountCorrelators

    if payload["kind"] == "counts":
        joint = {
            (tuple(int(b) for b in final), tuple(int(b) for b in cut)): count
            for final, cut, count in payload["counts"]
        }
        corr = _CountCorrelators(joint, payload["prep"], payload["k"])
    else:
        corr = np.asarray(payload["corr"], dtype=float)
    return corr, list(payload["original"])


# --------------------------------------------------------------------------- #
# Keys
# --------------------------------------------------------------------------- #


def _canonical(value):
    """A picklable, deterministic form of a nested value (floats as their exact
    bits, arrays as bytes): equal inputs give equal bytes."""
    if isinstance(value, np.ndarray):
        return ("nd", value.dtype.str, value.shape, value.tobytes())
    if isinstance(value, dict):
        return ("d", tuple(sorted((repr(k), _canonical(v)) for k, v in value.items())))
    if isinstance(value, (list, tuple)):
        return ("l", tuple(_canonical(v) for v in value))
    if isinstance(value, (set, frozenset)):
        return ("s", tuple(sorted(repr(v) for v in value)))
    if isinstance(value, (bool, int, str, type(None))):
        return value
    if isinstance(value, (float, complex, np.number)):
        return repr(complex(value) if isinstance(value, complex) else float(value))
    if hasattr(value, "__dict__"):
        return ("o", type(value).__name__, _canonical(vars(value)))
    return repr(value)


def template_fingerprint(template):
    """A hash of everything a fragment's circuits are built from: its gates, cut
    markers and edges, its cuts, and the term tables they use."""
    cached = getattr(template, "_fingerprint", None)
    if cached is not None:
        return cached
    graph = template.graph
    nodes = sorted((str(n), _canonical(dict(d))) for n, d in graph.nodes(data=True))
    edges = sorted((str(u), str(v)) for u, v in graph.edges())
    fragment = {k: v for k, v in template.fragment.items() if k != "nodes"}
    tables = [template.gate_tables[k - 1] for k, _ in template.fragment["gate_sides"]]
    groups = [template.groups[g] for g in template.fragment["had_groups"]]
    blob = pickle.dumps(
        _canonical(
            (
                nodes,
                edges,
                fragment,
                tables,
                groups,
                template.ancilla_base,
                template.cut_source_index,
                template.prepared_to_orig,
            )
        ),
        protocol=4,
    )
    template._fingerprint = hashlib.sha256(blob).hexdigest()
    return template._fingerprint


def backend_id(runner, shots):
    """What a result depends on besides the circuit: the backend (``"exact"`` for
    a statevector run, the same everywhere) and the shots."""
    if shots is None:
        return "exact"
    if isinstance(runner, str):
        name = runner
    else:
        config = getattr(runner, "config", None)
        extra = getattr(config, "extra", None) or {}
        name = ":".join(
            str(part)
            for part in (
                type(runner).__name__,
                getattr(runner, "name", ""),
                extra.get("backend_name", "") if isinstance(extra, dict) else "",
            )
        )
    return f"{name}|shots={int(shots)}"


def circuit_key(template, spec, runner, shots):
    """The exact key of one fragment circuit run (see the module)."""
    gate_config, states, setting = spec
    leg = template.leg(gate_config)
    blob = pickle.dumps(
        _canonical((leg, states, setting, backend_id(runner, shots))), protocol=4
    )
    return (
        f"qd-exact:{template_fingerprint(template)}:{hashlib.sha256(blob).hexdigest()}"
    )


class EngineCache:
    """One call's cache: the store, the key mode, and the circuits the call has
    already submitted (their futures are shared, so a circuit two terms of an
    observable need runs once, also with ``sync=False``, before any is stored)."""

    def __init__(self, store, key="exact", refresh=False):
        self.store = store
        self.key = key
        self.refresh = bool(refresh)
        #: key -> the (future) result of a circuit this call submitted.
        self._inflight = {}

    def run(self, circuit, shots, seed, runner=None):
        """The result of one fragment circuit: from the store or an earlier
        submission when it has one, else run (in a task that stores it).
        ``runner`` is where it runs (default: the circuit's own backend)."""
        from Qdislib.core.cutting_algorithms.mixed_cutting import _run_fragment

        runner = circuit.software if runner is None else runner
        key = circuit_key(circuit.template, circuit.spec, runner, shots)
        if key in self._inflight:  # this call already runs it (another term)
            return self._inflight[key]
        if not self.refresh:
            if self.key == "exact":
                payload = self._get(key)
                if payload is not None:
                    return from_payload(payload)
        request = CacheRequest(
            self.store, key if self.key == "exact" else None, self.refresh, shots
        )
        future = _run_fragment(
            circuit.template,
            circuit.spec,
            runner,
            shots,
            circuit.qpu,
            circuit.qpu_dict,
            seed,
            None,
            request,
        )
        self._inflight[key] = future
        return future

    def _get(self, key):
        from Qdislib.utils.circuit_cache import is_pending

        try:
            payload = self.store.get(key)
        except Exception as exc:  # noqa: BLE001 -- a store error is a miss
            logger.warning("cache read failed (%s); running the circuit", exc)
            return None
        if payload is None or is_pending(payload):
            return None
        return payload


class CacheRequest:
    """What a fragment task needs to use the cache: the store, the exact key
    (``None``: compute the semantic key in the task) and whether to refresh."""

    def __init__(self, store, key, refresh, shots):
        self.store = store
        self.key = key
        self.refresh = refresh
        self.shots = shots

    def lookup_key(self, dag, runner):
        """The key in the task: the exact one, or the circuit's semantic key
        (``None`` when it has none: then it is not cached)."""
        if self.key is not None:
            return self.key
        from Qdislib.utils.circuit_cache import semantic_key
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.utils.graph_qibo import _update_qubits

        renumbered, n = _update_qubits(dag.copy())
        if renumbered is None:
            return None
        built = dag_to_circuit_qiskit(renumbered, n)
        circuit = built[0] if isinstance(built, (list, tuple)) else built
        key = semantic_key(circuit, "qiskit")
        if key is None:
            return None
        return f"qd-semantic:{key}:{backend_id(runner, self.shots)}"

    def get(self, key):
        if self.refresh:
            return None
        try:
            payload = self.store.get(key)
        except Exception as exc:  # noqa: BLE001
            logger.warning("cache read failed (%s); running the circuit", exc)
            return None
        from Qdislib.utils.circuit_cache import is_pending

        return None if payload is None or is_pending(payload) else payload

    def put(self, key, result):
        payload = to_payload(result)
        if payload is None:
            return
        try:
            self.store.put(key, payload)
        except Exception as exc:  # noqa: BLE001 -- caching never breaks a run
            logger.warning("cache write failed (%s)", exc)
