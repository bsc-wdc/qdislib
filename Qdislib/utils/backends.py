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

"""
Qdislib backend-agnostic helpers.

Utilities to detect which quantum SDK (Qiskit, Qibo or CUDA-Q) a circuit comes
from and to convert it into Qdislib's internal DAG (``networkx``) representation.
Keeping this dispatch in one place avoids duplicating the type-to-converter map
across the cutting algorithms.
"""

import logging

from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag
from Qdislib.utils.graph_cudaq import circuit_cudaq_to_dag
from Qdislib.utils.runtime import task, compss_wait_on, pycompss_available
from Qdislib.utils.circuit_cache import (
    semantic_key,
    default_circuit_cache,
    is_pending,
    PendingResult,
    _PENDING,
)

# Each SDK is optional. When one is absent its circuit class is replaced by a
# placeholder sentinel so the type dispatch below stays valid (no real object is
# ever an instance of the sentinel, so it simply never matches).
try:
    from qiskit.circuit import QuantumCircuit
except ImportError:
    QuantumCircuit = type("QuantumCircuit", (), {})
try:
    from qibo.models import Circuit as QiboCircuit
except ImportError:
    QiboCircuit = type("QiboCircuit", (), {})
try:
    from cudaq.kernel.kernel_decorator import PyKernelDecorator
except ImportError:
    PyKernelDecorator = type("PyKernelDecorator", (), {})

logger = logging.getLogger(__name__)

__all__ = [
    "convert_to_dag",
    "iter_gates",
    "num_qubits",
    "expectation_value",
    "normalize_observables",
    "Subcircuit",
    "CutSubcircuits",
]

#: Mapping from a circuit's SDK class to the converter that turns it into a DAG.
_CIRCUIT_TO_DAG = {
    QuantumCircuit: circuit_qiskit_to_dag,
    QiboCircuit: circuit_qibo_to_dag,
    PyKernelDecorator: circuit_cudaq_to_dag,
}


def convert_to_dag(circuit, observables=None):
    """Convert a Qiskit/Qibo/CUDA-Q circuit into a Qdislib DAG.

    Dispatches on the circuit type. If ``circuit`` is already a DAG (or any
    unrecognised type) it is returned unchanged.

    ``CX``/``CNOT`` gates are kept as-is here (not rewritten to ``H·CZ·H``): the
    DAG must stay in step with ``find_cut``, which names cut nodes after the
    gates the user actually wrote, so a user never sees a phantom ``CZ`` cut on a
    circuit that has none. Gate cutting handles the ``CX = H·CZ·H`` identity
    itself, at the cut site (see ``_generate_cut``).

    :param circuit: A Qiskit ``QuantumCircuit``, a Qibo ``Circuit``, a CUDA-Q
        kernel, or an already-built ``networkx`` DAG.
    :param observables: Optional observable string; when it contains ``"I"`` the
        identity observable nodes are appended to the resulting DAG.
    :return: The corresponding ``networkx`` DAG, or ``circuit`` unchanged if the
        type is not recognised.
    """
    converter = next(
        (fn for t, fn in _CIRCUIT_TO_DAG.items() if isinstance(circuit, t)), None
    )
    if converter is None:
        logger.warning("Unsupported circuit type: %s", type(circuit).__name__)
        return circuit  # Unknown type (e.g. already a DAG), return unchanged.

    kwargs = {}
    if observables and "I" in observables:
        kwargs["obs_I"] = observables

    return converter(circuit, **kwargs)


def normalize_observables(observables):
    """Normalize an observable specification into ``[(coefficient, pauli_string), ...]``.

    Lets the cutting routines accept a **weighted sum of Pauli operators** (a
    Hamiltonian), not just a single per-qubit Pauli string, since
    ``<sum_i c_i P_i> = sum_i c_i <P_i>``. Accepted forms:

    - a single Pauli string, e.g. ``"ZZI"`` -> ``[(1.0, "ZZI")]``;
    - a mapping ``{pauli_string: coefficient}``, e.g. ``{"ZZ": 0.5, "XI": -1.0}``;
    - an iterable of ``(coefficient, pauli_string)`` pairs, e.g.
      ``[(0.5, "ZZ"), (-1.0, "XI")]``.

    :param observables: A Pauli string, mapping, or iterable of weighted terms.
    :return: A list of ``(float, str)`` ``(coefficient, pauli_string)`` terms.
    :raises ValueError: If a Pauli string contains a character other than I/X/Y/Z,
        or if the terms have inconsistent lengths.
    """
    if isinstance(observables, str):
        terms = [(1.0, observables)]
    elif isinstance(observables, dict):
        terms = [(float(coeff), pauli) for pauli, coeff in observables.items()]
    else:
        terms = [(float(coeff), str(pauli)) for coeff, pauli in observables]

    if not terms:
        raise ValueError("Empty observable specification.")
    lengths = {len(pauli) for _, pauli in terms}
    for _, pauli in terms:
        if any(ch not in "IXYZ" for ch in pauli):
            raise ValueError(
                f"Invalid Pauli string {pauli!r}: only 'I', 'X', 'Y', 'Z' are allowed."
            )
    if len(lengths) > 1:
        raise ValueError(
            f"All Pauli terms must act on the same number of qubits; got lengths {sorted(lengths)}."
        )
    return terms


def change_basis_dag(dag, observables):
    """Append measurement-basis-rotation gates to a circuit DAG (backend-neutral).

    To measure a non-Z Pauli, the qubit is rotated into the Z basis before
    measurement: ``X`` -> append ``H``; ``Y`` -> append ``SDG`` then ``H``. ``Z``
    and ``I`` qubits are left unchanged. The rotation is appended at the qubit's
    current sink (its last gate), so a subsequent Z-basis measurement yields the
    X/Y expectation value.

    This operates on the DAG rather than the source circuit so it works for any
    input backend; the SDK circuits (Qiskit ``QuantumCircuit``, CUDA-Q kernels)
    cannot be appended to generically the way a Qibo ``Circuit`` can.

    :param dag: A circuit DAG produced by :func:`convert_to_dag`.
    :param observables: Per-qubit Pauli string (e.g. ``"XZY"``).
    :return: The same DAG with rotation nodes appended (mutated in place).
    """
    if not observables:
        return dag

    rotations = {"X": ("H",), "Y": ("SDG", "H")}
    for qubit, pauli in enumerate(observables):
        chain = rotations.get(pauli)
        if not chain:
            continue
        for offset, gate_name in enumerate(chain):
            node_id = f"BASIS{pauli}{qubit}_{offset}"
            dag.add_node(node_id, gate=gate_name, qubits=(qubit,), parameters=())
            # Reproduce the dependency edges that ``circuit_*_to_dag`` would build
            # for this gate: a single blue edge to the nearest predecessor on the
            # qubit, and red edges to every earlier gate on that qubit. The red
            # edges matter for wire cutting's component generation.
            for pred in reversed(list(dag.nodes)):
                if pred != node_id and qubit in (dag.nodes[pred].get("qubits") or ()):
                    dag.add_edge(pred, node_id, color="blue")
                    break
            for pred in reversed(list(dag.nodes)):
                if (
                    pred != node_id
                    and qubit in (dag.nodes[pred].get("qubits") or ())
                    and not dag.has_edge(pred, node_id)
                ):
                    dag.add_edge(pred, node_id, color="red")

    return dag


def iter_gates(circuit):
    """Yield each gate of a circuit as a backend-neutral ``(name, qubits)`` pair.

    This normalises the very different gate APIs of the supported SDKs into a
    single ordered list, so callers don't have to special-case the backend:

    - ``name`` is the upper-cased gate name following each SDK's own convention,
      kept consistent with ``circuit_{qibo,qiskit,cudaq}_to_dag`` so the
      ``f"{name}_{index}"`` identifiers line up (e.g. Qibo ``CNOT`` vs Qiskit
      ``CX``).
    - ``qubits`` is a tuple of integer qubit indices the gate acts on.

    One entry is produced per operation, in circuit order, matching the 1-based
    indexing used by the DAG converters.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :raises TypeError: If the circuit type is not recognised.
    :return: List of ``(name, qubits)`` tuples.
    """
    if isinstance(circuit, QiboCircuit):
        return [
            (gate.__class__.__name__.upper(), tuple(gate.qubits))
            for gate in circuit.queue
        ]
    if isinstance(circuit, QuantumCircuit):
        return [
            (
                instruction.operation.name.upper(),
                tuple(q._index for q in instruction.qubits),
            )
            for instruction in circuit.data
        ]
    if isinstance(circuit, PyKernelDecorator):
        # CUDA-Q has no gate-iteration API; reuse the source parser that
        # circuit_cudaq_to_dag already relies on. Node insertion order is gate
        # order, so this stays consistent with that converter's indexing.
        dag = circuit_cudaq_to_dag(circuit)
        return [
            (data["gate"].upper(), tuple(data["qubits"]))
            for _, data in dag.nodes(data=True)
        ]
    raise TypeError(
        f"Unsupported circuit type for gate iteration: {type(circuit).__name__}"
    )


def num_qubits(circuit):
    """Return the number of qubits of a Qiskit/Qibo/CUDA-Q circuit.

    Uses the SDK's native attribute when available (so idle qubits are counted);
    for CUDA-Q kernels, which expose no such attribute, the count is inferred
    from the highest qubit index touched by any gate.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :return: The number of qubits as an ``int``.
    """
    if isinstance(circuit, QiboCircuit):
        return circuit.nqubits
    if isinstance(circuit, QuantumCircuit):
        return circuit.num_qubits
    indices = [qubit for _, qubits in iter_gates(circuit) for qubit in qubits]
    return max(indices) + 1 if indices else 0


def _materialize(dag, num_qubits, software):
    """Build the native backend circuit for a subcircuit DAG.

    Returns ``[circuit, observables]`` where ``circuit`` is a Qibo ``Circuit``, a
    Qiskit ``QuantumCircuit`` or a CUDA-Q kernel. Used to produce the serial
    subcircuits and to rebuild a :class:`Subcircuit` in the process that runs it.
    """
    from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
    from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
    from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel
    from Qdislib.utils.exceptions import UnsupportedBackendError

    builders = {
        "qibo": dag_to_circuit_qibo,
        "qiskit": dag_to_circuit_qiskit,
        "cudaq": dag_to_circuit_cudaq_kernel,
    }
    if software not in builders:
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(builders)}."
        )
    return builders[software](dag, num_qubits)


class Subcircuit:
    """A serializable cut subcircuit: a DAG plus the metadata to rebuild it.

    A native Qibo or CUDA-Q subcircuit cannot be serialized (the Qibo collapse
    measurement carries a symbol bound to its own circuit, and a CUDA-Q kernel is
    JIT-compiled from source), so it cannot be shipped to a remote PyCOMPSs task.
    This wrapper instead carries the backend-neutral ``dag`` (plain nodes/edges,
    fully picklable) and builds the native circuit **on demand, in the process
    that evaluates it** -- so nothing fragile ever crosses a task boundary.

    :func:`gate_cutting_subcircuits` / :func:`wire_cutting_subcircuits` return
    these under PyCOMPSs; in serial they return the native circuits directly.
    Both forms are accepted by :func:`expectation_value`.
    """

    def __init__(self, dag, num_qubits, software):
        self.dag = dag
        self.num_qubits = num_qubits
        self.software = software

    def to_circuit(self):
        """Materialize the native backend circuit ``[circuit, observables]``."""
        return _materialize(self.dag, self.num_qubits, self.software)

    def __repr__(self):
        return f"Subcircuit(software={self.software!r}, num_qubits={self.num_qubits})"


def make_subcircuit(dag, num_qubits, software):
    """Return the per-subcircuit object handed back by the ``*_subcircuits`` API.

    In serial this is the **native backend circuit** (convenient to inspect or
    draw). Under PyCOMPSs it is a serializable :class:`Subcircuit` wrapper, so it
    can be shipped to remote tasks and evaluated there.
    """
    from Qdislib.utils.exceptions import UnsupportedBackendError

    if software not in ("qibo", "qiskit", "cudaq"):
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            "Expected one of ['qibo', 'qiskit', 'cudaq']."
        )
    if pycompss_available:
        # The DAG/qubit count come from @task helpers (e.g. _update_qubits), so
        # under COMPSs they are futures; resolve them to concrete values before
        # wrapping. These dag-prep tasks are cheap -- the expensive circuit build
        # and execution still happen later, in the per-subcircuit eval task.
        dag = compss_wait_on(dag)
        num_qubits = compss_wait_on(num_qubits)
        return Subcircuit(dag, num_qubits, software)
    return _materialize(dag, num_qubits, software)


@task(returns=1)
def _make_subcircuit_task(component, software):
    """Prepare one subcircuit as a COMPSs task -- the three-step *creation* step,
    a separate task from the per-subcircuit evaluation so creation and execution
    are distinct, independently schedulable tasks (create -> evaluate, task -> task).

    Renumbers the raw component to a standalone ``0..n`` subcircuit and returns the
    object handed to :func:`expectation_value`: a serializable :class:`Subcircuit`
    under PyCOMPSs (so ``*_subcircuits`` returns futures -- ``compss_wait_on`` one to
    inspect it), a native circuit in serial. It runs on a worker with concrete
    values, so it wraps directly, without the master-side ``compss_wait_on`` that
    :func:`make_subcircuit` uses when its inputs are futures.
    """
    from Qdislib.utils.graph_qibo import _update_qubits

    s_new, highest_qubit = _update_qubits(component)
    if pycompss_available:
        return Subcircuit(s_new, highest_qubit, software)
    return _materialize(s_new, highest_qubit, software)


def _run_native(circuit, software, shots):
    """Evaluate an already-built native subcircuit with the matching evaluator."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expec_value_qibo,
        _expec_value_cudaq,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit

    if software == "qibo":
        return _expec_value_qibo(circuit, shots)
    if software == "qiskit":
        return _expec_value_qiskit(circuit, shots)
    return _expec_value_cudaq(circuit, shots)


def _unpack_native(subcircuit):
    """Split a native subcircuit into ``(circuit, observables, software)``.

    Accepts a bare circuit or a ``[circuit, observables]`` pair and infers the
    backend from the circuit type. ``observables`` defaults to an all-``Z`` string
    of the right width; ``software`` is ``"cudaq"`` for anything not recognised as
    Qibo/Qiskit (a built kernel or a legacy ``(dag, num_qubits)`` payload).
    """
    payload = subcircuit
    observables = None
    if (
        isinstance(subcircuit, (list, tuple))
        and len(subcircuit) == 2
        and isinstance(subcircuit[1], str)
    ):
        payload, observables = subcircuit[0], subcircuit[1]

    if isinstance(payload, QiboCircuit):
        software = "qibo"
        if observables is None:
            observables = "Z" * payload.nqubits
    elif isinstance(payload, QuantumCircuit):
        software = "qiskit"
        if observables is None:
            observables = "Z" * payload.num_qubits
    else:
        software = "cudaq"
    return payload, observables, software


def _simulate_counts(circuit, software, shots):
    """Simulate a native circuit and return its raw measurement counts."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _counts_qibo,
        _counts_cudaq,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _counts_qiskit

    if software == "qibo":
        return _counts_qibo(circuit, shots)
    if software == "qiskit":
        return _counts_qiskit(circuit, shots)
    return _counts_cudaq(circuit, shots)


def _expectation_from_counts(counts, observables, software, shots):
    """Derive an expectation value from cached counts for the given backend."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expectation_from_counts_qibo,
        _expectation_from_counts_cudaq,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _expectation_from_counts_qiskit,
    )

    if software == "qibo":
        return _expectation_from_counts_qibo(counts, observables, shots)
    if software == "qiskit":
        return _expectation_from_counts_qiskit(counts, observables, shots)
    return _expectation_from_counts_cudaq(counts, observables, shots)


class _SimContext:
    """A subcircuit plus the simulation primitives a ``CacheStrategy`` may call.

    Passed to the strategy hooks so they can decide what to compute on a miss
    (the expectation directly, raw counts, or a custom payload) without knowing
    the backend.
    """

    def __init__(self, circuit, subcircuit, observables, software, shots):
        self.circuit = circuit
        self.subcircuit = subcircuit
        self.observables = observables
        self.software = software
        self.shots = shots

    def run_native(self):
        """Expectation value via the backend evaluator (the non-cached path)."""
        return _run_native(self.subcircuit, self.software, self.shots)

    def counts(self):
        """Raw measurement counts for the circuit."""
        return _simulate_counts(self.circuit, self.software, self.shots)

    def expectation_from_counts(self, counts):
        """Expectation value derived from ``counts`` for this observable."""
        return _expectation_from_counts(
            counts, self.observables, self.software, self.shots
        )


def _resolve_strategy(cache_mode):
    """Turn a ``cache_mode`` (string or ``CacheStrategy``) into a strategy."""
    from Qdislib.utils.circuit_cache import (
        CacheStrategy,
        CountsStrategy,
        ExpectationStrategy,
    )

    if isinstance(cache_mode, CacheStrategy):
        return cache_mode
    if cache_mode == "expectation":
        return ExpectationStrategy()
    if cache_mode == "counts":
        return CountsStrategy()
    raise ValueError(
        f"Unknown cache_mode {cache_mode!r}; expected 'expectation', 'counts', "
        "or a CacheStrategy instance."
    )


def _safe_get(cache, key):
    """Read from the cache, treating any backend error as a miss."""
    try:
        return cache.get(key)
    except Exception as exc:  # a cache fault must never break the computation
        logger.warning("Cache lookup failed for key %s (%s).", key, exc)
        return None


def _safe_put(cache, key, value):
    """Store in the cache, swallowing (and logging) any backend error."""
    try:
        cache.put(key, value)
    except Exception as exc:  # a cache fault must never break the computation
        logger.warning("Cache store failed for key %s (%s).", key, exc)


def _safe_delete(cache, key):
    """Remove a key from the cache, swallowing (and logging) any backend error."""
    try:
        cache.delete(key)
    except Exception as exc:
        logger.warning("Cache delete failed for key %s (%s).", key, exc)


def _cut_point_fingerprint(circuit, software):
    """Fingerprint the mid-circuit *collapse* measurements a gate cut inserts.

    A gate cut measures its cut qubit mid-circuit and lets it keep evolving; that
    collapse outcome carries the quasiprobability reconstruction sign (recovered by
    ``_expec_value_*``). The ZX/QASM ``semantic_key`` hashes only the *unitary*, so
    it drops the collapse -- two subcircuits with the same unitary but opposite
    cut-point sign would otherwise share a cache entry (a wrong hit that flips the
    reconstructed value). This returns a string identifying the collapse
    measurements so those subcircuits get distinct keys.

    Only **gate cutting** produces such measurements: wire-cut fragments differ in
    real (unitary) prepare/basis gates -- already captured by ``semantic_key`` -- and
    uncut circuits have none, so both yield an empty fingerprint and an unchanged key.
    """
    try:
        if software == "qibo":
            import qibo.gates as qg

            marks = sorted(
                tuple(gate.qubits)
                for gate in circuit.queue
                if isinstance(gate, qg.M) and getattr(gate, "collapse", False)
            )
            return ";".join(str(m) for m in marks)
        if software == "qiskit":
            # A gate-cut Qiskit subcircuit measures the cut qubit into its own
            # register mid-circuit (before the final measure_all is appended at run
            # time), so any measurement present here is a cut point.
            marks = sorted(
                tuple(circuit.find_bit(q).index for q in inst.qubits)
                for inst in circuit.data
                if inst.operation.name == "measure"
            )
            return ";".join(str(m) for m in marks)
    except Exception:  # a fingerprint helper must never break the computation
        return ""
    return ""


def _cached_expectation(subcircuit, shots, cache, cache_mode):
    """Get-or-simulate one subcircuit's expectation value through ``cache``.

    Computes the semantic key from the circuit and, on a hit, returns the stored
    result without simulating. The ``cache_mode`` strategy decides what is stored
    and how a hit is turned back into an expectation (see
    :class:`~Qdislib.utils.circuit_cache.CacheStrategy`); the built-in
    ``"expectation"`` and ``"counts"`` modes, or any custom ``CacheStrategy``.

    Under PyCOMPSs a miss first reserves the key with a placeholder so other
    concurrent tasks on the same fragment defer (returning a
    :class:`~Qdislib.utils.circuit_cache.PendingResult`, resolved at
    reconstruction) instead of all re-simulating -- without this, the hit rate
    collapses under high concurrency. Caching never changes the result: a circuit
    with no computable key (PyZX absent, an unparsable gate, a CUDA-Q kernel that
    won't translate) or a store error falls back to a direct simulation.
    """
    circuit, observables, software = _unpack_native(subcircuit)

    key = semantic_key(circuit, software)
    if key is None:
        return _run_native(subcircuit, software, shots)

    # A gate cut's mid-circuit collapse measurement carries the reconstruction sign
    # but is invisible to the (unitary) semantic hash; fold it into the key so two
    # subcircuits with the same unitary but opposite cut-point sign don't collide.
    # Empty for wire cutting / uncut circuits (no collapse) -> key unchanged.
    fingerprint = _cut_point_fingerprint(circuit, software)
    if fingerprint:
        key = f"{key}#{fingerprint}"

    strategy = _resolve_strategy(cache_mode)
    sim = _SimContext(circuit, subcircuit, observables, software, shots)
    cache_key = strategy.make_key(key, sim)

    cached = _safe_get(cache, cache_key)
    if cached is not None:
        if is_pending(cached):
            # Another task reserved this key and is still computing it; defer and
            # resolve at reconstruction, when COMPSs guarantees it has finished.
            return PendingResult(
                cache, cache_key, strategy, observables, software, shots
            )
        try:
            return strategy.load(cached, sim)
        except Exception as exc:  # corrupt/incompatible entry -> recompute
            logger.warning(
                "Cached value for %s unusable (%s); recomputing.", cache_key, exc
            )

    # Miss. Under PyCOMPSs, reserve the key with a placeholder first so other
    # concurrent tasks on the same fragment defer instead of all re-simulating.
    reserved = pycompss_available
    if reserved:
        _safe_put(cache, cache_key, _PENDING)
    try:
        value_to_store, result = strategy.store(sim)
    except Exception:
        if reserved:
            _safe_delete(cache, cache_key)  # release the reservation for a retry
        raise
    _safe_put(cache, cache_key, value_to_store)
    return result


@task(returns=1)
def _evaluate_subcircuit(subcircuit, shots, cache=None, cache_mode="expectation"):
    """Build and evaluate a :class:`Subcircuit` in a single task.

    Materialising and running in one task keeps the fragile native circuit
    (Qibo symbol / CUDA-Q kernel) inside one process; only the picklable DAG
    crossed the task boundary to get here. When a ``cache`` is supplied the
    lookup happens here too, on the worker, so a hit skips the build+simulate.
    """
    native = subcircuit.to_circuit()
    if cache is None:
        return _run_native(native, subcircuit.software, shots)
    return _cached_expectation(native, shots, cache, cache_mode)


def expectation_value(subcircuit, shots=1024, cache=None, cache_mode="expectation"):
    """Evaluate one cut subcircuit's expectation value, dispatching on its type.

    Convenience wrapper for the recommended workflow::

        subcircuits = qd.gate_cutting_subcircuits(circuit, cut, software="qibo")
        values = [qd.expectation_value(s) for s in subcircuits]
        result = subcircuits.reconstruct(values)

    Accepts both forms produced by the ``*_subcircuits`` functions: a serializable
    :class:`Subcircuit` (returned under PyCOMPSs), which is built and run in a
    single task, or a native backend circuit (returned in serial).

    Pass ``cache`` to skip re-simulating semantically-equivalent subcircuits via
    the ZX/Weisfeiler-Lehman cache (see :mod:`Qdislib.utils.circuit_cache`)::

        cache = qd.default_circuit_cache()      # PSCO/Redis under PyCOMPSs
        values = [qd.expectation_value(s, cache=cache) for s in subcircuits]

    :param subcircuit: A :class:`Subcircuit`, a native circuit, or
        ``[circuit, observables]``.
    :param shots: Number of measurement shots.
    :param cache: A :class:`~Qdislib.utils.circuit_cache.CircuitCache` to enable
        semantic caching, or ``True`` to use
        :func:`~Qdislib.utils.circuit_cache.default_circuit_cache`. ``None`` (the
        default) disables caching. Caching never changes the result.
    :param cache_mode: ``"expectation"`` (default) caches per-subcircuit
        expectations keyed by semantic hash + observable; ``"counts"`` caches raw
        measurement counts keyed by the semantic hash and re-derives the
        expectation (maximal reuse, matching the paper).
    :return: The subcircuit's expectation value (a PyCOMPSs future under COMPSs).
    """
    if subcircuit is None:  # empty component produced by a cut; contributes nothing
        return None

    if cache is True:
        cache = default_circuit_cache()

    if isinstance(subcircuit, Subcircuit):
        return _evaluate_subcircuit(subcircuit, shots, cache, cache_mode)

    # Under PyCOMPSs the ``*_subcircuits`` items are *futures* of Subcircuits
    # (creation is a task) -- a future can't be isinstance-checked on the master.
    # Native circuits only occur in serial, so under COMPSs anything that is not a
    # concrete native circuit is such a future: hand it to the eval task, which
    # resolves it (create -> evaluate, task -> task).
    if pycompss_available and not isinstance(
        subcircuit, (QiboCircuit, QuantumCircuit, list, tuple)
    ):
        return _evaluate_subcircuit(subcircuit, shots, cache, cache_mode)

    if cache is not None:
        return _cached_expectation(subcircuit, shots, cache, cache_mode)

    # Native circuit (serial path): detect the backend from the payload type.
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expec_value_qibo,
        _expec_value_cudaq,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit

    payload = subcircuit
    if (
        isinstance(subcircuit, (list, tuple))
        and len(subcircuit) == 2
        and isinstance(subcircuit[1], str)
    ):
        payload = subcircuit[0]

    if isinstance(payload, QiboCircuit):
        return _expec_value_qibo(subcircuit, shots)
    if isinstance(payload, QuantumCircuit):
        return _expec_value_qiskit(subcircuit, shots)
    # CUDA-Q: a built kernel or a legacy ``(dag, num_qubits)`` payload.
    return _expec_value_cudaq(subcircuit, shots)


class CutSubcircuits(list):
    """The subcircuits produced by a cut, that also know how to be recombined.

    This behaves exactly like the plain ``list`` of subcircuits returned before
    (iterate, index, ``len`` as usual), but additionally records the metadata
    needed to recombine the per-subcircuit expectation values and exposes a
    convenient :meth:`reconstruct`::

        subcircuits = qd.wire_cutting_subcircuits(circuit, cut, software="qibo")
        values = [qd.expectation_value(s) for s in subcircuits]
        result = subcircuits.reconstruct(values)        # no manual bookkeeping

    The standalone ``*_subcircuit_reconstruction`` functions remain available for
    callers that prefer to pass ``number_cuts``/``number_components`` explicitly.
    """

    def __init__(self, subcircuits, *, number_cuts, number_components, reconstructor):
        super().__init__(subcircuits)
        #: Number of (single-wire / single-gate) cuts applied.
        self.number_cuts = number_cuts
        #: Number of independent components each cut configuration splits into.
        self.number_components = number_components
        self._reconstructor = reconstructor

    def reconstruct(self, expected_values):
        """Recombine per-subcircuit expectation values into the original one.

        :param expected_values: Expectation values, one per subcircuit, in the
            same order as the subcircuits in this object.
        :return: The reconstructed expectation value of the uncut circuit.
        """
        return self._reconstructor(
            expected_values, self.number_cuts, self.number_components
        )
