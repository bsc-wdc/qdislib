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

from Qdislib.utils.sdks import sdk_names as _sdk_names
from Qdislib.utils.graph_qibo import _final_measurements, circuit_qibo_to_dag
from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag
from Qdislib.utils.graph_cudaq import circuit_cudaq_to_dag
from Qdislib.utils.graph_pennylane import circuit_pennylane_to_dag, PennylaneCircuit
from Qdislib.utils.runtime import task, compss_wait_on, pycompss_available
from Qdislib.utils.circuit_cache import (
    semantic_key,
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
try:
    from pennylane.tape import QuantumScript as PennylaneScript
    from pennylane import QNode as PennylaneQNode
except ImportError:
    PennylaneScript = type("PennylaneScript", (), {})
    PennylaneQNode = type("PennylaneQNode", (), {})

logger = logging.getLogger(__name__)

__all__ = [
    "convert_to_dag",
    "iter_gates",
    "num_qubits",
    "expectation_value",
    "normalize_observables",
    "Subcircuit",
    "CutSubcircuits",
    "CombinedCutSubcircuits",
]

#: Mapping from a circuit's SDK class to the converter that turns it into a DAG.


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
    from Qdislib.utils.sdks import sdk_of

    entry = sdk_of(circuit)
    if entry is None:
        import networkx

        if not isinstance(circuit, networkx.Graph):
            logger.warning("Unsupported circuit type: %s", type(circuit).__name__)
        return circuit  # already a DAG (or an unknown type): unchanged

    kwargs = {}
    if observables and "I" in observables:
        kwargs["obs_I"] = observables

    return entry.call("to_dag", circuit, **kwargs)


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
    indexing used by the DAG converters. What those converters skip while keeping
    the numbering -- a Qiskit ``barrier``, and a measurement that ends its qubits
    (see :func:`~Qdislib.utils.graph_qibo._final_measurements`) -- keeps its slot
    but acts on no qubit, so it is never offered as a cut.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :raises TypeError: If the circuit type is not recognised.
    :return: List of ``(name, qubits)`` tuples.
    """
    from Qdislib.utils.sdks import sdk_of

    entry = sdk_of(circuit)
    if entry is None:
        raise TypeError(
            f"Unsupported circuit type for gate iteration: {type(circuit).__name__}"
        )
    listed = entry.call("gates", circuit)
    final = _final_measurements([(m, qubits) for _, qubits, m in listed])
    return [
        (name, () if position in final else qubits)
        for position, (name, qubits, _) in enumerate(listed)
    ]


def num_qubits(circuit):
    """Return the number of qubits of a Qiskit/Qibo/CUDA-Q circuit.

    Uses the SDK's native attribute when available (so idle qubits are counted);
    for CUDA-Q kernels, which expose no such attribute, the count is inferred
    from the highest qubit index touched by any gate.

    :param circuit: A Qibo ``Circuit``, Qiskit ``QuantumCircuit`` or CUDA-Q kernel.
    :return: The number of qubits as an ``int``.
    """
    from Qdislib.utils.sdks import sdk_of

    entry = sdk_of(circuit)
    if entry is not None and entry.num_qubits:
        return entry.call("num_qubits", circuit)
    indices = [qubit for _, qubits in iter_gates(circuit) for qubit in qubits]
    return max(indices) + 1 if indices else 0


def _materialize(dag, num_qubits, software):
    """Build the native backend circuit for a subcircuit DAG.

    Returns ``[circuit, observables]`` where ``circuit`` is a Qibo ``Circuit``, a
    Qiskit ``QuantumCircuit`` or a CUDA-Q kernel. Used to produce the serial
    subcircuits and to rebuild a :class:`Subcircuit` in the process that runs it.
    """
    from Qdislib.utils.exceptions import UnsupportedBackendError
    from Qdislib.utils.sdks import sdk, sdk_names

    entry = sdk(software) if isinstance(software, str) else None
    if entry is None:
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(sdk_names())}."
        )
    if not entry.installed():
        # say so here, not as an AttributeError deep inside the converter
        raise ImportError(
            f"Building {software} circuits needs {entry.package}, which is not "
            f"installed in this Python. Install it with:  pip install "
            f"qdislib[{entry.name}]"
        )
    return entry.call("from_dag", dag, num_qubits)


def _first(built):
    """The circuit of a converter's result (some return ``[circuit, ...]``)."""
    return built[0] if isinstance(built, (list, tuple)) else built


class Subcircuit:
    """One subcircuit of a cut -- the same object from every ``*_subcircuits`` call.

    - ``circuit``: the native circuit in ``software`` (a qibo ``Circuit``, Qiskit
      ``QuantumCircuit``, CUDA-Q kernel or PennyLane circuit), built on access;
    - ``software``: the SDK it is written in;
    - ``num_qubits``: its width -- e.g. to pick a device per subcircuit;
    - ``observable``: the Z/I string read at its end, or ``None`` when every qubit
      is read in Z (distribution lists, and fragments read for several Paulis).

    Run it with :func:`expectation_value` (``*_subcircuits`` lists) or
    :func:`~Qdislib.sample_subcircuit` (``*_sampling_subcircuits`` lists), choosing
    ``backend`` there, and hand the results to the list's ``reconstruct``. It holds
    the backend-neutral DAG, so it pickles across PyCOMPSs tasks; the native circuit
    is only built in the process that needs it.

    ``s[0]`` / ``s[1]`` are the circuit and the observable, as the 1.x
    ``[circuit, observable]`` pairs were.
    """

    #: What the run step returns for it: "expectation" or "distribution".
    kind = "expectation"

    def __init__(self, dag, num_qubits, software):
        self.dag = dag
        self.num_qubits = num_qubits
        self.software = software

    def _renumbered(self):
        """``(dag, width)`` on qubits ``0..width-1`` (a copy)."""
        from Qdislib.utils.graph_qibo import _update_qubits

        return _update_qubits(self.dag.copy())

    def to_circuit(self):
        """``[circuit, observable]`` in ``software``. Built from a copy: the
        converters consume the DAG's ``Observable I`` markers."""
        return _materialize(self.dag.copy(), self.num_qubits, self.software)

    @property
    def circuit(self):
        """The native circuit (a new one on every access: a qibo circuit with
        mid-circuit measurements can only be run once)."""
        return _first(self.to_circuit())

    @property
    def observable(self):
        built = self.to_circuit()
        if isinstance(built, (list, tuple)) and len(built) > 1:
            if isinstance(built[1], str):
                return built[1]
        return "Z" * self.num_qubits

    def __getitem__(self, index):
        return [self.circuit, self.observable][index]

    def __iter__(self):
        return iter([self.circuit, self.observable])

    def __len__(self):
        return 2

    def _expectation(self, shots, backend, seed, cache, cache_mode):
        if self.dag is None:  # an empty component: contributes nothing
            return None
        if shots is None:
            return _exact_expectation_task(self, backend)
        if cache is not None:
            if backend is not None:
                raise ValueError("cache= runs on the local simulator: drop backend=")
            return _evaluate_subcircuit(self, shots, cache, cache_mode)
        return _run_subcircuit(
            self.dag,
            self.num_qubits,
            self.software if backend is None else backend,
            shots,
            seed,
        )

    def _exact_expectation(self, backend):
        """``<observable>`` from the statevector (mid-circuit measurements deferred)."""
        from Qdislib.core.cutting_algorithms.mixed_cutting import (
            _exact_fragment,
            _exact_runner,
        )

        runner = self.software if backend is None else backend
        if runner not in _sdk_names():
            raise ValueError("shots=None (exact) needs a simulator backend")
        corr, original = _exact_fragment(self.dag, _exact_runner(runner))
        observable = self.observable
        mask = sum(1 << i for i, q in enumerate(original) if observable[q] != "I")
        return float(corr[mask])

    def _sample(self, shots, backend, seed):
        raise TypeError(
            "This subcircuit comes from a *_subcircuits list (expectation value): run "
            "it with expectation_value(). sample_subcircuit() is for the lists that "
            "*_sampling_subcircuits return."
        )

    def __repr__(self):
        return (
            f"{type(self).__name__}(software={self.software!r}, "
            f"num_qubits={self.num_qubits})"
        )


class SampledSubcircuit(Subcircuit):
    """A subcircuit of a ``*_sampling_subcircuits`` list: run it with
    :func:`~Qdislib.sample_subcircuit`, which returns its Z correlators.

    ``dag`` keeps the circuit's own qubit numbering (the reconstruction needs it);
    ``circuit`` is renumbered to ``0..num_qubits-1``.
    """

    kind = "distribution"

    def __init__(self, dag, software):
        qubits = {q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()}
        super().__init__(dag, len(qubits), software)

    @property
    def component(self):  # the 1.x name
        return self.dag

    def to_circuit(self):
        dag, width = self._renumbered()
        return _materialize(dag, width, self.software)

    @property
    def observable(self):
        return None

    def _expectation(self, shots, backend, seed, cache, cache_mode):
        raise TypeError(
            "This subcircuit comes from a *_sampling_subcircuits list (distribution): "
            "run it with sample_subcircuit()."
        )

    def _sample(self, shots, backend, seed):
        if shots is None:  # statevector, mid-circuit measurements deferred
            runner = self.software if backend is None else backend
            if runner not in _sdk_names():
                raise ValueError("shots=None (exact) needs a simulator backend")
            return _exact_correlators_task(self.dag, runner)
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _sample_gate_subcircuit,
        )

        return _sample_gate_subcircuit(
            self.dag,
            self.software if backend is None else backend,
            shots,
            False,
            None,
            seed,
        )


@task(returns=1)
def _run_subcircuit(dag, num_qubits, backend, shots, seed):
    """Run one subcircuit's DAG on ``backend`` (an SDK, a registered backend name or
    a backend object) and return its expectation value."""
    from Qdislib.core.backend.registry import build_execution_backend
    from Qdislib.core.cutting_algorithms.gate_cutting import _cut_measurements

    executor = build_execution_backend(backend, num_qubits=num_qubits)
    dag = dag.copy()  # the converters consume its Observable I markers
    return executor.evaluate_expectation(
        dag, num_qubits, _cut_measurements(dag), shots, seed=seed
    )


def make_subcircuit(dag, num_qubits, software):
    """The :class:`Subcircuit` handed back by the ``*_subcircuits`` API."""
    from Qdislib.utils.exceptions import UnsupportedBackendError

    if software not in _sdk_names():
        raise UnsupportedBackendError(
            f"Unsupported software backend: '{software}'. "
            f"Expected one of {list(_sdk_names())}."
        )
    if pycompss_available:
        # The DAG/qubit count come from @task helpers (e.g. _update_qubits), so
        # under COMPSs they are futures; resolve them to concrete values before
        # wrapping.
        dag = compss_wait_on(dag)
        num_qubits = compss_wait_on(num_qubits)
    return Subcircuit(dag, num_qubits, software)


@task(returns=1)
def _make_subcircuit_task(component, software):
    """Prepare one subcircuit as a COMPSs task -- the three-step *creation* step,
    a separate task from the per-subcircuit evaluation (create -> evaluate, task ->
    task). Renumbers the raw component and returns its :class:`Subcircuit` (a future
    under PyCOMPSs: ``compss_wait_on`` one to inspect it)."""
    from Qdislib.utils.graph_qibo import _update_qubits

    s_new, highest_qubit = _update_qubits(component)
    return Subcircuit(s_new, highest_qubit, software)


def _run_native(circuit, software, shots):
    """Evaluate an already-built native subcircuit with the matching evaluator."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expec_value_qibo,
        _expec_value_cudaq,
        _expec_value_pennylane,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit

    if software == "qibo":
        return _expec_value_qibo(circuit, shots)
    if software == "qiskit":
        return _expec_value_qiskit(circuit, shots)
    if software == "pennylane":
        return _expec_value_pennylane(circuit, shots)
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
    elif isinstance(payload, PennylaneCircuit):
        software = "pennylane"
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
        _counts_pennylane,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import _counts_qiskit

    if software == "qibo":
        return _counts_qibo(circuit, shots)
    if software == "qiskit":
        return _counts_qiskit(circuit, shots)
    if software == "pennylane":
        return _counts_pennylane(circuit, shots)
    return _counts_cudaq(circuit, shots)


def _expectation_from_counts(counts, observables, software, shots):
    """Derive an expectation value from cached counts for the given backend."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expectation_from_counts_qibo,
        _expectation_from_counts_cudaq,
        _expectation_from_counts_pennylane,
    )
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _expectation_from_counts_qiskit,
    )

    if software == "qibo":
        return _expectation_from_counts_qibo(counts, observables, shots)
    if software == "qiskit":
        return _expectation_from_counts_qiskit(counts, observables, shots)
    if software == "pennylane":
        return _expectation_from_counts_pennylane(counts, observables, shots)
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
        if software == "pennylane":
            # A gate-cut PennylaneCircuit carries the collapse points as MEASURE ops.
            marks = sorted(
                tuple(qubits)
                for gate, qubits, _params in circuit.ops
                if gate.upper() == "MEASURE"
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
    if shots is None and cache is None:  # exact: its statevector path, inline here
        return subcircuit._expectation(None, None, None, None, cache_mode)
    native = subcircuit.to_circuit()
    if cache is None:
        return _run_native(native, subcircuit.software, shots)
    return _cached_expectation(native, shots, cache, cache_mode)


def expectation_value(
    subcircuit,
    shots=1024,
    *,
    backend=None,
    seed=None,
    cache=None,
    cache_mode="expectation",
):
    """Run one subcircuit of a ``*_subcircuits`` list (step two of three)::

        subs = qd.gate_cutting_subcircuits(circuit, cuts, "ZZZZ")
        values = [qd.expectation_value(s, shots=4096) for s in subs]
        result = subs.reconstruct(values)

    :param subcircuit: A :class:`Subcircuit` (or, as in 1.x, a native circuit or a
        ``[circuit, observable]`` pair).
    :param shots: Shots; ``None`` is exact for the subcircuits that support it.
    :param backend: Where to run it: an SDK (its local simulator), a registered
        backend name (``"IBM_BACKEND"``, ``"IQM_BACKEND"``, ...) or a
        :func:`~Qdislib.build_backend` object; ``None`` is the subcircuit's own SDK.
        Choosing per subcircuit routes them by size::

            dev = qd.build_backend("IQM_BACKEND", qd.BackendConfig(quantum_computer="garnet"))
            values = [
                qd.expectation_value(s, shots=4000, backend=dev if s.num_qubits <= 20 else None)
                for s in subs
            ]
    :param seed: Seed for a simulator's sampling.
    :param cache: Reuse a result already computed: what ``cache=`` takes on
        :func:`~Qdislib.circuit_cutting` (``True``, ``"memory"``, ``"lmdb:<path>"``,
        ``"redis[://host:port/db]"`` or a
        :class:`~Qdislib.utils.circuit_cache.CircuitCache`); ``None`` or ``False``
        is off.
    :param cache_mode: ``"expectation"``, ``"counts"`` or a
        :class:`~Qdislib.utils.circuit_cache.CacheStrategy`, for the 1.x native
        forms; a ``circuit_cutting`` fragment's cache keeps its whole result.
    :raises ValueError: For ``shots=None`` with a backend that is no simulator.
    :return: What the list's ``reconstruct`` takes for this subcircuit (its
        expectation value; Z correlators for ``circuit_cutting`` fragments).
    """
    if subcircuit is None:  # empty component produced by a cut; contributes nothing
        return None
    if cache is False:
        cache = None
    elif cache is True or isinstance(cache, str):
        # The shared store: a new one per call would never hold anything to reuse.
        from Qdislib.core.cutting_algorithms.engine_cache import _store

        cache = _store(cache)

    if isinstance(subcircuit, Subcircuit):
        return subcircuit._expectation(shots, backend, seed, cache, cache_mode)

    # Under PyCOMPSs the ``*_subcircuits`` items are *futures* of Subcircuits
    # (creation is a task) -- a future can't be isinstance-checked on the master.
    # Native circuits only occur in serial, so under COMPSs anything that is not a
    # concrete native circuit is such a future: hand it to a task that resolves it.
    if pycompss_available and not isinstance(
        subcircuit, (QiboCircuit, QuantumCircuit, list, tuple)
    ):
        if backend is None and seed is None:
            return _evaluate_subcircuit(subcircuit, shots, cache, cache_mode)
        return _evaluate_subcircuit_on(subcircuit, shots, backend, seed)

    # A native circuit (the 1.x serial form).
    if backend is not None:
        if hasattr(backend, "evaluate_subcircuit"):
            return backend.evaluate_subcircuit(subcircuit, shots=shots)
        raise TypeError(
            "backend= on a native circuit needs a backend object with "
            "evaluate_subcircuit; pass the Subcircuit from *_subcircuits instead."
        )
    if cache is not None:
        return _cached_expectation(subcircuit, shots, cache, cache_mode)

    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expec_value_qibo,
        _expec_value_cudaq,
        _expec_value_pennylane,
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
        return _expec_value_qiskit(subcircuit, shots, seed=seed)
    if isinstance(payload, PennylaneCircuit):
        return _expec_value_pennylane(subcircuit, shots, seed=seed)
    # CUDA-Q: a built kernel or a legacy ``(dag, num_qubits)`` payload.
    return _expec_value_cudaq(subcircuit, shots)


@task(returns=1)
def _exact_expectation_task(subcircuit, backend):
    """A subcircuit's exact ``<observable>`` (statevector), as a task."""
    return subcircuit._exact_expectation(backend)


@task(returns=1)
def _exact_correlators_task(dag, runner):
    """A sampling subcircuit's exact Z correlators (statevector), as a task."""
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        _exact_fragment,
        _exact_runner,
    )

    corr, original = _exact_fragment(dag, _exact_runner(runner))
    return dict(enumerate(corr)), original


@task(returns=1)
def _evaluate_subcircuit_on(subcircuit, shots, backend, seed):
    """Run a (future of a) :class:`Subcircuit` on ``backend`` inside one task."""
    from Qdislib.core.backend.registry import build_execution_backend
    from Qdislib.core.cutting_algorithms.gate_cutting import _cut_measurements

    if shots is None:  # exact: the subcircuit's own statevector path, inline here
        return subcircuit._expectation(None, backend, seed, None, "expectation")

    executor = build_execution_backend(
        subcircuit.software if backend is None else backend,
        num_qubits=subcircuit.num_qubits,
    )
    dag = subcircuit.dag.copy()  # the converters consume its Observable I markers
    return executor.evaluate_expectation(
        dag, subcircuit.num_qubits, _cut_measurements(dag), shots, seed=seed
    )


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


class CombinedCutSubcircuits(CutSubcircuits):
    """Several independently cut subcircuit sets, recombined into one value.

    Behaves as the flat list of every part's subcircuits, in order, so the usual
    three-step loop is unchanged. :meth:`reconstruct` hands each part its own slice
    of the values and combines the parts' results:

    - ``coefficients=None``: a **product**, for disconnected pieces of one circuit
      (each cut, or uncut, on its own): ``<A (x) B> = <A><B>``.
    - ``coefficients=[c_1, ...]``: a **weighted sum**, for the terms of a Pauli sum
      ``sum_k c_k P_k``, each term cut separately.

    Grouping the flat list by its length instead (``len // 6**n``) mixes up the
    pieces as soon as they are not all alike, e.g. one cut piece and one uncut.
    """

    def __init__(self, parts, coefficients=None):
        list.__init__(self, [sub for part in parts for sub in part])
        self.parts = list(parts)
        self.coefficients = None if coefficients is None else list(coefficients)
        if self.coefficients is not None and len(self.coefficients) != len(self.parts):
            raise ValueError("need one coefficient per part")
        #: Cuts applied: across all pieces, or per term for a weighted sum.
        counts = [part.number_cuts for part in self.parts]
        self.number_cuts = sum(counts) if self.coefficients is None else max(counts)
        #: Number of parts recombined.
        self.number_components = len(self.parts)

    @property
    def parameters(self):
        """The circuit's free parameters (names), as each part's."""
        return self.parts[0].parameters

    def bind(self, values):
        """Every part with the circuit's free parameters bound (see
        ``CutCircuits.bind``)."""
        return CombinedCutSubcircuits(
            [part.bind(values) for part in self.parts], self.coefficients
        )

    def reconstruct_counts(self, counts, bit_order="q0_first"):
        """``reconstruct`` from the counts of :meth:`circuits` run with your own
        tools (see ``CutCircuits.reconstruct_counts``)."""
        from Qdislib.core.cutting_algorithms.mixed_cutting import result_from_counts

        counts = list(counts)
        if len(counts) != len(self):
            raise ValueError(
                f"expected {len(self)} counts, one per subcircuit, got {len(counts)}"
            )
        return self.reconstruct(
            [result_from_counts(s, c, bit_order) for s, c in zip(self, counts)]
        )

    def circuits(self, software=None, *, defer_measurements=False):
        """The native circuits of every part, in order (see
        ``CutCircuits.circuits``)."""
        return [
            circuit
            for part in self.parts
            for circuit in part.circuits(
                software, defer_measurements=defer_measurements
            )
        ]

    def describe(self, values=None):
        """A table of every part's subcircuits (see ``CutCircuits.describe``)."""
        from Qdislib.core.cutting_algorithms.describe import describe

        return describe(self, values)

    def reconstruct(self, expected_values):
        """Recombine one value per subcircuit, in this object's order."""
        values = list(expected_values)
        if len(values) != len(self):
            raise ValueError(
                f"expected {len(self)} values, one per subcircuit, got {len(values)}"
            )
        results, start = [], 0
        for part in self.parts:
            results.append(part.reconstruct(values[start : start + len(part)]))
            start += len(part)
        if self.coefficients is None:
            product = 1.0
            for result in results:
                product *= result
            return product
        return sum(c * r for c, r in zip(self.coefficients, results))
