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

"""Backend abstraction for Qdislib.

The cutting algorithms are backend-agnostic: they produce subcircuits and ask a
backend to *sample* them (counts) or *estimate* an observable (expectation
value). A backend is anything implementing :class:`BaseBackend` -- a simulator
(Aer, Qibo, CUDA-Q) or real hardware (IQM, IBM, a BSC QPU). New backends are
registered in :mod:`~Qdislib.core.backend.registry` and can be added by users
without touching the algorithms.
"""

import threading
import uuid
import logging
from abc import ABC, abstractmethod
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Literal

from Qdislib.core.backend.circuit_converter import CircuitConverter

#: A circuit in any of the supported source formats (Qiskit / Qibo / a CUDA-Q
#: kernel or its ``(dag, num_qubits)`` payload). Annotated as ``Any`` so importing
#: the backend layer never pulls in a specific SDK -- ``import Qdislib`` must work
#: with only one backend installed (Qibo-only / Qiskit-only / CUDA-Q-only).
SupportedCircuit = Any
logger = logging.getLogger(__name__)


def marginalize_identity_qubits(counts: dict, observables) -> dict:
    """Marginalize out any ``"I"``-observed qubit positions from a sampling result.

    ``dag_to_circuit_*`` (e.g. :func:`~Qdislib.utils.graph_qibo.dag_to_circuit_qibo`,
    and its Qiskit/CUDA-Q/PennyLane counterparts) returns ``[circuit, observables]``,
    where ``observables`` is a per-qubit ``"Z"``/``"I"`` string -- an ``"I"`` marks a
    qubit that :func:`~Qdislib.utils.backends.change_basis_dag` inserted an
    Observable-I terminal node for (the user's observable was identity there, so its
    measured value shouldn't affect the result). ``evaluate_expectation`` honours
    this by construction (it's handed the whole ``[circuit, observables]`` pair). A
    raw-sampling ``evaluate_samples`` must call this explicitly on its result, since
    a bare ``{bitstring: count}`` dict has nowhere to carry which positions are
    identity-observed: forcing those positions to a constant character and merging
    the resulting duplicate keys recovers exactly the marginal distribution over the
    remaining bits (the "don't care about this qubit's value" semantics of an
    identity observable) -- equivalent to XOR-ing a constant into a downstream
    Z-parity computation, i.e. a no-op there, which is what "identity" means.

    :param counts: ``{qubit-0-first bitstring: count}``, before masking.
    :param observables: The per-qubit ``"Z"``/``"I"`` string from ``dag_to_circuit_*``,
        or ``None``/empty (no-op) when the source has no observable metadata (e.g. a
        bare native circuit with no ``[circuit, observables]`` wrapping).
    :return: ``counts`` unchanged if there's no ``"I"`` in ``observables``, else the
        marginalized ``{bitstring: count}`` dict (same key width, duplicates merged).
    """
    if not observables or "I" not in observables:
        return counts
    positions = {i for i, ch in enumerate(observables) if ch == "I"}
    merged: dict = {}
    for key, count in counts.items():
        masked = "".join("0" if i in positions else ch for i, ch in enumerate(key))
        merged[masked] = merged.get(masked, 0) + count
    return merged


@dataclass(init=False)
class BackendConfig:
    """Static configuration for a backend: its width and its own options.

    A backend's options are named directly -- they are what you actually
    configure, so they need no wrapper::

        BackendConfig(url=..., token=..., quantum_computer="garnet")   # IQM
        BackendConfig(use_gpu=True)                                    # Aer
        BackendConfig(platform="numpy")                                # Qibo

    Each backend documents the ones it reads (and ignores the rest), which is why
    they aren't declared here. They are collected into :attr:`extra`, so a caller
    holding a ready-made dict can still pass ``extra={...}`` instead.

    ``num_qubits`` is the width the backend is allowed to use. Omit it and the
    backend resolves its own -- a real device reports what it physically has, and
    anything that can't resolve one falls back to ``DEFAULT_NUM_QUBITS`` in
    :func:`~Qdislib.core.backend.registry.build_backend`. State it only to *cap* a
    backend below its natural width.
    """

    num_qubits: int = 0
    #: The backend's own options, however they were passed.
    extra: dict = field(default_factory=dict)

    def __init__(self, num_qubits: int = 0, extra: dict | None = None, **options):
        self.num_qubits = num_qubits
        self.extra = {**(extra or {}), **options}


@dataclass
class Job:
    """A unit of work tracked in a backend's FIFO queue."""

    job_id: str
    kind: Literal["samples", "estimator"]
    n_qubits: int
    shots: int | None
    submitted_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))


@dataclass
class SampleJobResult:
    """Result of a sampling run: measurement counts plus provenance metadata."""

    counts: dict[str, int]
    backend_name: str
    shots: int
    job_id: str
    metadata: dict = field(default_factory=dict)


@dataclass
class ExpectationJobResult:
    """Result of an estimator run: an expectation value plus provenance metadata."""

    expectation_value: float
    shots: int
    backend_name: str
    job_id: str
    obs: str = field(default="")
    metadata: dict = field(default_factory=dict)


class BaseBackend(ABC):
    """Common interface for all Qdislib backends.

    FIFO queue management is centralised here via the Template Method pattern:
    subclasses implement ``_execute_samples`` / ``_execute_estimator`` / ``transpile``
    (the hooks), while ``run_samples`` / ``run_estimator`` orchestrate the full
    enqueue -> execute -> dequeue lifecycle. The queue length is what
    :class:`~Qdislib.core.backend.unified.QdBackEnd` uses to load-balance.
    """

    def __init__(self, name: str, config: BackendConfig) -> None:
        self.name = name
        self.config = config
        self._queue: deque[Job] = deque()
        self._lock = threading.Lock()  # guards queue access in multi-threaded contexts

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def num_qubits(self) -> int:
        return self.config.num_qubits

    @property
    def queue_length(self) -> int:
        with self._lock:
            return len(self._queue)

    @property
    def queue_snapshot(self) -> list[Job]:
        """Return an immutable copy of the current queue (useful for logging / monitoring)."""
        with self._lock:
            return list(self._queue)

    def can_run(self, circuit: SupportedCircuit) -> bool:
        if not self.num_qubits:  # width never resolved: don't claim to know
            return True
        return CircuitConverter.num_qubits(circuit) <= self.num_qubits

    # ------------------------------------------------------------------
    # Public API -- manage the queue and delegate to the abstract hooks
    # ------------------------------------------------------------------

    def run_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int
    ) -> SampleJobResult:
        """Sample the circuit, updating the FIFO queue before and after execution."""
        job = self._enqueue(
            kind="samples",
            n_qubits=CircuitConverter.num_qubits(circuit),
            shots=shots,
        )
        try:
            return self._execute_samples(circuit, transpile, shots, job)
        finally:
            self._dequeue(job)

    def run_estimator(
        self, circuit: SupportedCircuit, transpile: bool, observable: Any, shots: int
    ) -> ExpectationJobResult:
        """Compute the expectation value of an observable, updating the FIFO queue."""
        job = self._enqueue(
            kind="estimator",
            n_qubits=CircuitConverter.num_qubits(circuit),
            shots=shots,
        )
        try:
            return self._execute_estimator(circuit, transpile, shots, observable, job)
        finally:
            self._dequeue(job)

    # ------------------------------------------------------------------
    # Cutting hook -- evaluate a Qdislib subcircuit given as a backend-neutral DAG
    # ------------------------------------------------------------------

    def evaluate_expectation(
        self,
        dag: Any,
        num_qubits: int,
        prep_qubits: Any = (),
        shots: int = 1024,
        method: str = "automatic",
        seed: int = None,
    ) -> float:
        """Evaluate the expectation of a Qdislib subcircuit DAG on this backend.

        This is the seam the cutting algorithms delegate to instead of branching on
        the backend name: they hand over the backend-neutral DAG (its Pauli
        observable rides along as ``Observable I`` markers) and the set of
        ``prep_qubits`` (cut qubits measured for the decomposition, whose sign is
        recovered from the backend's own mid-circuit outcome). The backend builds
        its native circuit, executes it, and reconstructs the value. Backends that
        only support :meth:`run_samples` / :meth:`run_estimator` may leave this
        unimplemented.

        :param dag: A renumbered connected-component DAG of the cut circuit.
        :param num_qubits: Qubit count of the component.
        :param prep_qubits: Renumbered cut qubits measured for the decomposition.
        :param shots: Number of measurement shots.
        :param method: Backend-specific simulation method hint.
        :param seed: Seed for a simulator's sampling (ignored by hardware).
        :return: The subcircuit's expectation value.
        """
        raise NotImplementedError(
            f"{self.name} does not implement evaluate_expectation (circuit cutting)"
        )

    def evaluate_samples(
        self, dag: Any, num_qubits: int, shots: int = 1024, seed: int = None
    ) -> dict[str, int]:
        """Sample a Qdislib subcircuit DAG and return counts in a standard convention.

        Used by the sampling reconstruction of gate and wire cutting (including the
        optimized fragment-separated wire-cutting path). Every backend returns keys
        as **qubit-0-first** bitstrings (``key[i]`` is qubit ``i``), regardless of the
        SDK's native endianness, so the reconstruction combines components uniformly.

        Implementations: ``dag_to_circuit_*`` returns ``[circuit, observables]`` when
        ``dag`` carries Observable-I markers (an ``"I"`` per-qubit observable
        component); pass that ``observables`` string through
        :func:`marginalize_identity_qubits` on the raw counts before returning, so an
        identity-observed qubit doesn't leak into a caller's Z-parity extraction as
        if it were a real ``"Z"`` observation.

        :param dag: A renumbered connected-component DAG of the cut circuit.
        :param num_qubits: Qubit count of the component.
        :param shots: Number of measurement shots.
        :return: ``{qubit-0-first bitstring: count}``, with any Observable-I
            positions already marginalized.
        """
        raise NotImplementedError(
            f"{self.name} does not implement evaluate_samples (sampling reconstruction)"
        )

    def evaluate_joint_counts(
        self,
        dag: Any,
        num_qubits: int,
        cut_point_qubits: Any = (),
        shots: int = 1024,
        seed: int = None,
    ) -> dict:
        """Sample a gate-cut subcircuit, capturing the JOINT final/cut-point outcomes.

        The seam the *sampling* reconstruction delegates to -- the distribution
        counterpart of :meth:`evaluate_expectation`. Gate-cut sampling needs, per
        shot, the final measurement **and** the mid-circuit cut-point outcomes of
        ``cut_point_qubits`` kept correlated, so a cut qubit that keeps evolving
        after the cut still reconstructs correctly (its final bit no longer equals
        its cut-point bit). Returns ``{(final_bits, cut_point_bits): count}`` with
        ``final_bits`` qubit-0-first over the component and ``cut_point_bits`` over
        ``cut_point_qubits`` in order.

        Backends whose execution yields only aggregated *final* counts (most real
        QPUs) cannot report the mid-circuit outcome and leave this unimplemented.

        :param dag: A renumbered connected-component DAG of the cut circuit.
        :param num_qubits: Qubit count of the component.
        :param cut_point_qubits: Renumbered cut qubits measured for the decomposition.
        :param shots: Number of measurement shots.
        :return: ``{(final_bits, cut_point_bits): count}``.
        """
        raise NotImplementedError(
            f"{self.name} does not implement evaluate_joint_counts (sampling reconstruction)"
        )

    # ------------------------------------------------------------------
    # Abstract hooks -- to be implemented by subclasses
    # ------------------------------------------------------------------

    @abstractmethod
    def transpile(self, circuit: SupportedCircuit) -> Any:
        """Transpile the circuit into the backend's native format."""

    @abstractmethod
    def _execute_samples(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        job: Job,
    ) -> SampleJobResult:
        """Backend-specific sampling logic. Queue management is handled by the base class."""

    @abstractmethod
    def _execute_estimator(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        observable: Any,
        job: Job,
    ) -> ExpectationJobResult:
        """Backend-specific estimator logic. Raise ``NotImplementedError`` if unsupported."""

    # ------------------------------------------------------------------
    # Internal queue management
    # ------------------------------------------------------------------

    def _enqueue(
        self,
        kind: Literal["samples", "estimator"],
        n_qubits: int,
        shots: int | None,
    ) -> Job:
        """Create a Job, append it to the FIFO queue, and return it."""
        job = Job(
            job_id=str(uuid.uuid4()),
            kind=kind,
            n_qubits=n_qubits,
            shots=shots,
        )
        with self._lock:
            self._queue.append(job)
        logger.debug(
            "[%s] Job %s enqueued (kind=%s, queue_len=%d)",
            self.name,
            job.job_id,
            kind,
            len(self._queue),
        )
        return job

    def _dequeue(self, job: Job) -> None:
        """Remove a completed (or failed) job from the queue."""
        with self._lock:
            try:
                self._queue.remove(job)
            except ValueError:
                logger.warning(
                    "[%s] Job %s not found in queue during dequeue",
                    self.name,
                    job.job_id,
                )
        logger.debug(
            "[%s] Job %s removed (queue_len=%d)",
            self.name,
            job.job_id,
            len(self._queue),
        )
