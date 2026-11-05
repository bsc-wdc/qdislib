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

"""Hardware-QPU backends: the BSC cross-cluster bridge and IBM Quantum.

These wrap the existing execution paths (``_expec_value_qibo_qpu`` for the BSC
Blue/Red QPUs via the shared-filesystem bridge, ``_expec_value_qiskit_qpu`` for
IBM Quantum via qiskit-ibm-runtime) behind the :class:`BaseBackend` interface, so
the cutting algorithms select them uniformly instead of branching on ``qpu``.
The cutting code uses :meth:`evaluate_expectation` on the 1.x paths and
:meth:`evaluate_joint_counts` (from final counts on BSC, which has no mid-circuit
read-out; mid-circuit on IBM).
"""

from Qdislib.core.backend.base_backend import BaseBackend, BackendConfig


class BSCQpuBackend(BaseBackend):
    """A BSC QPU (Blue or Red) reached through the shared-filesystem bridge.

    Delegates to :func:`~Qdislib.core.cutting_algorithms.gate_cutting._expec_value_qibo_qpu`
    (which ships the subcircuit via :mod:`Qdislib.utils.qpu_bridge` and reads back
    the result). Configure the bridge with the ``QDISLIB_*`` / ``RUNCARD_*`` env
    vars documented in :mod:`Qdislib.utils.qpu_bridge`.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("BSC_QPU_BACKEND", config)

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="numpy", seed=None
    ):
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qibo_qpu

        subcirc = dag_to_circuit_qibo(dag, num_qubits)
        return _expec_value_qibo_qpu(subcirc, shots=shots, method=method)

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Final counts from a BSC QPU, through the shared-filesystem bridge
        (qubit-0-first). ``seed`` is ignored (hardware)."""
        return self._final_counts(dag, num_qubits)

    @staticmethod
    def _final_counts(dag, num_qubits):
        """Final counts of a DAG run through the bridge (``{}`` for an empty one)."""
        import qibo
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        from Qdislib.utils.qpu_bridge import submit_and_wait

        subcirc = dag_to_circuit_qibo(dag, num_qubits)
        if subcirc is None:
            return {}
        if type(subcirc) != qibo.models.circuit.Circuit:
            subcirc = subcirc[0]  # drop the "Observable I" marker; sampling reads bits
        return submit_and_wait(subcirc).counts()


#: IBM connections of this process, by device (a task's backend reconnects once).
_CONNECTIONS = {}


def _ibm_setting(extra: dict, key: str, env: str, default=""):
    """Read an IBM setting from explicit config, falling back to the env var."""
    import os

    value = extra.get(key)
    if value is not None:
        return value
    return os.environ.get(env, default)


class IBMBackend(BaseBackend):
    """IBM Quantum hardware via qiskit-ibm-runtime.

    Preferred setup is **explicit configuration** -- no environment variables
    required::

        QdBackEnd({"IBM_BACKEND": {
            "num_qubits": 127,
            "token": "<token>",
            "channel": "ibm_cloud",      # or "ibm_quantum"
            "instance": "<instance>",    # ibm_cloud only
            "backend_name": "ibm_marrakesh",  # omit -> least-busy real backend
            "mode": "job",               # "batch" / "session" only if your plan allows
            "options": {"twirling": {"enable_gates": True}},  # SamplerV2's options
            "pack": True,                # several circuits in each job: see below
            "max_time": 60,              # batch / session only
        }})

    By default every subcircuit is its own plain job (``SamplerV2(mode=backend)``):
    some plans allow no Batch or Session. ``mode="batch"`` or ``mode="session"``
    groups the jobs in a runtime Batch or Session instead.

    With ``pack=True`` (or a dict of :func:`~Qdislib.pack` options) the circuits of a
    group run side by side in one device-wide circuit, so ``circuits_per_task=`` many
    subcircuits take one job. Choose ``circuits_per_task`` near the number of groups of
    qubits the device has for them (``qd.pack(...).regions``).

    Under PyCOMPSs this backend goes to a task as its device's name and settings, not
    as a connection or a token (task arguments are files in the job's directory):
    each worker connects again, reading ``IBM_QUANTUM_TOKEN`` from its environment. An
    opened Batch or Session cannot go to a task.

    Each ``BackendConfig.extra`` key falls back to its ``IBM_QUANTUM_*`` env var
    if omitted, so existing env-var-based setups keep working. You may also hand
    in an already-authenticated runtime object via ``extra["service"]`` (a
    ``QiskitRuntimeService``) or a live ``extra["batch"]`` + ``extra["backend"]``
    (the pair the cutting code threads to its per-subcircuit tasks; in job mode
    ``batch`` is the backend itself).

    The runtime connection and backend selection (the "check") run **once**, in
    the constructor -- not per subcircuit.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("IBM_BACKEND", config)
        extra = config.extra
        self._options = extra.get("options")
        self._pack = extra.get("pack")
        # Task side: a session + backend were established upstream and threaded in.
        if extra.get("batch") is not None and extra.get("backend") is not None:
            self._session = extra["batch"]
            self._ibm_backend = extra["backend"]
        else:
            self._session, self._ibm_backend = self._connect(extra)

    @staticmethod
    def _connect(extra: dict):
        """Authenticate and select the backend (named or least-busy): ``(mode,
        backend)``, where ``mode`` is what the jobs run in -- the backend itself
        (plain jobs, the default) or an opened Batch / Session."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService
        except ImportError as exc:
            raise ImportError(
                "IBMBackend requires qiskit-ibm-runtime (the `qiskit` extra)."
            ) from exc

        service = extra.get("service")
        if service is None:
            token = _ibm_setting(extra, "token", "IBM_QUANTUM_TOKEN")
            instance = _ibm_setting(extra, "instance", "IBM_QUANTUM_INSTANCE")
            channel = _ibm_setting(extra, "channel", "IBM_QUANTUM_CHANNEL")
            proxies = extra.get("proxies")
            if channel == "ibm_quantum":
                service = QiskitRuntimeService(
                    channel=channel, token=token, proxies=proxies
                )
            elif channel == "ibm_cloud":
                service = QiskitRuntimeService(
                    channel=channel, token=token, instance=instance
                )
            else:
                raise ValueError(
                    f"Unsupported IBM Runtime channel: {channel!r}. Set channel= and "
                    "token= (on a PyCOMPSs worker: IBM_QUANTUM_CHANNEL and "
                    "IBM_QUANTUM_TOKEN in its environment, e.g. with --env_script)"
                )

        backend_name = _ibm_setting(extra, "backend_name", "IBM_QUANTUM_QPU_NAME")
        backend = (
            service.backend(backend_name)
            if backend_name
            else service.least_busy(operational=True, simulator=False)
        )
        mode = str(_ibm_setting(extra, "mode", "IBM_QUANTUM_MODE", "job")).lower()
        if mode == "job":
            return backend, backend
        if mode not in ("batch", "session"):
            raise ValueError(f"IBM mode is 'job', 'batch' or 'session', got {mode!r}")
        from qiskit_ibm_runtime import Batch, Session

        max_time = int(_ibm_setting(extra, "max_time", "IBM_QUANTUM_MAX_TIME", 60))
        opened = (Batch if mode == "batch" else Session)(
            backend=backend, max_time=max_time
        )
        return opened, backend

    def __getstate__(self):
        if self._session is not self._ibm_backend:
            raise TypeError(
                "IBMBackend: an opened Batch or Session cannot go to a PyCOMPSs task; "
                "use plain jobs (the default)"
            )
        state = super().__getstate__()
        extra = {
            key: value
            for key, value in self.config.extra.items()
            if key not in ("service", "batch", "backend", "token")
        }
        if self._ibm_backend is not None:
            extra["backend_name"] = self._ibm_backend.name
        state["config"] = BackendConfig(self.config.num_qubits, extra=extra)
        state["_session"] = state["_ibm_backend"] = None
        return state

    def _connected(self):
        """Connect, if this copy reached a task without a connection (once per
        device in a worker's process)."""
        if self._ibm_backend is None:
            extra = self.config.extra
            key = tuple(
                str(extra.get(k)) for k in ("backend_name", "channel", "instance")
            )
            if key not in _CONNECTIONS:
                _CONNECTIONS[key] = self._connect(extra)
            self._session, self._ibm_backend = _CONNECTIONS[key]
        return self._session, self._ibm_backend

    @property
    def packs(self):
        return bool(self._pack)

    @property
    def session(self):
        """What the jobs run in: the backend (plain jobs) or the opened Batch /
        Session (the cutting code threads this to its tasks)."""
        return self._connected()[0]

    @property
    def ibm_backend(self):
        """The selected IBM backend."""
        return self._connected()[1]

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """Joint final/cut-point counts on the IBM device (the cut-point
        measurements run mid-circuit). ``seed`` is ignored (hardware)."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _joint_counts_qiskit_qpu,
        )

        return _joint_counts_qiskit_qpu(
            dag,
            num_qubits,
            sorted(cut_point_qubits),
            shots,
            self.session,
            self.ibm_backend,
            self._options,
        )

    def evaluate_joint_counts_many(self, requests):
        """With ``pack``: every circuit of the same shots side by side, in as few
        plain jobs as the device's groups of qubits allow; submitted, then read."""
        if not self._pack or len(requests) < 2:
            return super().evaluate_joint_counts_many(requests)
        from Qdislib.core.cutting_algorithms import gate_cutting as gc
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.utils.packing import pack

        options = self._pack if isinstance(self._pack, dict) else {}
        out = [None] * len(requests)
        by_shots = {}
        for n, (dag, num_qubits, cuts, shots, _) in enumerate(requests):
            built = dag_to_circuit_qiskit(dag, num_qubits)
            if built is not None:
                base = built[0] if isinstance(built, (list, tuple)) else built
                by_shots.setdefault(shots, []).append((n, base, num_qubits, cuts))
        for shots, group in by_shots.items():
            circuits = []
            for _, base, _, _ in group:
                measured = base.copy()
                measured.measure_all()
                circuits.append(measured)
            packed = pack(circuits, self.ibm_backend, **options)
            sampler = gc.Sampler(mode=self.session, options=self._options)
            jobs = [sampler.run([c], shots=shots) for c in packed.circuits]
            counts = packed.unpack([job.result()[0] for job in jobs])
            for (n, base, num_qubits, cuts), one in zip(group, counts):
                out[n] = gc._joint_from_qiskit_counts(
                    one, base, num_qubits, sorted(cuts)
                )
        return out

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Final counts on the IBM device (qubit-0-first). ``seed`` is ignored."""
        return {
            "".join(str(b) for b in final): count
            for (final, _), count in self.evaluate_joint_counts(
                dag, num_qubits, (), shots
            ).items()
        }

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic", seed=None
    ):
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit_qpu

        subcirc = dag_to_circuit_qiskit(dag, num_qubits)
        return _expec_value_qiskit_qpu(
            subcirc,
            shots=shots,
            batch=self.session,
            backend=self.ibm_backend,
            prep_qubits=list(prep_qubits) or None,
            options=self._options,
        )
