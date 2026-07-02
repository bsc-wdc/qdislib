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
The cutting code uses :meth:`evaluate_expectation` (both backends) and, for the
sampling reconstruction, :meth:`evaluate_joint_counts` (BSC, for terminal cuts);
the generic sample/estimator hooks are left unimplemented.
"""

from Qdislib.core.backend.base_backend import BaseBackend, BackendConfig


class _QpuBackendBase(BaseBackend):
    """Shared boilerplate: the generic hooks are not used on the cut QPU path."""

    def transpile(self, circuit):
        raise NotImplementedError(f"{self.name}: use evaluate_expectation")

    def _execute_samples(self, *args, **kwargs):
        raise NotImplementedError(f"{self.name}: use evaluate_expectation")

    def _execute_estimator(self, *args, **kwargs):
        raise NotImplementedError(f"{self.name}: use evaluate_expectation")


class BSCQpuBackend(_QpuBackendBase):
    """A BSC QPU (Blue or Red) reached through the shared-filesystem bridge.

    Delegates to :func:`~Qdislib.core.cutting_algorithms.gate_cutting._expec_value_qibo_qpu`
    (which ships the subcircuit via :mod:`Qdislib.utils.qpu_bridge` and reads back
    the result). Configure the bridge with the ``QDISLIB_*`` / ``RUNCARD_*`` env
    vars documented in :mod:`Qdislib.utils.qpu_bridge`.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("BSC_QPU_BACKEND", config)

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="numpy"
    ):
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qibo_qpu

        subcirc = dag_to_circuit_qibo(dag, num_qubits)
        return _expec_value_qibo_qpu(subcirc, shots=shots, method=method)

    def evaluate_joint_counts(self, dag, num_qubits, cut_point_qubits=(), shots=1024):
        """Joint final/cut-point counts from a BSC QPU (sampling reconstruction).

        Ships the subcircuit through the shared-filesystem bridge and reads back the
        hardware counts (``result.counts()``, qubit-0-first). A real QPU reports only
        aggregated *final* counts -- the mid-circuit cut-point outcome isn't returned
        -- so each cut qubit's cut-point bit is taken to equal its final bit. That is
        exact for a **terminal** cut (no gate on a cut qubit after the cut, as in the
        image-cutting tutorial); a cut qubit that keeps evolving afterwards cannot be
        decoupled from aggregated hardware counts (the same limitation that guards the
        QPU expectation path).
        """
        import qibo
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        from Qdislib.utils.qpu_bridge import submit_and_wait

        subcirc = dag_to_circuit_qibo(dag, num_qubits)
        if subcirc is None:
            return None
        if type(subcirc) != qibo.models.circuit.Circuit:
            subcirc = subcirc[0]  # drop the "Observable I" marker; sampling reads bits

        freq = submit_and_wait(subcirc).counts()
        cut_points = list(cut_point_qubits)
        joint = {}
        for key, count in freq.items():
            final = tuple(int(key[q]) for q in range(num_qubits))
            cut = tuple(
                final[q] for q in cut_points
            )  # terminal cut: cut-point == final
            joint[(final, cut)] = joint.get((final, cut), 0) + count
        return joint


def _ibm_setting(extra: dict, key: str, env: str, default=""):
    """Read an IBM setting from explicit config, falling back to the env var."""
    import os

    value = extra.get(key)
    if value is not None:
        return value
    return os.environ.get(env, default)


class IBMBackend(_QpuBackendBase):
    """IBM Quantum hardware via qiskit-ibm-runtime.

    Preferred setup is **explicit configuration** -- no environment variables
    required::

        QdBackEnd({"IBM_BACKEND": {
            "num_qubits": 127,
            "token": "<token>",
            "channel": "ibm_cloud",      # or "ibm_quantum"
            "instance": "<instance>",    # ibm_cloud only
            "backend_name": "ibm_marrakesh",  # omit -> least-busy real backend
            "max_time": 60,
        }})

    Each ``BackendConfig.extra`` key falls back to its ``IBM_QUANTUM_*`` env var
    if omitted, so existing env-var-based setups keep working. You may also hand
    in an already-authenticated runtime object via ``extra["service"]`` (a
    ``QiskitRuntimeService``) or a live ``extra["batch"]`` + ``extra["backend"]``
    (the pair the cutting code threads to its per-subcircuit tasks).

    The runtime connection and backend selection (the "check") run **once**, in
    the constructor -- not per subcircuit.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("IBM_BACKEND", config)
        extra = config.extra
        # Task side: a session + backend were established upstream and threaded in.
        if extra.get("batch") is not None and extra.get("backend") is not None:
            self._session = extra["batch"]
            self._ibm_backend = extra["backend"]
        else:
            self._session, self._ibm_backend = self._connect(extra)

    @staticmethod
    def _connect(extra: dict):
        """Authenticate, select the backend (named or least-busy), open a Batch."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService, Batch
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
                raise ValueError(f"Unsupported IBM Runtime channel: {channel!r}")

        backend_name = _ibm_setting(extra, "backend_name", "IBM_QUANTUM_QPU_NAME")
        backend = (
            service.backend(backend_name)
            if backend_name
            else service.least_busy(operational=True, simulator=False)
        )
        max_time = int(_ibm_setting(extra, "max_time", "IBM_QUANTUM_MAX_TIME", 60))
        return Batch(backend=backend, max_time=max_time), backend

    @property
    def session(self):
        """The open runtime Batch/Session (the cutting code threads this to its tasks)."""
        return self._session

    @property
    def ibm_backend(self):
        """The selected IBM backend."""
        return self._ibm_backend

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic"
    ):
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit_qpu

        subcirc = dag_to_circuit_qiskit(dag, num_qubits)
        return _expec_value_qiskit_qpu(
            subcirc, shots=shots, batch=self._session, backend=self._ibm_backend
        )
