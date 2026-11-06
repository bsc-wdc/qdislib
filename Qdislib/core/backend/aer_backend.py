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

"""Qiskit Aer simulator backend (CPU or cuStateVec GPU)."""

import logging
from typing import Any

from qiskit import QuantumCircuit, transpile as qiskit_transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import EstimatorV2 as AerEstimator

from Qdislib.core.backend.base_backend import (
    BaseBackend,
    BackendConfig,
    Job,
    SupportedCircuit,
    SampleJobResult,
    ExpectationJobResult,
)
from Qdislib.core.backend.circuit_converter import CircuitConverter

logger = logging.getLogger(__name__)

# cuStateVec is the recommended statevector backend for NVIDIA GPUs.
_GPU_METHOD = "statevector"
_GPU_DEVICE = "GPU"
_GPU_CUSTATEVEC_FLAG = True
_CPU_METHOD = "automatic"


class AerBackend(BaseBackend):
    """Backend wrapping Qiskit Aer's ``AerSimulator``.

    Configuration keys (``BackendConfig.extra``)
    ---------------------------------------------
    use_gpu : bool, optional
        Run on GPU via cuStateVec (needs ``qiskit-aer-gpu`` + an NVIDIA GPU).
    noise_model : qiskit_aer.noise.NoiseModel, optional
        A pre-built Aer noise model.
    precision : {"single", "double"}, optional
        Statevector precision (default ``"double"``).

    A cut runs with these options. On a node with no GPU, ``use_gpu`` fails when the
    first circuit runs, with a :class:`RuntimeError`. With a noise model each circuit
    is first written in the model's gates, as its errors only touch a gate of their
    name; a fragment runs on qubits ``0..n-1``, so errors set on particular qubits
    apply to those numbers, not to where the fragment was in the circuit.
    """

    @property
    def uses_gpu(self) -> bool:
        return self._use_gpu

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("AER_BACKEND", config)
        self._use_gpu: bool = config.extra.get("use_gpu", False)
        self._noise_model = config.extra.get("noise_model", None)
        self._precision: str = config.extra.get("precision", "double")
        self._simulator_cache: AerSimulator | None = None
        self._estimator = AerEstimator()
        logger.info(
            "AerBackend initialised (gpu=%s, precision=%s, noise=%s)",
            self._use_gpu,
            self._precision,
            self._noise_model is not None,
        )

    @property
    def _simulator(self) -> AerSimulator:
        """The Aer simulator, built lazily (evaluate_expectation reuses the cutting
        code's own simulator, so a subcircuit-only backend never builds one here)."""
        if self._simulator_cache is None:
            self._simulator_cache = self._build_simulator(
                self._use_gpu, self._noise_model, self._precision
            )
        return self._simulator_cache

    @property
    def _plain(self) -> bool:
        """True with the default options: the cutting code's own CPU simulator gives
        the same result. Any option has to run on :attr:`_simulator`, which has it."""
        return (
            not self._use_gpu
            and self._noise_model is None
            and self._precision == "double"
        )

    def _counts(self, circuit: QuantumCircuit, shots: int, seed) -> dict:
        """``get_counts()`` of ``circuit`` with every qubit measured at the end, on
        the backend's own simulator (noise, device and precision applied)."""
        measured = circuit.copy()
        measured.measure_all()
        return self._run_for_counts(measured, shots, seed)

    def _run_for_counts(self, circuit: QuantumCircuit, shots: int, seed) -> dict:
        """``get_counts()`` of a measured ``circuit`` on the backend's simulator, in
        the noise model's gates if it has one: an error only touches a gate of its
        name, so gates under other names would run silently noiseless."""
        if self._noise_model is not None:
            circuit = qiskit_transpile(
                circuit, backend=self._simulator, optimization_level=0
            )
        job = self._simulator.run([circuit], shots=shots, seed_simulator=seed)
        return job.result().get_counts()

    def transpile(self, circuit: SupportedCircuit) -> QuantumCircuit:
        """Transpile the circuit for AerSimulator (accepts Qiskit or Qibo input)."""
        return qiskit_transpile(
            CircuitConverter.to_qiskit(circuit), backend=self._simulator
        )

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic", seed=None
    ):
        """Evaluate a gate-cut subcircuit via Qiskit Aer (CPU, or GPU when configured)."""
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _expec_value_qiskit,
            _expec_value_qiskit_gpu,
            _expectation_from_counts_qiskit,
        )

        subcirc = dag_to_circuit_qiskit(dag, num_qubits)
        if self._noise_model is not None or self._precision != "double":
            circuit, observables = subcirc
            return _expectation_from_counts_qiskit(
                self._counts(circuit, shots, seed),
                observables,
                shots,
                prep_qubits=prep_qubits,
            )
        if self._use_gpu:
            # _expec_value_qiskit_gpu recovers the cut (prep) qubits from the
            # subcircuit's own mid-circuit measurements and runs the same decoupled
            # reconstruction as CPU, so identity-on-a-cut-qubit works on GPU too --
            # the passed prep_qubits aren't needed here.
            return _expec_value_qiskit_gpu(subcirc, shots=shots, method=method)
        return _expec_value_qiskit(
            subcirc, shots=shots, method=method, prep_qubits=prep_qubits, seed=seed
        )

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Sample a subcircuit; return qubit-0-first counts."""
        from Qdislib.core.backend.base_backend import sample_with_qiskit

        return sample_with_qiskit(
            dag,
            num_qubits,
            lambda qc: self._run_for_counts(qc, shots, seed),
        )

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """Joint final/cut-point counts via Qiskit Aer (reuses the cut-measure
        register). With a noise model the cut points stay mid-circuit, as the noisy
        device would run them; on a noiseless GPU they are deferred, so the state is
        simulated once instead of once per shot."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _defer_cut_measurements,
            _deferred_joint_counts,
            _joint_counts_qiskit,
            _joint_from_qiskit_counts,
        )
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        cut_points = list(cut_point_qubits)
        if self._plain:
            return _joint_counts_qiskit(dag, num_qubits, cut_points, shots, seed=seed)
        if self._noise_model is None:
            deferred = _defer_cut_measurements(dag, num_qubits, cut_points)
            if deferred is not None:
                return _deferred_joint_counts(
                    self, deferred, num_qubits, len(cut_points), shots, seed
                )
        circuit = dag_to_circuit_qiskit(dag, num_qubits)[0]
        return _joint_from_qiskit_counts(
            self._counts(circuit, shots, seed), circuit, num_qubits, cut_points
        )

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        prepared = (
            self.transpile(circuit)
            if transpile
            else CircuitConverter.to_qiskit(circuit)
        )
        # AerSimulator.run() needs measurements present.
        if not prepared.count_ops().get("measure"):
            prepared = prepared.copy()
            prepared.measure_all()

        aer_job = self._simulator.run([prepared], shots=shots)
        counts = aer_job.result().get_counts()
        return SampleJobResult(
            counts=counts,
            backend_name=self.name,
            shots=shots,
            job_id=job.job_id,
            metadata={"aer_job_id": aer_job.job_id()},
        )

    def _execute_estimator(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        observable: SparsePauliOp | str,
        job: Job,
    ) -> ExpectationJobResult:
        prepared = (
            self.transpile(circuit)
            if transpile
            else CircuitConverter.to_qiskit(circuit)
        )
        if isinstance(observable, str):
            observable = SparsePauliOp(CircuitConverter.to_qiskit_pauli(observable))

        aer_job = self._estimator.run([(prepared, observable)])
        expectation_value = float(aer_job.result()[0].data.evs.item())
        return ExpectationJobResult(
            expectation_value=expectation_value,
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=str(observable),
            metadata={"aer_job_id": aer_job.job_id()},
        )

    # ------------------------------------------------------------------
    # Simulator construction
    # ------------------------------------------------------------------

    def _build_simulator(
        self, use_gpu: bool, noise_model: Any | None, precision: str
    ) -> AerSimulator:
        options: dict[str, Any] = {"precision": precision}
        if use_gpu:
            self._validate_gpu()
            options.update(
                {
                    "method": _GPU_METHOD,
                    "device": _GPU_DEVICE,
                    "cuStateVec_enable": _GPU_CUSTATEVEC_FLAG,
                }
            )
        else:
            options["method"] = _CPU_METHOD
        if noise_model is not None:
            options["noise_model"] = noise_model
        return AerSimulator(**options)

    @staticmethod
    def _validate_gpu() -> None:
        try:
            available = AerSimulator().available_devices()
        except Exception as exc:
            raise RuntimeError(
                "Could not query AerSimulator devices. Is qiskit-aer-gpu installed?"
            ) from exc
        if _GPU_DEVICE not in available:
            raise RuntimeError(
                f"GPU requested but '{_GPU_DEVICE}' is not available: {available}. "
                "Install qiskit-aer-gpu and ensure CUDA drivers are present."
            )
