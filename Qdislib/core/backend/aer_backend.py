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
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("AER_BACKEND", config)
        self._use_gpu: bool = config.extra.get("use_gpu", False)
        self._noise_model = config.extra.get("noise_model", None)
        self._precision: str = config.extra.get("precision", "double")
        self._simulator_cache: AerSimulator | None = None
        self._estimator = AerEstimator()
        logger.info(
            "AerBackend initialised (gpu=%s, precision=%s, noise=%s)",
            self._use_gpu, self._precision, self._noise_model is not None,
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

    def transpile(self, circuit: SupportedCircuit) -> QuantumCircuit:
        """Transpile the circuit for AerSimulator (accepts Qiskit or Qibo input)."""
        return qiskit_transpile(CircuitConverter.to_qiskit(circuit), backend=self._simulator)

    def evaluate_expectation(self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic"):
        """Evaluate a gate-cut subcircuit via Qiskit Aer (CPU, or GPU when configured)."""
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _expec_value_qiskit,
            _expec_value_qiskit_gpu,
        )

        subcirc = dag_to_circuit_qiskit(dag, num_qubits)
        if self._use_gpu:
            # The GPU path does not implement the cut-point sign recovery; the
            # cut guard forbids identity-on-a-cut-qubit for GPU, so prep_qubits
            # is empty here.
            return _expec_value_qiskit_gpu(subcirc, shots=shots, method=method)
        return _expec_value_qiskit(
            subcirc, shots=shots, method=method, prep_qubits=prep_qubits
        )

    def evaluate_samples(self, dag, num_qubits, shots=1024):
        """Sample a subcircuit; return qubit-0-first counts (reversed from Qiskit's little-endian)."""
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        built = dag_to_circuit_qiskit(dag, num_qubits)
        qc = built[0] if isinstance(built, (list, tuple)) else built
        qc = qc.copy()
        qc.measure_all()
        counts = self._simulator.run([qc], shots=shots).result().get_counts()
        out: dict[str, int] = {}
        for key, value in counts.items():
            # keep the measure_all register (first) and reverse to qubit-0-first
            bitstring = key.split()[0][::-1]
            out[bitstring] = out.get(bitstring, 0) + value
        return out

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        prepared = self.transpile(circuit) if transpile else CircuitConverter.to_qiskit(circuit)
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
        prepared = self.transpile(circuit) if transpile else CircuitConverter.to_qiskit(circuit)
        if isinstance(observable, str):
            observable = SparsePauliOp(observable)

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

    def _build_simulator(self, use_gpu: bool, noise_model: Any | None, precision: str) -> AerSimulator:
        options: dict[str, Any] = {"precision": precision}
        if use_gpu:
            self._validate_gpu()
            options.update(
                {"method": _GPU_METHOD, "device": _GPU_DEVICE, "cuStateVec_enable": _GPU_CUSTATEVEC_FLAG}
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
            raise RuntimeError("Could not query AerSimulator devices. Is qiskit-aer-gpu installed?") from exc
        if _GPU_DEVICE not in available:
            raise RuntimeError(
                f"GPU requested but '{_GPU_DEVICE}' is not available: {available}. "
                "Install qiskit-aer-gpu and ensure CUDA drivers are present."
            )
