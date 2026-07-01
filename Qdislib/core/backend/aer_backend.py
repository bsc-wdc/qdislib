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
        use_gpu: bool = config.extra.get("use_gpu", False)
        noise_model = config.extra.get("noise_model", None)
        precision: str = config.extra.get("precision", "double")

        self._simulator = self._build_simulator(use_gpu, noise_model, precision)
        self._estimator = AerEstimator()
        logger.info(
            "AerBackend initialised (gpu=%s, precision=%s, noise=%s)",
            use_gpu, precision, noise_model is not None,
        )

    def transpile(self, circuit: SupportedCircuit) -> QuantumCircuit:
        """Transpile the circuit for AerSimulator (accepts Qiskit or Qibo input)."""
        return qiskit_transpile(CircuitConverter.to_qiskit(circuit), backend=self._simulator)

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
