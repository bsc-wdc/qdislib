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

"""CUDA-Q backend (GPU when available, CPU ``qpp-cpu`` otherwise).

Accepts a built CUDA-Q kernel, or the ``(dag, num_qubits)`` payload the cutting
code carries for CUDA-Q (which it builds into a kernel on demand). Qiskit/Qibo
circuits are not auto-converted -- CUDA-Q has no stable QASM round-trip, so the
cutting code feeds it via the DAG.
"""

import logging
from typing import Any

from Qdislib.core.backend.base_backend import (
    BaseBackend,
    BackendConfig,
    Job,
    SupportedCircuit,
    SampleJobResult,
    ExpectationJobResult,
)

logger = logging.getLogger(__name__)


def select_cudaq_target(target=None):
    """Set CUDA-Q's target: ``target`` when given (an error if it cannot be set),
    else the ``nvidia`` GPU, or ``qpp-cpu`` with a warning saying why not."""
    import cudaq

    if target is not None:
        cudaq.set_target(target)
        return
    try:
        cudaq.set_target("nvidia")
        logger.debug("CUDA-Q using GPU target 'nvidia'")
    except Exception as exc:
        from Qdislib.utils.logging_config import warn_once

        cudaq.set_target("qpp-cpu")
        warn_once(
            f"CUDA-Q cannot use its GPU target 'nvidia' ({exc}): running on the CPU "
            "target 'qpp-cpu'. To choose it yourself, pass target='qpp-cpu' to the "
            "backend."
        )


class CudaqBackend(BaseBackend):
    """Backend wrapping CUDA-Q sampling / observation.

    Configuration keys (``BackendConfig.extra``)
    ---------------------------------------------
    target : str, optional
        CUDA-Q target (default: the ``"nvidia"`` GPU, or ``"qpp-cpu"`` with a warning
        when it cannot be used). Under PyCOMPSs, name ``"nvidia"`` here to run in tasks
        that ask for a GPU: COMPSs assigns a GPU only to those.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("CUDAQ_BACKEND", config)
        import importlib.util

        if importlib.util.find_spec("cudaq") is None:
            raise ImportError(
                "CudaqBackend requires CUDA-Q. Install with:  pip install qdislib[cudaq]"
            )
        self._target: str | None = config.extra.get("target")

    @property
    def uses_gpu(self) -> bool:
        return str(self._target or "").startswith("nvidia")

    def _set_target(self) -> None:
        select_cudaq_target(self._target)

    def _to_kernel(self, circuit: SupportedCircuit):
        """Return a built CUDA-Q kernel from a kernel or a ``(dag, num_qubits)`` payload."""
        if (
            isinstance(circuit, (list, tuple))
            and len(circuit) == 2
            and hasattr(circuit[0], "nodes")
        ):
            from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel

            kernel, _obs = dag_to_circuit_cudaq_kernel(circuit[0], circuit[1])
            return kernel
        return circuit

    def transpile(self, circuit: SupportedCircuit) -> Any:
        return self._to_kernel(circuit)

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic", seed=None
    ):
        """Evaluate a gate-cut subcircuit via CUDA-Q."""
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq
        from Qdislib.core.cutting_algorithms.wire_cutting import _expec_value_cudaq

        subcirc = dag_to_circuit_cudaq(dag, num_qubits)
        if seed is not None:
            import cudaq

            cudaq.set_random_seed(seed)
        return _expec_value_cudaq(
            subcirc, shots, prep_qubits=prep_qubits, target=self._target
        )

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Sample a subcircuit via CUDA-Q; return qubit-0-first counts."""
        import cudaq
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel
        from Qdislib.core.backend.base_backend import marginalize_identity_qubits

        self._set_target()
        if seed is not None:
            cudaq.set_random_seed(seed)
        built = dag_to_circuit_cudaq_kernel(dag, num_qubits)
        kernel, observables = (
            built[:2] if isinstance(built, (list, tuple)) else (built, None)
        )
        result = cudaq.sample(kernel, shots_count=shots)
        counts = {str(k): int(v) for k, v in result.items()}
        return marginalize_identity_qubits(counts, observables)

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """Joint final/cut-point counts via CUDA-Q (cudaq.run collapses mid-circuit mz)."""
        from Qdislib.core.cutting_algorithms.gate_cutting import _joint_counts_cudaq

        return _joint_counts_cudaq(
            dag,
            num_qubits,
            list(cut_point_qubits),
            shots,
            seed=seed,
            target=self._target,
        )

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        import cudaq

        self._set_target()
        kernel = self._to_kernel(circuit)
        result = cudaq.sample(kernel, shots_count=shots)
        counts = {str(k): int(v) for k, v in result.items()}
        return SampleJobResult(
            counts=counts, backend_name=self.name, shots=shots, job_id=job.job_id
        )

    def _execute_estimator(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        observable: str,
        job: Job,
    ) -> ExpectationJobResult:
        import cudaq

        self._set_target()
        kernel = self._to_kernel(circuit)

        # Build a CUDA-Q spin operator from the per-qubit Pauli string.
        pauli = {"X": cudaq.spin.x, "Y": cudaq.spin.y, "Z": cudaq.spin.z}
        term = None
        for qubit, op in enumerate(observable):
            if op == "I":
                continue
            factor = pauli[op](qubit)
            term = factor if term is None else term * factor
        if term is None:  # all-identity observable
            return ExpectationJobResult(
                1.0, shots, self.name, job.job_id, obs=observable
            )

        expectation = cudaq.observe(kernel, term, shots_count=shots).expectation()
        return ExpectationJobResult(
            expectation_value=float(expectation),
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=observable,
        )
