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

"""IQM quantum-hardware backend."""

from typing import Any

from qiskit.quantum_info import SparsePauliOp

from Qdislib.core.backend.base_backend import (
    BaseBackend,
    BackendConfig,
    Job,
    SupportedCircuit,
    SampleJobResult,
    ExpectationJobResult,
)
from Qdislib.core.backend.circuit_converter import CircuitConverter


class IQMBackend(BaseBackend):
    """Backend for IQM quantum hardware.

    Configuration keys (``BackendConfig.extra``)
    ---------------------------------------------
    url : str
        IQM server URL.
    tokens_path : str
        Path to the IQM authentication tokens file.
    project : str, optional
        Project name, attached to transpiled circuit metadata.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("IQM_BACKEND", config)
        try:
            import os
            from iqm.qiskit_iqm import IQMProvider
        except ImportError as exc:
            raise ImportError(
                "IQMBackend requires the IQM extra. Install with:  pip install qdislib[iqm]"
            ) from exc

        os.environ["IQM_TOKENS_FILE"] = config.extra["tokens_path"]
        provider = IQMProvider(config.extra["url"], quantum_computer="default")
        self.backend = provider.get_backend()

    def transpile(self, circuit: SupportedCircuit) -> Any:
        from qiskit import transpile as qiskit_transpile

        transpiled = qiskit_transpile(
            CircuitConverter.to_qiskit(circuit), backend=self.backend
        )
        transpiled.metadata["project"] = self.config.extra.get("project", "")
        return transpiled

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        prepared = self.transpile(circuit) if transpile else circuit
        counts = self.backend.run(prepared, shots=shots).result().get_counts()
        return SampleJobResult(
            counts=counts, backend_name=self.name, shots=shots, job_id=job.job_id
        )

    def _execute_estimator(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        observable: SparsePauliOp | str,
        job: Job,
    ) -> ExpectationJobResult:
        prepared = self.transpile(circuit) if transpile else circuit
        counts = self.backend.run(prepared, shots=shots).result().get_counts()

        if isinstance(observable, str):
            paulis, coeffs = [CircuitConverter.to_qiskit_pauli(observable)], [1.0]
        elif isinstance(observable, SparsePauliOp):
            simplified = observable.simplify()
            paulis = list(simplified.paulis)
            coeffs = list(simplified.coeffs)
        else:
            raise TypeError(
                f"observable must be SparsePauliOp or str, got {type(observable)}"
            )

        expectation_value = 0.0
        for pauli, coeff in zip(paulis, coeffs):
            pauli_str = pauli if isinstance(pauli, str) else pauli.to_label()
            if "X" in pauli_str or "Y" in pauli_str:
                raise ValueError(
                    f"Pauli term '{pauli_str}' contains X or Y, which need basis-change "
                    "gates not yet supported by IQMBackend."
                )
            term_ev = 0.0
            for bitstring, count in counts.items():
                eigenvalue = 1.0
                for pauli_char, bit in zip(pauli_str, bitstring):
                    if pauli_char == "Z":
                        eigenvalue *= 1.0 if bit == "0" else -1.0
                term_ev += (count / shots) * eigenvalue
            expectation_value += coeff * term_ev

        return ExpectationJobResult(
            expectation_value=float(expectation_value.real),
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=str(observable),
        )
