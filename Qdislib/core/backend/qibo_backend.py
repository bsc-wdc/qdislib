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

"""Qibo simulator backend."""

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
from Qdislib.core.backend.circuit_converter import CircuitConverter

logger = logging.getLogger(__name__)


class QiboBackend(BaseBackend):
    """Backend wrapping the Qibo simulator.

    Configuration keys (``BackendConfig.extra``)
    ---------------------------------------------
    platform : str, optional
        Qibo backend/platform name (default ``"numpy"``).
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("QIBO_BACKEND", config)
        self._platform: str = config.extra.get("platform", "numpy")
        logger.info("QiboBackend initialised (platform=%s)", self._platform)

    def transpile(self, circuit: SupportedCircuit) -> Any:
        """Qibo runs the circuit directly; just normalise to a Qibo circuit."""
        return CircuitConverter.to_qibo(circuit)

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        import qibo
        from qibo import gates

        qibo.set_backend(self._platform)
        qc = CircuitConverter.to_qibo(circuit)
        if not qc.measurements:
            qc.add(gates.M(*range(qc.nqubits)))
        result = qc(nshots=shots)
        counts = {k: int(v) for k, v in result.frequencies(binary=True).items()}
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
        """Shot-based Pauli expectation. ``observable`` is a per-qubit string over I/X/Y/Z."""
        import qibo
        from qibo import gates

        qibo.set_backend(self._platform)
        qc = CircuitConverter.to_qibo(circuit).copy(deep=True)

        # Rotate X/Y into the Z basis so a computational-basis measurement yields <O>.
        for qubit, pauli in enumerate(observable):
            if pauli == "X":
                qc.add(gates.H(qubit))
            elif pauli == "Y":
                qc.add(gates.SDG(qubit))
                qc.add(gates.H(qubit))
        qc.add(gates.M(*range(qc.nqubits)))

        frequencies = qc(nshots=shots).frequencies(binary=True)
        measured = [q for q, pauli in enumerate(observable) if pauli != "I"]
        expectation = 0.0
        for bitstring, count in frequencies.items():
            parity = sum(int(bitstring[q]) for q in measured) % 2
            expectation += ((-1) ** parity) * count
        expectation /= shots
        return ExpectationJobResult(
            expectation_value=float(expectation),
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=observable,
        )
