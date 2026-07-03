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

"""PennyLane backend.

Executes subcircuits with PennyLane's own engine -- a :class:`~pennylane.QNode`
bound to a :func:`~pennylane.device`. Defaults to the reference ``default.qubit``
device (pure Python, always available, full mid-circuit-measurement support);
set ``extra={"device": "lightning.qubit"}`` to use the compiled Lightning device.

The cutting code feeds this backend via Qdislib's backend-neutral DAG (built into
a :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit`), exactly as it does
for CUDA-Q; the actual sampling / expectation always runs through a PennyLane
QNode built on the device.
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


def _ops_and_width(circuit: SupportedCircuit):
    """Return ``(ops, num_qubits)`` for a PennyLane payload.

    Accepts a :class:`~Qdislib.utils.graph_pennylane.PennylaneCircuit`, a
    ``[PennylaneCircuit, observables]`` pair, or a ``(dag, num_qubits)`` payload
    (which it materialises).
    """
    from Qdislib.utils.graph_pennylane import PennylaneCircuit, dag_to_circuit_pennylane

    payload = circuit
    if (
        isinstance(circuit, (list, tuple))
        and len(circuit) == 2
        and isinstance(circuit[1], str)
    ):
        payload = circuit[0]
    if isinstance(payload, PennylaneCircuit):
        return payload.ops, payload.num_qubits
    if (
        isinstance(payload, (list, tuple))
        and len(payload) == 2
        and hasattr(payload[0], "nodes")
    ):
        built = dag_to_circuit_pennylane(payload[0], payload[1])
        return built[0].ops, built[0].num_qubits
    raise TypeError(
        f"PennylaneBackend cannot execute {type(circuit).__name__}; "
        "expected a PennylaneCircuit or (dag, num_qubits) payload."
    )


class PennylaneBackend(BaseBackend):
    """Backend wrapping PennyLane sampling / expectation on a QNode + device.

    Configuration keys (``BackendConfig.extra``)
    ---------------------------------------------
    device : str, optional
        PennyLane device short name (default ``"default.qubit"``); e.g.
        ``"lightning.qubit"`` for the compiled state-vector simulator.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("PENNYLANE_BACKEND", config)
        import importlib.util

        if importlib.util.find_spec("pennylane") is None:
            raise ImportError(
                "PennylaneBackend requires PennyLane. "
                "Install with:  pip install qdislib[pennylane]"
            )
        self._device: str = config.extra.get("device", "default.qubit")

    def transpile(self, circuit: SupportedCircuit) -> Any:
        """Materialise a PennyLane subcircuit (a ``PennylaneCircuit``)."""
        from Qdislib.utils.graph_pennylane import PennylaneCircuit

        ops, num_qubits = _ops_and_width(circuit)
        return PennylaneCircuit(ops, num_qubits)

    # ------------------------------------------------------------------
    # Cutting hooks -- evaluate a Qdislib subcircuit given as a DAG
    # ------------------------------------------------------------------

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic"
    ):
        """Evaluate a cut subcircuit via a PennyLane QNode."""
        from Qdislib.utils.graph_pennylane import dag_to_circuit_pennylane
        from Qdislib.core.cutting_algorithms.wire_cutting import _expec_value_pennylane

        subcirc = dag_to_circuit_pennylane(dag, num_qubits)
        return _expec_value_pennylane(
            subcirc, shots, prep_qubits=list(prep_qubits), device=self._device
        )

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Sample a subcircuit via PennyLane; return qubit-0-first counts."""
        from Qdislib.utils.graph_pennylane import dag_to_circuit_pennylane
        from Qdislib.core.cutting_algorithms.wire_cutting import _counts_pennylane

        built = dag_to_circuit_pennylane(dag, num_qubits)
        circuit = built[0] if isinstance(built, (list, tuple)) else built
        return _counts_pennylane(circuit, shots, device=self._device, seed=seed)

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """Joint final/cut-point counts via PennyLane mid-circuit measurements."""
        from Qdislib.core.cutting_algorithms.gate_cutting import _joint_counts_pennylane

        return _joint_counts_pennylane(
            dag, num_qubits, list(cut_point_qubits), shots, seed=seed, device=self._device
        )

    # ------------------------------------------------------------------
    # Generic FIFO-queue hooks (standalone backend use, not the cut path)
    # ------------------------------------------------------------------

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        from Qdislib.core.cutting_algorithms.wire_cutting import _counts_pennylane
        from Qdislib.utils.graph_pennylane import PennylaneCircuit

        ops, num_qubits = _ops_and_width(circuit)
        counts = _counts_pennylane(
            PennylaneCircuit(ops, num_qubits), shots, device=self._device
        )
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
        import pennylane as qml
        from Qdislib.utils.graph_pennylane import _apply_pennylane_gate

        ops, num_qubits = _ops_and_width(circuit)

        # Build the tensor-product Pauli observable from the per-qubit string.
        factory = {"X": qml.PauliX, "Y": qml.PauliY, "Z": qml.PauliZ}
        obs_op = None
        for qubit, pauli in enumerate(observable):
            if pauli == "I":
                continue
            term = factory[pauli](qubit)
            obs_op = term if obs_op is None else obs_op @ term
        if obs_op is None:  # all-identity observable
            return ExpectationJobResult(
                1.0, shots, self.name, job.job_id, obs=observable
            )

        dev = qml.device(self._device, wires=num_qubits)

        def qfunc():
            for gate, qubits, params in ops:
                _apply_pennylane_gate(gate, qubits, params)
            return qml.expval(obs_op)

        qnode = qml.set_shots(qml.QNode(qfunc, dev), shots=shots)
        expectation = float(qnode())
        return ExpectationJobResult(
            expectation_value=expectation,
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=observable,
        )
