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

import logging
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

logger = logging.getLogger(__name__)


class IQMBackend(BaseBackend):
    """Backend for IQM quantum hardware.

    Configuration keys (named directly on :class:`BackendConfig`)
    -------------------------------------------------------------
    url : str
        Base URL of the IQM server, e.g. ``"https://resonance.iqm.tech"``. The
        quantum computer's name does *not* belong here -- pass it as
        ``quantum_computer``. Falls back to the ``IQM_SERVER_URL`` env var.
    token : str, optional
        Long-lived authentication token (IQM Resonance / cloud). Falls back to
        the ``IQM_TOKEN`` env var.
    tokens_path : str, optional
        Path to a tokens *file* (the on-premises Cortex-CLI flow), used only when
        no ``token`` is given. Falls back to the ``IQM_TOKENS_FILE`` env var.
    quantum_computer : str, optional
        Alias of the quantum computer to run on (e.g. ``"garnet"``) when the
        server hosts more than one. ``None`` (default) means the server's default.
    project : str, optional
        Free-form label attached to each transpiled circuit's metadata and echoed
        back with the job. The IQM client does not route or bill on it.

    ``BackendConfig.num_qubits`` can be omitted: the device is asked for its own
    width on connection. State it only to *restrict* a device to smaller fragments.
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

        # Auth is passed explicitly rather than through the environment: a plain
        # cloud token and a Cortex-CLI tokens file are mutually exclusive in
        # IQMClient, and anything omitted here falls back to IQM_TOKEN /
        # IQM_TOKENS_FILE / IQM_SERVER_URL, which the client reads itself.
        auth: dict[str, str] = {}
        if config.extra.get("token"):
            auth["token"] = config.extra["token"]
        elif config.extra.get("tokens_path"):
            auth["tokens_file"] = config.extra["tokens_path"]

        provider = IQMProvider(
            config.extra.get("url") or os.environ.get("IQM_SERVER_URL"),
            quantum_computer=config.extra.get("quantum_computer"),
            **auth,
        )
        self.backend = provider.get_backend()

        # The device knows its own width, so a caller who didn't state one gets
        # the real number rather than a guess that has to be kept in sync. An
        # explicit num_qubits still wins -- that's how you restrict a large
        # device to smaller fragments when QdBackEnd routes across backends.
        if not config.num_qubits:
            config.num_qubits = self.backend.num_qubits

    def transpile(self, circuit: SupportedCircuit) -> Any:
        from qiskit import transpile as qiskit_transpile

        transpiled = qiskit_transpile(
            CircuitConverter.to_qiskit(circuit), backend=self.backend
        )
        transpiled.metadata["project"] = self.config.extra.get("project", "")
        return transpiled

    def _counts(
        self, circuit: SupportedCircuit, transpile: bool, shots: int
    ) -> dict[str, int]:
        """Submit one circuit to the device and return its raw Qiskit counts."""
        prepared = self.transpile(circuit) if transpile else circuit
        job = self.backend.run(prepared, shots=shots)
        logger.info("IQM job submitted: %s", getattr(job, "job_id", lambda: "?")())
        return job.result().get_counts()

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        counts = self._counts(circuit, transpile, shots)
        return SampleJobResult(
            counts=counts, backend_name=self.name, shots=shots, job_id=job.job_id
        )

    # ------------------------------------------------------------------
    # Circuit-cutting hooks -- what lets a cut subcircuit run on the device
    # ------------------------------------------------------------------

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic", seed=None
    ):
        """Evaluate a gate-cut subcircuit DAG on the IQM device (all-in-one path)."""
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        return self.evaluate_subcircuit(
            dag_to_circuit_qiskit(dag, num_qubits), shots=shots, prep_qubits=prep_qubits
        )

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """Joint final/cut-point counts on the IQM device (the cut-point
        measurements run mid-circuit). ``seed`` is ignored (hardware)."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _joint_from_qiskit_counts,
        )
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        built = dag_to_circuit_qiskit(dag, num_qubits)
        if built is None:
            return None
        base = built[0] if isinstance(built, (list, tuple)) else built
        circuit = base.copy()
        circuit.measure_all()
        counts = self._counts(circuit, True, shots)
        return _joint_from_qiskit_counts(
            counts, base, num_qubits, sorted(cut_point_qubits)
        )

    def evaluate_subcircuit(self, subcircuit, shots: int = 1024, prep_qubits=None):
        """Evaluate one already-built cut subcircuit on the device (3-step workflow).

        Takes what ``gate_cutting_subcircuits`` / ``wire_cutting_subcircuits``
        return -- a ``[circuit, observables]`` pair, a bare circuit, or a
        serializable :class:`~Qdislib.utils.backends.Subcircuit` -- appends the
        final measurement, submits one job, and reconstructs the expectation with
        the same Z/I + cut-point-sign rule the simulators use
        (:func:`~Qdislib.core.cutting_algorithms.gate_cutting._expectation_from_counts_qiskit`),
        so a cut evaluated here recombines exactly as it does on Aer.

        ``prep_qubits`` (the cut qubits measured for the decomposition) is
        recovered from the subcircuit's own mid-circuit measurements when not given.
        """
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _expectation_from_counts_qiskit,
        )

        if subcircuit is None:  # empty component of a cut: contributes nothing
            return None
        if hasattr(subcircuit, "to_circuit"):  # serializable Subcircuit (PyCOMPSs)
            subcircuit = subcircuit.to_circuit()

        if isinstance(subcircuit, (list, tuple)):
            circuit, observables = subcircuit[0], subcircuit[1]
        else:
            circuit, observables = subcircuit, None

        circuit = CircuitConverter.to_qiskit(circuit).copy()
        if observables is None:
            observables = "Z" * circuit.num_qubits
        if prep_qubits is None:
            # Same recovery as the CPU evaluators: a mid-circuit measurement in a
            # cut subcircuit *is* the cut point, and its outcome carries the
            # decomposition sign.
            prep_qubits = [
                circuit.find_bit(instr.qubits[0]).index
                for instr in circuit.data
                if instr.operation.name == "measure"
            ]

        if any(reg.name == "meas" for reg in circuit.cregs):
            # The simulator evaluators append measure_all *in place*, so a
            # subcircuit that already carries a "meas" register has been
            # evaluated before. Submitting it would read the wrong register and
            # return a plausible-but-wrong number -- on hardware, after paying
            # for the shots. Cut a fresh set of subcircuits per evaluation.
            raise ValueError(
                "This subcircuit has already been evaluated (it carries a 'meas' "
                "register). Cut the circuit again before evaluating on hardware."
            )

        circuit.measure_all()
        counts = self._counts(circuit, True, shots)
        return _expectation_from_counts_qiskit(
            counts, observables, shots, prep_qubits=prep_qubits
        )

    def _execute_estimator(
        self,
        circuit: SupportedCircuit,
        transpile: bool,
        shots: int,
        observable: SparsePauliOp | str,
        job: Job,
    ) -> ExpectationJobResult:
        counts = self._counts(circuit, transpile, shots)

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
