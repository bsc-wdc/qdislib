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

# -*- coding: utf-8 -*-

"""QiliSDK backend: subcircuits on Qilimanjaro's QiliSim simulator.

The cutting code hands it a Qdislib DAG, built into a QiliSDK ``Circuit``
(:mod:`Qdislib.utils.graph_qilisdk`) and sampled on QiliSim. QiliSim samples the
final state, so a measure-and-prepare gate cut's mid-circuit measurement is deferred
(a CNOT onto a fresh qubit, read at the end; the same channel).

Configuration keys (``BackendConfig.extra``)
---------------------------------------------
gpu : bool, optional
    Run QiliSim on the GPU (default ``False``).
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
    marginalize_identity_qubits,
)

logger = logging.getLogger(__name__)


def _renumbered(dag):
    """``(dag, width, identity)``: the DAG on qubits ``0..width-1`` without its
    ``Observable I`` markers, and the qubits those markers named."""
    dag = dag.copy()
    identity = set()
    for node in list(dag.nodes):
        if (dag.nodes[node].get("gate") or "").upper() == "OBSERVABLE I":
            identity.add(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)
    used = sorted({q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()})
    index = {q: i for i, q in enumerate(used)}
    for node in dag.nodes:
        dag.nodes[node]["qubits"] = tuple(index[q] for q in dag.nodes[node]["qubits"])
    return dag, len(used), {index[q] for q in identity if q in index}


class QiliSDKBackend(BaseBackend):
    """Backend running Qdislib subcircuits on QiliSDK's QiliSim."""

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("QILISDK_BACKEND", config)
        import importlib.util

        if importlib.util.find_spec("qilisdk") is None:
            raise ImportError(
                "QiliSDKBackend requires QiliSDK, on Python 3.11 or later. Install with:  pip install qilisdk"
            )
        self._gpu: bool = bool(config.extra.get("gpu", False))

    def transpile(self, circuit: SupportedCircuit) -> Any:
        return _qilisdk_circuit(circuit)

    # ------------------------------------------------------------------
    # Cutting hooks -- a Qdislib subcircuit given as a DAG
    # ------------------------------------------------------------------

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        """Sample a subcircuit; qubit-0-first counts, identity-observed qubits
        marginalised out."""
        from Qdislib.utils.graph_qilisdk import dag_to_circuit_qilisdk, sample_qilisdk

        circuit, observables = dag_to_circuit_qilisdk(dag, num_qubits)
        counts = sample_qilisdk(circuit, shots, seed, gpu=self._gpu)
        return marginalize_identity_qubits(counts, observables)

    def evaluate_joint_counts(
        self, dag, num_qubits, cut_point_qubits=(), shots=1024, seed=None
    ):
        """``{(final_bits, cut_point_bits): count}``, the cut-point measurements
        deferred to the end."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _defer_cut_measurements,
            _deferred_joint_counts,
        )

        prep = sorted(cut_point_qubits)
        if not prep:
            counts = self.evaluate_samples(dag, num_qubits, shots, seed)
            return {(tuple(int(b) for b in k), ()): c for k, c in counts.items()}
        deferred = _defer_cut_measurements(dag, num_qubits, prep)
        if deferred is None:
            raise NotImplementedError(
                "QiliSDK: a qubit measured twice mid-circuit cannot be deferred"
            )
        return _deferred_joint_counts(
            self, deferred, num_qubits, len(prep), shots, seed
        )

    def evaluate_expectation(
        self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic", seed=None
    ):
        """``<O>`` of a cut subcircuit, ``O`` the Z/I string its ``Observable I``
        markers give, with the measure-and-prepare sign rule for cut qubits."""
        from Qdislib.core.cutting_algorithms.gate_cutting import _z_correlators

        dag, width, identity = _renumbered(dag)
        prep = sorted(
            d["qubits"][0]
            for _, d in dag.nodes(data=True)
            if (d.get("gate") or "").lower() == "measure"
        )
        joint = self.evaluate_joint_counts(dag, width, prep, shots, seed)
        mask = sum(1 << q for q in range(width) if q not in identity)
        return _z_correlators(joint, prep, width, [mask])[mask]

    # ------------------------------------------------------------------
    # Generic queue hooks (standalone backend use)
    # ------------------------------------------------------------------

    def _execute_samples(
        self, circuit: SupportedCircuit, transpile: bool, shots: int, job: Job
    ) -> SampleJobResult:
        from Qdislib.utils.graph_qilisdk import sample_qilisdk

        counts = sample_qilisdk(_qilisdk_circuit(circuit), shots, gpu=self._gpu)
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
        from Qdislib.utils.graph_qilisdk import sample_qilisdk

        if set(observable) - {"Z", "I"}:
            raise NotImplementedError("QiliSDKBackend estimates Z/I observables")
        counts = sample_qilisdk(_qilisdk_circuit(circuit), shots, gpu=self._gpu)
        total = sum(counts.values()) or 1
        value = 0.0
        for bits, count in counts.items():
            parity = sum(int(b) for b, p in zip(bits, observable) if p == "Z") % 2
            value += (-1) ** parity * count
        return ExpectationJobResult(
            expectation_value=value / total,
            shots=shots,
            backend_name=self.name,
            job_id=job.job_id,
            obs=observable,
        )


def _qilisdk_circuit(circuit):
    """A QiliSDK ``Circuit`` from a circuit, a ``[circuit, observables]`` pair or a
    ``(dag, num_qubits)`` payload."""
    if isinstance(circuit, (list, tuple)) and len(circuit) == 2:
        first = circuit[0]
        if hasattr(first, "nodes"):
            from Qdislib.utils.graph_qilisdk import dag_to_circuit_qilisdk

            return dag_to_circuit_qilisdk(first, circuit[1])[0]
        return first
    return circuit
