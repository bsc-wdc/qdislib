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
Only :meth:`evaluate_expectation` is used by the cutting code; the generic
sample/estimator hooks are left unimplemented.
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

    def evaluate_expectation(self, dag, num_qubits, prep_qubits=(), shots=1024, method="numpy"):
        from Qdislib.utils.graph_qibo import dag_to_circuit_qibo
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qibo_qpu

        subcirc = dag_to_circuit_qibo(dag, num_qubits)
        return _expec_value_qibo_qpu(subcirc, shots=shots, method=method)


class IBMBackend(_QpuBackendBase):
    """IBM Quantum hardware via qiskit-ibm-runtime.

    Delegates to :func:`~Qdislib.core.cutting_algorithms.gate_cutting._expec_value_qiskit_qpu`.
    ``BackendConfig.extra`` carries the runtime ``batch`` (Session/Batch) and
    ``backend`` (the IBM backend) objects.
    """

    def __init__(self, config: BackendConfig) -> None:
        super().__init__("IBM_BACKEND", config)

    def evaluate_expectation(self, dag, num_qubits, prep_qubits=(), shots=1024, method="automatic"):
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit_qpu

        subcirc = dag_to_circuit_qiskit(dag, num_qubits)
        return _expec_value_qiskit_qpu(
            subcirc,
            shots=shots,
            batch=self.config.extra.get("batch"),
            backend=self.config.extra.get("backend"),
        )
