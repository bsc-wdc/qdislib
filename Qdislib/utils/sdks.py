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

"""The SDKs Qdislib reads circuits from and writes subcircuits in, in one table.

An :class:`SDK` says how to recognise its circuits, how to turn one into a Qdislib
DAG and a DAG back into one, and which registered backend runs it. Everything that
dispatches on the SDK (``convert_to_dag``, the subcircuit builders, the
``backend=`` / ``software=`` names, ...) asks this table, so supporting an SDK is a
converter module and one :func:`register_sdk` call::

    register_sdk(SDK(
        name="mysdk",
        package="mysdk",
        circuit_types=("mysdk.circuit:Circuit",),
        to_dag="my_adapter:circuit_to_dag",
        from_dag="my_adapter:dag_to_circuit",
        backend="MYSDK_BACKEND",
        gates="my_adapter:gates",
        num_qubits="my_adapter:num_qubits",
    ))

Functions are given as ``"module:attribute"`` strings, imported on first use, so an
SDK that is not installed costs nothing.
"""

import dataclasses
import importlib
import importlib.util

__all__ = ["SDK", "register_sdk", "sdk_names", "sdk", "sdk_of"]


@dataclasses.dataclass(frozen=True)
class SDK:
    """One SDK: what its circuits are and how they convert."""

    #: The name users pass as ``backend=`` / ``software=``.
    name: str
    #: The importable package (to tell whether it is installed).
    package: str
    #: ``"module:Class"`` of every circuit type it reads.
    circuit_types: tuple
    #: ``"module:function"``: ``(circuit, obs_I=None) -> DAG``.
    to_dag: str
    #: ``"module:function"``: ``(dag, num_qubits) -> [circuit, observables]``.
    from_dag: str
    #: The registered backend (see ``Qdislib.core.backend.registry``) that runs it.
    backend: str
    #: ``"module:function"``: ``circuit -> [(label, qubits, is_measurement), ...]``
    #: in circuit order (``label`` as the DAG names the gate).
    gates: str = None
    #: ``"module:function"``: ``circuit -> number of qubits``.
    num_qubits: str = None
    #: ``"module:function"``: ``circuit -> exact state``, qubit 0 the most
    #: significant bit (optional; exact mode falls back to Qiskit or qibo).
    statevector: str = None

    def installed(self):
        return importlib.util.find_spec(self.package) is not None

    def call(self, field, *args, **kwargs):
        module, attribute = getattr(self, field).split(":")
        return getattr(importlib.import_module(module), attribute)(*args, **kwargs)

    def owns(self, circuit):
        """True if ``circuit`` is one of this SDK's circuits (without importing an
        SDK that ``circuit`` cannot come from)."""
        module_of = type(circuit).__module__.split(".")[0]
        for spec in self.circuit_types:
            module, attribute = spec.split(":")
            if module.split(".")[0] != module_of:
                continue
            try:
                cls = getattr(importlib.import_module(module), attribute)
            except (ImportError, AttributeError):
                continue
            if isinstance(circuit, cls):
                return True
        return False


_SDKS = {}


def register_sdk(entry):
    """Add (or replace) an SDK."""
    _SDKS[entry.name] = entry


def sdk_names():
    """The SDK names, in registration order."""
    return tuple(_SDKS)


def sdk(name):
    """The :class:`SDK` called ``name``, or ``None``."""
    return _SDKS.get(name)


def sdk_of(circuit):
    """The :class:`SDK` a circuit comes from, or ``None``."""
    return next((entry for entry in _SDKS.values() if entry.owns(circuit)), None)


_U = "Qdislib.utils"
for _entry in (
    SDK(
        name="qibo",
        package="qibo",
        circuit_types=("qibo.models.circuit:Circuit",),
        to_dag=f"{_U}.graph_qibo:circuit_qibo_to_dag",
        from_dag=f"{_U}.graph_qibo:dag_to_circuit_qibo",
        backend="QIBO_BACKEND",
        gates=f"{_U}.graph_qibo:qibo_gate_list",
        num_qubits=f"{_U}.graph_qibo:qibo_num_qubits",
        statevector=f"{_U}.graph_qibo:qibo_statevector",
    ),
    SDK(
        name="qiskit",
        package="qiskit",
        circuit_types=("qiskit.circuit:QuantumCircuit",),
        to_dag=f"{_U}.graph_qiskit:circuit_qiskit_to_dag",
        from_dag=f"{_U}.graph_qiskit:dag_to_circuit_qiskit",
        backend="AER_BACKEND",
        gates=f"{_U}.graph_qiskit:qiskit_gates",
        num_qubits=f"{_U}.graph_qiskit:qiskit_num_qubits",
        statevector=f"{_U}.graph_qiskit:qiskit_statevector",
    ),
    SDK(
        name="cudaq",
        package="cudaq",
        circuit_types=("cudaq.kernel.kernel_decorator:PyKernelDecorator",),
        to_dag=f"{_U}.graph_cudaq:circuit_cudaq_to_dag",
        from_dag=f"{_U}.graph_cudaq:dag_to_circuit_cudaq_kernel",
        backend="CUDAQ_BACKEND",
        gates=f"{_U}.graph_cudaq:cudaq_gates",
        num_qubits=f"{_U}.graph_cudaq:cudaq_num_qubits",
        statevector=f"{_U}.graph_cudaq:cudaq_statevector",
    ),
    SDK(
        name="pennylane",
        package="pennylane",
        circuit_types=("pennylane.tape:QuantumScript", "pennylane:QNode"),
        to_dag=f"{_U}.graph_pennylane:circuit_pennylane_to_dag",
        from_dag=f"{_U}.graph_pennylane:dag_to_circuit_pennylane",
        backend="PENNYLANE_BACKEND",
        gates=f"{_U}.graph_pennylane:pennylane_gates",
        num_qubits=f"{_U}.graph_pennylane:pennylane_num_qubits",
        statevector=f"{_U}.graph_pennylane:pennylane_statevector",
    ),
    SDK(
        name="qilisdk",
        package="qilisdk",
        circuit_types=("qilisdk.digital.circuit:Circuit",),
        to_dag=f"{_U}.graph_qilisdk:circuit_qilisdk_to_dag",
        from_dag=f"{_U}.graph_qilisdk:dag_to_circuit_qilisdk",
        backend="QILISDK_BACKEND",
        gates=f"{_U}.graph_qilisdk:qilisdk_gates",
        num_qubits=f"{_U}.graph_qilisdk:qilisdk_num_qubits",
        statevector=f"{_U}.graph_qilisdk:qilisdk_statevector",
    ),
):
    register_sdk(_entry)
