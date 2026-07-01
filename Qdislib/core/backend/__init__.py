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

"""Pluggable, user-swappable execution backends for Qdislib.

Public surface::

    from Qdislib import QdBackEnd, register_backend

    be = QdBackEnd({"AER_BACKEND": {"num_qubits": 30}})
    result = be.sample_circuit(circuit, shots=4096)      # -> SampleJobResult
    result = be.get_expectation_value(circuit, "ZZI")    # -> ExpectationJobResult

Add your own backend by subclassing :class:`BaseBackend` and registering it::

    register_backend("MY_BACKEND", "my_pkg.my_backend", "MyBackend")

Concrete backend modules import their SDKs lazily, so only the extras you use
are required.
"""

from Qdislib.core.backend.base_backend import (
    BaseBackend,
    BackendConfig,
    Job,
    SampleJobResult,
    ExpectationJobResult,
    SupportedCircuit,
)
from Qdislib.core.backend.circuit_converter import CircuitConverter
from Qdislib.core.backend.registry import (
    register_backend,
    build_backend,
    available_backends,
    DEFAULT_BACKEND_NAME,
    DEFAULT_NUM_QUBITS,
)
from Qdislib.core.backend.unified import QdBackEnd

__all__ = [
    "QdBackEnd",
    "BaseBackend",
    "BackendConfig",
    "Job",
    "SampleJobResult",
    "ExpectationJobResult",
    "SupportedCircuit",
    "CircuitConverter",
    "register_backend",
    "build_backend",
    "available_backends",
    "DEFAULT_BACKEND_NAME",
    "DEFAULT_NUM_QUBITS",
]
