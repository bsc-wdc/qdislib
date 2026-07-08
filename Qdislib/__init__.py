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

"""Main Quantum Distributed Library package."""

import logging

__version__ = "1.0.0"

# Attach a NullHandler so the library never emits "no handlers" warnings and
# stays silent unless the host application configures logging itself. Diagnostic
# output is opt-in via qd.set_log_level(...) or the QDISLIB_LOG_LEVEL env var.
logging.getLogger(__name__).addHandler(logging.NullHandler())

from Qdislib.utils import logging_config as _logging_config

_logging_config._init_from_env()

from Qdislib.api import *

# Curated public surface. This is what `from Qdislib import *`, tab-completion and
# the generated docs advertise; internal plumbing (DAG converters, job/cache
# strategy objects, the QPU listener, module aliases) stays reachable via explicit
# attribute access but is kept out of the public namespace. Keep grouped by role.
__all__ = [
    # Cutting algorithms -- all-in-one.
    "gate_cutting",
    "gate_cutting_hadamard",
    "wire_cutting",
    "gate_wire_cutting",
    # Cutting -- explicit subcircuits + reconstruction (the three-step workflow).
    "gate_cutting_subcircuits",
    "gate_cutting_subcircuit_reconstruction",
    "gate_cutting_hadamard_subcircuits",
    "gate_cutting_hadamard_subcircuit_reconstruction",
    "wire_cutting_subcircuits",
    "wire_cutting_subcircuit_reconstruction",
    # Cutting -- sampling-based reconstruction.
    "gate_cutting_sampling",
    "gate_cutting_sampling_subcircuits",
    "gate_cutting_sampling_reconstruction",
    "gate_cutting_hadamard_sampling",
    "gate_cutting_hadamard_sampling_subcircuits",
    "wire_cutting_sampling",
    "wire_cutting_sampling_subcircuits",
    # Subcircuit containers.
    "Subcircuit",
    "CutSubcircuits",
    "HadamardCutSubcircuits",
    "SampledCutSubcircuits",
    "WireSampledCutSubcircuits",
    # Executing subcircuits / evaluating observables.
    "expectation_value",
    "sample_subcircuit",
    "analytical_solution",
    # Finding and inspecting cuts.
    "find_cut",
    "draw_cut",
    "cut_cost",
    # Cost estimation.
    "estimate",
    # Execution backends.
    "available_backends",
    "build_backend",
    "build_execution_backend",
    "register_backend",
    "QdBackEnd",
    "BaseBackend",
    "BackendConfig",
    # Semantic circuit cache.
    "default_circuit_cache",
    "semantic_key",
    "CircuitCache",
    "InMemoryCircuitCache",
    "LMDBCircuitCache",
    "RedisCircuitCache",
    # Circuit I/O.
    "from_qasm",
    # Graph helpers.
    "generate_random_graph",
    # Logging.
    "set_log_level",
    # Exceptions.
    "QdislibException",
    "UnsupportedBackendError",
    "InvalidCutError",
    "QPUExecutionError",
    "QPUTimeoutError",
]
