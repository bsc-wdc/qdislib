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

"""Lazy, user-extensible registry of Qdislib backends.

Backends are registered by name -> (module path, class name) and imported only
on first use, so a backend's optional third-party dependencies (qiskit-aer,
iqm, ...) are not required until that backend is actually built.
"""

import importlib
import logging

from Qdislib.core.backend.base_backend import BaseBackend, BackendConfig

logger = logging.getLogger(__name__)

DEFAULT_BACKEND_NAME = "AER_BACKEND"
DEFAULT_NUM_QUBITS = 30

_BACKEND_REGISTRY: dict[str, tuple[str, str]] = {
    "AER_BACKEND": ("Qdislib.core.backend.aer_backend", "AerBackend"),
    "QIBO_BACKEND": ("Qdislib.core.backend.qibo_backend", "QiboBackend"),
    "CUDAQ_BACKEND": ("Qdislib.core.backend.cudaq_backend", "CudaqBackend"),
    "IQM_BACKEND": ("Qdislib.core.backend.iqm_backend", "IQMBackend"),
}


def register_backend(name: str, module_path: str, class_name: str) -> None:
    """Register a backend without importing it.

    Intended for user-defined backends or third-party extensions::

        register_backend("MY_BACKEND", "my_package.my_backend", "MyBackend")
    """
    _BACKEND_REGISTRY[name] = (module_path, class_name)
    logger.debug("Backend '%s' registered (imported on first use)", name)


def available_backends() -> list[str]:
    """Return the names of all registered backends."""
    return list(_BACKEND_REGISTRY)


#: Maps the algorithms' ``software=`` strings to registered backend names, so a
#: cutting call can also name any user-registered backend directly.
_SOFTWARE_ALIAS = {
    "qibo": "QIBO_BACKEND",
    "qiskit": "AER_BACKEND",
    "cudaq": "CUDAQ_BACKEND",
}


def build_execution_backend(
    software: str | None = None, use_gpu: bool = False, num_qubits: int = 1, **extra
) -> "BaseBackend":
    """Build the backend that evaluates a cut subcircuit for a ``software`` choice.

    This is the single place mapping the cutting algorithms' ``software=`` knob to
    a concrete backend, so the algorithms never branch on the backend name and a
    user can pass the name of any backend they registered via
    :func:`register_backend`. Extra keyword arguments (e.g. ``platform`` for Qibo)
    flow into the backend's :class:`BackendConfig`.
    """
    size = max(int(num_qubits or 1), 1)
    if use_gpu:
        extra = {**extra, "use_gpu": True}
        name = "AER_BACKEND"
    else:
        name = _SOFTWARE_ALIAS.get(software, software or DEFAULT_BACKEND_NAME)
    return build_backend(name, BackendConfig(num_qubits=size, extra=extra))


def build_backend(name: str, config: BackendConfig) -> BaseBackend:
    """Lazily import and instantiate a registered backend.

    The backend module -- and its third-party dependencies -- are imported only
    at this point, not at registry import time.
    """
    if name not in _BACKEND_REGISTRY:
        raise ValueError(
            f"Backend '{name}' is not registered. Available: {available_backends()}"
        )

    module_path, class_name = _BACKEND_REGISTRY[name]
    try:
        module = importlib.import_module(module_path)
    except ImportError as exc:
        raise ImportError(
            f"Could not import backend '{name}' from '{module_path}'. "
            f"Its optional dependency may be missing."
        ) from exc

    cls: type[BaseBackend] = getattr(module, class_name)
    return cls(config)
