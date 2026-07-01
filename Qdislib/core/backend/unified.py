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

"""The unified backend: one object the algorithms ask to sample / estimate.

:class:`QdBackEnd` holds one or more :class:`~Qdislib.core.backend.base_backend.BaseBackend`
instances and routes each job to the most suitable one -- shortest queue first,
qubit-count compatible, with fallback to the next candidate on failure. This is
what decouples the cutting algorithms from any specific SDK or device.
"""

import logging

from Qdislib.core.backend.base_backend import (
    BaseBackend,
    BackendConfig,
    SampleJobResult,
    ExpectationJobResult,
    SupportedCircuit,
)
from Qdislib.core.backend.registry import (
    build_backend,
    DEFAULT_BACKEND_NAME,
    DEFAULT_NUM_QUBITS,
)

try:
    from pycompss.api.task import task
except ImportError:

    def task(*args, **kwargs):
        def decorator(func):
            return func

        return decorator


logger = logging.getLogger(__name__)


class QdBackEnd:
    """Unified, load-balancing backend for Qdislib.

    Parameters
    ----------
    qpu_dict:
        Maps backend names to their config, e.g.::

            QdBackEnd({
                "AER_BACKEND": {"num_qubits": 30},
                "IQM_BACKEND": {"num_qubits": 20, "url": "...", "tokens_path": "..."},
            })

        With no argument, a single default Aer backend is used.
    """

    def __init__(self, qpu_dict: dict | None = None) -> None:
        self._backends: list[BaseBackend] = []
        if qpu_dict is None:
            cfg = BackendConfig(num_qubits=DEFAULT_NUM_QUBITS, extra={})
            backend = build_backend(DEFAULT_BACKEND_NAME, cfg)
            self._backends.append(backend)
            logger.info("No backend specified; using default: %s", backend.name)
        else:
            for name, data in qpu_dict.items():
                extra = {k: v for k, v in data.items() if k != "num_qubits"}
                cfg = BackendConfig(
                    num_qubits=data.get("num_qubits", DEFAULT_NUM_QUBITS), extra=extra
                )
                self._backends.append(build_backend(name, cfg))
                logger.info("Registered backend: %s (%d qubits)", name, cfg.num_qubits)

    @property
    def backends(self) -> list[BaseBackend]:
        return list(self._backends)

    @task(returns=1)
    def sample_circuit(
        self,
        circuit: SupportedCircuit,
        shots: int = 1024,
        transpile: bool = True,
        backend_name: str | None = None,
    ) -> SampleJobResult:
        """Sample the circuit on the most suitable backend (or a named one)."""
        last_exc: Exception | None = None
        for backend in self._select_backends(circuit, backend_name):
            try:
                logger.debug("Sampling on '%s'", backend.name)
                return backend.run_samples(circuit, transpile, shots)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Sampling failed on '%s': %s", backend.name, exc)
                last_exc = exc
        raise RuntimeError(f"No backend could sample the circuit. Last error: {last_exc}")

    @task(returns=1)
    def get_expectation_value(
        self,
        circuit: SupportedCircuit,
        observable: str,
        shots: int = 1024,
        transpile: bool = True,
        backend_name: str | None = None,
    ) -> ExpectationJobResult:
        """Estimate ``<observable>``; falls back past backends without estimator support."""
        last_exc: Exception | None = None
        for backend in self._select_backends(circuit, backend_name):
            try:
                logger.debug("Estimating on '%s'", backend.name)
                return backend.run_estimator(circuit, transpile, observable, shots)
            except NotImplementedError:
                logger.debug("'%s' has no estimator, skipping", backend.name)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Estimation failed on '%s': %s", backend.name, exc)
                last_exc = exc
        raise RuntimeError(f"No backend could estimate the observable. Last error: {last_exc}")

    def _select_backends(
        self, circuit: SupportedCircuit, backend_name: str | None
    ) -> list[BaseBackend]:
        """Candidate backends, shortest queue first (or the single named one)."""
        if backend_name is not None:
            matches = [b for b in self._backends if b.name == backend_name]
            if not matches:
                raise ValueError(f"Backend '{backend_name}' not found")
            return matches

        compatible = [b for b in self._backends if b.can_run(circuit)]
        if not compatible:
            raise ValueError("No registered backend is large enough for this circuit")
        return sorted(compatible, key=lambda b: b.queue_length)
