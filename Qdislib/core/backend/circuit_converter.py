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

"""Convert circuits between the supported source formats (Qiskit / Qibo).

Qiskit and Qibo interconvert cleanly through OpenQASM 2.0. CUDA-Q kernels have
no stable QASM round-trip, so Qdislib feeds CUDA-Q via its backend-neutral DAG
instead; :meth:`CircuitConverter.num_qubits` still understands a CUDA-Q kernel
and the ``(dag, num_qubits)`` payload the cutting code uses for it.
"""

import logging
from typing import Any

logger = logging.getLogger(__name__)


def _qiskit_circuit():
    """The Qiskit ``QuantumCircuit`` class, or ``None`` if Qiskit isn't installed.

    Imported lazily so ``import Qdislib`` works with only one backend SDK present
    (a Qibo-only / Qiskit-only / CUDA-Q-only install). Conversions that need the
    missing SDK then raise a clear ``TypeError`` instead of failing at import.
    """
    try:
        from qiskit import QuantumCircuit

        return QuantumCircuit
    except ImportError:
        return None


def _qibo_circuit():
    """The Qibo ``Circuit`` class, or ``None`` if Qibo isn't installed."""
    try:
        from qibo import Circuit

        return Circuit
    except ImportError:
        return None


class CircuitConverter:
    """Converts circuits between the supported formats and reports their size."""

    @staticmethod
    def to_qiskit(circuit: Any):
        qk, qb = _qiskit_circuit(), _qibo_circuit()
        if qk and isinstance(circuit, qk):
            return circuit
        if qk and qb and isinstance(circuit, qb):
            from qiskit.qasm2 import loads

            return loads(circuit.to_qasm())
        raise TypeError(f"Cannot convert {type(circuit).__name__} to a Qiskit circuit")

    @staticmethod
    def to_qibo(circuit: Any):
        qk, qb = _qiskit_circuit(), _qibo_circuit()
        if qb and isinstance(circuit, qb):
            return circuit
        if qb and qk and isinstance(circuit, qk):
            from qiskit.qasm2 import dumps

            return qb.from_qasm(dumps(circuit))
        raise TypeError(f"Cannot convert {type(circuit).__name__} to a Qibo circuit")

    @staticmethod
    def num_qubits(circuit: Any) -> int:
        """Qubit count for any supported circuit, including CUDA-Q payloads."""
        qk, qb = _qiskit_circuit(), _qibo_circuit()
        if qk and isinstance(circuit, qk):
            return circuit.num_qubits
        if qb and isinstance(circuit, qb):
            return circuit.nqubits
        # The CUDA-Q ``(dag, num_qubits)`` payload (second element is the count).
        if (
            isinstance(circuit, (list, tuple))
            and len(circuit) == 2
            and isinstance(circuit[1], int)
        ):
            return circuit[1]
        # A ``[circuit, observable]`` wrapper (as returned by the dag_to_circuit_* helpers).
        if (
            isinstance(circuit, (list, tuple))
            and circuit
            and not isinstance(circuit[0], int)
        ):
            return CircuitConverter.num_qubits(circuit[0])
        # Fall back to a common attribute (covers CUDA-Q kernels exposing one).
        for attr in ("num_qubits", "nqubits"):
            value = getattr(circuit, attr, None)
            if isinstance(value, int):
                return value
        # A built CUDA-Q kernel has no reliably-introspectable qubit count; the
        # cutting code feeds CUDA-Q via the (dag, num_qubits) payload above, so
        # this branch is only hit for a raw kernel. Report "size unknown" (0) so
        # queue bookkeeping / can_run don't block it, rather than raising.
        logger.debug(
            "num_qubits: unknown circuit type %s; assuming size 0",
            type(circuit).__name__,
        )
        return 0
