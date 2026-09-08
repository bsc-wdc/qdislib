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

"""Guards for round-tripping circuits through OpenQASM 2.

OpenQASM 2 has no way to express a measurement that is followed by further
operations on the same wire, and the SDK exporters do not all say so. Qibo's
``Circuit.to_qasm()`` emits a *terminal* measurement as ``creg`` + ``measure``
but **silently drops a mid-circuit one**::

    c.add(gates.RY(0, 0.8)); c.add(gates.CZ(0, 1))
    c.add(gates.M(1))                      # <- disappears
    c.add(gates.RZ(1, -math.pi / 2))
    c.to_qasm()                            # identical to the circuit without M(1)

That matters here because gate cutting produces exactly such fragments: two of
the six quasiprobability configurations measure the cut qubit mid-circuit and
re-prepare it. Exported to QASM they become indistinguishable from the
configuration that does neither -- so anything keyed on, or reconstructed from,
that QASM is silently wrong rather than loudly broken.

:func:`qasm_drops_measurements` compares what went in with what came out instead
of trying to predict each exporter's rules, so it stays correct as they change.
"""

import logging
import re

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_measurement_count",
    "qasm_measurement_count",
    "qasm_drops_measurements",
]

#: ``measure q[0] -> c[0];`` and the whole-register ``measure q -> c;`` form.
_MEASURE_RE = re.compile(r"^\s*measure\b", re.MULTILINE)


def circuit_measurement_count(circuit, software="qibo"):
    """Number of measured qubits in ``circuit``, or ``None`` if we cannot tell.

    ``None`` means "no opinion" -- the caller should not treat it as zero.
    """
    try:
        if software == "qibo":
            queue = getattr(circuit, "queue", None)
            if queue is None:
                return None
            return sum(len(gate.qubits) for gate in queue if type(gate).__name__ == "M")
        if software == "qiskit":
            return sum(
                len(instruction.qubits)
                for instruction in circuit.data
                if instruction.operation.name == "measure"
            )
    except Exception:  # a guard must never be the thing that breaks a run
        logger.debug("Could not count measurements on the circuit", exc_info=True)
    return None


def qasm_measurement_count(qasm):
    """Number of ``measure`` statements in a QASM string."""
    return len(_MEASURE_RE.findall(qasm))


def qasm_drops_measurements(circuit, qasm, software="qibo"):
    """True if exporting ``circuit`` to ``qasm`` lost a measurement.

    Compares the measured-qubit count before and after rather than inspecting
    gate order, so it catches whatever an exporter drops. Returns ``False`` when
    the circuit's measurements cannot be counted -- an unknown SDK must not make
    a working path start failing.

    :param circuit: The source circuit.
    :param qasm: Its OpenQASM-2 serialisation.
    :param software: SDK the circuit belongs to (``"qibo"``, ``"qiskit"``, ...).
    :return: ``True`` when the QASM has fewer measurements than the circuit.
    """
    expected = circuit_measurement_count(circuit, software)
    if expected is None:
        return False
    # A whole-register `measure q -> c;` covers every qubit in one statement, so
    # only a strict shortfall of statements is evidence of loss.
    return expected > 0 and qasm_measurement_count(qasm) == 0
