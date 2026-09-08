#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""OpenQASM 2 cannot carry a measurement that has operations after it.

Qibo's exporter drops such a measurement silently, which makes two different
gate-cut fragments serialise identically. These tests pin the guard that turns
that silent corruption into a skipped cache entry or a clear error.
"""

import importlib.util
import math
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None


def _circuit(mid_measure):
    from qibo import gates, models

    circuit = models.Circuit(2)
    circuit.add(gates.RY(0, 0.8))
    circuit.add(gates.CZ(0, 1))
    if mid_measure:
        circuit.add(gates.M(1))
        circuit.add(gates.RZ(1, -math.pi / 2))
    else:
        circuit.add(gates.RZ(1, -math.pi / 2))
    return circuit


class QasmCountingTest(BaseTimedTestCase):
    """Counting works on strings alone, so it needs no SDK."""

    def test_counts_measure_statements(self):
        from Qdislib.utils.qasm import qasm_measurement_count

        self.assertEqual(qasm_measurement_count("qreg q[2];\nry(0.8) q[0];\n"), 0)
        self.assertEqual(
            qasm_measurement_count("measure q[0] -> c[0];\nmeasure q[1] -> c[1];\n"), 2
        )

    def test_unknown_sdk_is_not_reported_as_lossy(self):
        """A guard must not make an unrecognised circuit start failing."""
        from Qdislib.utils.qasm import qasm_drops_measurements

        self.assertFalse(qasm_drops_measurements(object(), "qreg q[1];", "qibo"))


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class QiboExportTest(BaseTimedTestCase):
    def test_mid_circuit_measurement_is_detected_as_dropped(self):
        from Qdislib.utils.qasm import qasm_drops_measurements

        circuit = _circuit(mid_measure=True)
        qasm = circuit.to_qasm()
        # The bug this guards: the export is byte-identical to the circuit that
        # never had the measurement at all.
        self.assertEqual(qasm, _circuit(mid_measure=False).to_qasm())
        self.assertTrue(qasm_drops_measurements(circuit, qasm, "qibo"))

    def test_terminal_measurement_survives_and_is_not_flagged(self):
        from qibo import gates, models

        from Qdislib.utils.qasm import qasm_drops_measurements

        circuit = models.Circuit(2)
        circuit.add(gates.RY(0, 0.8))
        circuit.add(gates.M(1))
        qasm = circuit.to_qasm()
        self.assertIn("measure", qasm)
        self.assertFalse(qasm_drops_measurements(circuit, qasm, "qibo"))

    def test_circuit_without_measurements_is_not_flagged(self):
        from Qdislib.utils.qasm import qasm_drops_measurements

        circuit = _circuit(mid_measure=False)
        self.assertFalse(qasm_drops_measurements(circuit, circuit.to_qasm(), "qibo"))

    def test_semantic_key_refuses_to_key_a_lossy_circuit(self):
        """A colliding key would serve one fragment's value for another."""
        from Qdislib.utils.circuit_cache import semantic_key

        self.assertIsNone(semantic_key(_circuit(mid_measure=True), "qibo"))


@unittest.skipUnless(QIBO_AVAILABLE and QISKIT_AVAILABLE, "qibo and qiskit both needed")
class ConverterGuardTest(BaseTimedTestCase):
    def test_conversion_refuses_a_dropped_measurement(self):
        from Qdislib.core.backend.circuit_converter import CircuitConverter
        from Qdislib.utils.exceptions import InvalidCutError

        with self.assertRaises(InvalidCutError):
            CircuitConverter.to_qiskit(_circuit(mid_measure=True))

    def test_clean_circuit_still_converts(self):
        from Qdislib.core.backend.circuit_converter import CircuitConverter

        converted = CircuitConverter.to_qiskit(_circuit(mid_measure=False))
        self.assertEqual(converted.num_qubits, 2)


if __name__ == "__main__":
    unittest.main()
