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

"""Tests for the public subcircuit workflow, backend helpers and error paths."""

import importlib.util
import random
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None
QISKIT_AVAILABLE = importlib.util.find_spec("qiskit") is not None

# Expectation value of <ZZ> for H(0), CZ(0,1), RY(0, 0.8); see analytical_solution.
_TOLERANCE = 0.2


def _example_circuit():
    """Build the small reference circuit used across these tests."""
    from qibo import models, gates

    circuit = models.Circuit(2)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.RY(0, 0.8))
    return circuit


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class PublicWorkflowTest(BaseTimedTestCase):
    """The recommended generate -> evaluate -> reconstruct workflow."""

    def test_gate_cutting_workflow_qibo(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        circuit = _example_circuit()
        analytical = qd.analytical_solution(circuit, "ZZ")

        subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")
        values = [qd.expectation_value(s, shots=20000) for s in subcircuits]
        reconstruction = subcircuits.reconstruct(values)

        self.assertEqual(subcircuits.number_cuts, 1)
        self.assertLess(abs(reconstruction - analytical), _TOLERANCE)

    def test_wire_cutting_workflow_qibo(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        circuit = _example_circuit()
        analytical = qd.analytical_solution(circuit, "ZZ")

        subcircuits = qd.wire_cutting_subcircuits(
            circuit, [("H_1", "CZ_2")], software="qibo"
        )
        values = [qd.expectation_value(s, shots=20000) for s in subcircuits]
        reconstruction = subcircuits.reconstruct(values)

        self.assertEqual(subcircuits.number_cuts, 1)
        self.assertLess(abs(reconstruction - analytical), _TOLERANCE)

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_gate_cutting_workflow_qiskit(self):
        import qibo
        import Qdislib as qd

        qibo.set_backend("numpy")
        circuit = _example_circuit()
        analytical = qd.analytical_solution(circuit, "ZZ")

        subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qiskit")
        values = [qd.expectation_value(s, shots=20000) for s in subcircuits]
        reconstruction = subcircuits.reconstruct(values)

        self.assertLess(abs(reconstruction - analytical), _TOLERANCE)

    @unittest.skipUnless(QISKIT_AVAILABLE, "qiskit not installed")
    def test_cnot_reconstructs_into_qiskit(self):
        """A Qibo ``CNOT`` must round-trip correctly into a Qiskit circuit.

        Guards the CNOT/CX backend-agnosticism fix in ``dag_to_circuit_qiskit``.
        """
        import qibo
        from qibo import models, gates
        import Qdislib as qd
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        qibo.set_backend("numpy")
        circuit = models.Circuit(2)
        circuit.add(gates.H(0))
        circuit.add(gates.CNOT(0, 1))
        circuit.add(gates.RY(0, 0.6))
        analytical = qd.analytical_solution(circuit, "ZZ")

        # Qibo DAG (gate name "CNOT") reconstructed as a Qiskit circuit.
        qiskit_subcircuit = dag_to_circuit_qiskit(circuit_qibo_to_dag(circuit), 2)
        value = qd.expectation_value(qiskit_subcircuit, shots=20000)

        self.assertLess(abs(value - analytical), _TOLERANCE)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class BackendHelpersTest(BaseTimedTestCase):
    """Backend-agnostic helpers in Qdislib.utils.backends."""

    def test_convert_to_dag_and_counts(self):
        import Qdislib as qd

        circuit = _example_circuit()
        dag = qd.convert_to_dag(circuit)

        # Three gates -> three DAG nodes; two qubits.
        self.assertEqual(dag.number_of_nodes(), 3)
        self.assertEqual(qd.num_qubits(circuit), 2)

    def test_iter_gates(self):
        import Qdislib as qd

        gates_list = qd.iter_gates(_example_circuit())
        names = [name for name, _ in gates_list]
        self.assertEqual(names, ["H", "CZ", "RY"])
        self.assertEqual(gates_list[1][1], (0, 1))  # CZ acts on qubits 0 and 1


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class ErrorPathTest(BaseTimedTestCase):
    """Negative tests: invalid input should fail clearly, not silently."""

    def test_invalid_software_raises(self):
        import Qdislib as qd
        from Qdislib import UnsupportedBackendError, QdislibException

        with self.assertRaises(UnsupportedBackendError):
            qd.gate_cutting_subcircuits(_example_circuit(), ["CZ_2"], software="bogus")
        # Backward compatible: it is also a ValueError and a QdislibException.
        try:
            qd.gate_cutting_subcircuits(_example_circuit(), ["CZ_2"], software="bogus")
        except UnsupportedBackendError as exc:
            self.assertIsInstance(exc, ValueError)
            self.assertIsInstance(exc, QdislibException)

    def test_invalid_cut_format_raises(self):
        import Qdislib as qd
        from Qdislib import InvalidCutError

        with self.assertRaises(InvalidCutError):
            qd.gate_cutting_subcircuits(_example_circuit(), "CZ_2")  # not a list

    def test_find_cut_requires_a_mode(self):
        import Qdislib as qd

        with self.assertRaises(ValueError):
            qd.find_cut(_example_circuit(), gate_cut=False, wire_cut=False)


class RandomGraphTest(BaseTimedTestCase):
    """The random graph generator must be deterministic and side-effect free."""

    def test_deterministic_and_no_global_state(self):
        import Qdislib as qd

        # Same seed -> identical edge sets.
        g1 = qd.generate_random_graph(8, probability=0.3, seed=7)
        g2 = qd.generate_random_graph(8, probability=0.3, seed=7)
        self.assertEqual(set(g1.edges()), set(g2.edges()))
        self.assertEqual(g1.number_of_nodes(), 8)

        # Generating a graph must not perturb the global RNG state.
        random.seed(123)
        expected = random.random()
        random.seed(123)
        qd.generate_random_graph(20, seed=99)
        self.assertEqual(random.random(), expected)
