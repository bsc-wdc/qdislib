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

"""Coverage for runtime shim, backend helpers, parsers and find_cut modes."""

import importlib.util
import os
import tempfile
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None


class ExceptionsTest(BaseTimedTestCase):
    def test_a_pickled_exception_keeps_its_message(self):
        # what happens to an exception raised in a PyCOMPSs task
        import pickle

        from Qdislib.utils.exceptions import InvalidCutError

        error = InvalidCutError("bad cut")
        again = pickle.loads(pickle.dumps(error))
        self.assertIsInstance(again, InvalidCutError)
        self.assertEqual(str(again), str(error))


class LogLevelFromEnvironmentTest(BaseTimedTestCase):
    def test_a_bad_level_does_not_stop_the_import(self):
        import warnings
        from unittest import mock

        from Qdislib.utils import logging_config

        with mock.patch.dict(os.environ, {"QDISLIB_LOG_LEVEL": "VERBOSE"}):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                logging_config._WARNED.clear()
                logging_config._init_from_env()
        self.assertTrue(any("VERBOSE" in str(w.message) for w in caught))

    def test_set_log_level_rejects_a_bad_level(self):
        import Qdislib as qd

        with self.assertRaises(ValueError):
            qd.set_log_level("VERBOSE")


class NativeFormsTest(BaseTimedTestCase):
    @unittest.skipUnless(importlib.util.find_spec("pennylane"), "pennylane required")
    def test_a_native_pennylane_circuit_with_pycompss_installed(self):
        # with PyCOMPSs installed it was taken for the future of a Subcircuit
        from unittest import mock

        import Qdislib as qd
        from Qdislib.utils.graph_pennylane import PennylaneCircuit

        circuit = PennylaneCircuit([("X", (0,), [])], 1)
        with mock.patch("Qdislib.utils.backends.pycompss_available", True):
            self.assertEqual(qd.expectation_value(circuit, shots=10, seed=1), -1.0)


class RuntimeShimTest(BaseTimedTestCase):
    """The PyCOMPSs shim must provide working no-op fallbacks when absent."""

    def test_shim_exports(self):
        from Qdislib.utils import runtime

        self.assertIsInstance(runtime.pycompss_available, bool)

        # task() is a decorator factory; in serial mode it returns the function.
        @runtime.task(returns=1)
        def add(a, b):
            return a + b

        self.assertEqual(add(2, 3), 5)
        # compss_wait_on is identity in serial mode (and a real sync under COMPSs).
        if not runtime.pycompss_available:
            self.assertEqual(runtime.compss_wait_on(42), 42)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class BackendHelpersTest(BaseTimedTestCase):
    """convert_to_dag / iter_gates / num_qubits / change_basis_dag."""

    def _circuit(self):
        from qibo import models, gates

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.RX(2, 0.5))
        return c

    def test_convert_unknown_type_returns_unchanged(self):
        from Qdislib.utils.backends import convert_to_dag

        sentinel = object()
        self.assertIs(convert_to_dag(sentinel), sentinel)

    def test_iter_gates_unsupported_type_raises(self):
        from Qdislib.utils.backends import iter_gates

        with self.assertRaises(TypeError):
            iter_gates(12345)

    def test_num_qubits(self):
        from Qdislib.utils.backends import num_qubits

        self.assertEqual(num_qubits(self._circuit()), 3)

    def test_iter_gates_names(self):
        from Qdislib.utils.backends import iter_gates

        names = [name for name, _ in iter_gates(self._circuit())]
        self.assertEqual(names, ["H", "CZ", "RX"])

    def test_change_basis_dag_appends_rotations(self):
        from Qdislib.utils.backends import convert_to_dag, change_basis_dag

        dag = convert_to_dag(self._circuit())
        before = dag.number_of_nodes()
        change_basis_dag(dag, "XYZ")  # X -> +1 node (H), Y -> +2 nodes (SDG, H)
        self.assertEqual(dag.number_of_nodes(), before + 3)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class QsimParserTest(BaseTimedTestCase):
    """parse_qsim builds the file's unitary, with qsim's own gate definitions."""

    def _parse(self, content, **kwargs):
        from Qdislib.utils.circuit import parse_qsim

        with tempfile.NamedTemporaryFile("w", suffix=".qsim", delete=False) as handle:
            handle.write(content)
            path = handle.name
        try:
            return parse_qsim(path, **kwargs)
        finally:
            os.remove(path)

    def _assert_same_gate(self, line, width, expected):
        import numpy as np

        unitary = np.array(self._parse(f"{width}\n{line}\n").unitary())
        overlap = abs(np.trace(unitary.conj().T @ expected)) / 2**width
        self.assertAlmostEqual(overlap, 1.0, places=9, msg=line)

    def test_one_qubit_gates_are_qsims(self):
        import numpy as np

        # the matrices of qsim's lib/gates_qsim.h
        h, r = (1 + 1j) / 2, 1 / np.sqrt(2)
        cases = {
            "0 x_1_2 0": [[h, (1 - 1j) / 2], [(1 - 1j) / 2, h]],
            "0 y_1_2 0": [[h, -h], [h, h]],
            "0 hz_1_2 0": [[h, -1j * r], [r, h]],
            "0 rz 0 0.7": [[np.exp(-0.35j), 0], [0, np.exp(0.35j)]],
            "0 rx 0 0.7": [
                [np.cos(0.35), -1j * np.sin(0.35)],
                [-1j * np.sin(0.35), np.cos(0.35)],
            ],
            "0 ry 0 0.7": [
                [np.cos(0.35), -np.sin(0.35)],
                [np.sin(0.35), np.cos(0.35)],
            ],
            "0 h 0": [[r, r], [r, -r]],
            "0 t 0": [[1, 0], [0, np.exp(0.25j * np.pi)]],
            "0 s 0": [[1, 0], [0, 1j]],
        }
        for line, matrix in cases.items():
            self._assert_same_gate(line, 1, np.array(matrix))
        theta, phi = 0.4, 0.9
        c, s = np.cos(phi / 2), np.sin(phi / 2)
        axis = np.cos(theta) * np.array([[0, 1], [1, 0]]) + np.sin(theta) * np.array(
            [[0, -1j], [1j, 0]]
        )
        self._assert_same_gate(
            f"0 rxy 0 {theta} {phi}", 1, c * np.eye(2) - 1j * s * axis
        )

    def test_two_qubit_gates_are_qsims(self):
        import numpy as np

        theta, phi = 0.6, 0.3
        c, s = np.cos(theta), np.sin(theta)
        cases = {
            f"0 fs 0 1 {theta} {phi}": [
                [1, 0, 0, 0],
                [0, c, -1j * s, 0],
                [0, -1j * s, c, 0],
                [0, 0, 0, np.exp(-1j * phi)],
            ],
            f"0 cp 0 1 {phi}": np.diag([1, 1, 1, np.exp(-1j * phi)]),
            "0 cz 0 1": np.diag([1, 1, 1, -1]),
            "0 cnot 0 1": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
            "0 sw 0 1": [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
            "0 is 0 1": [[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]],
            "0 id2 0 1": np.eye(4),
        }
        for line, matrix in cases.items():
            self._assert_same_gate(line, 2, np.array(matrix))

    def test_depth_is_by_time_not_by_line(self):
        c = self._parse("2\n1 h 0\n5 x 1\n0 h 1\n", depth=2)
        self.assertEqual([g.__class__.__name__ for g in c.queue], ["H", "H"])

    def test_what_it_cannot_read_is_refused(self):
        for content in (
            "two\n0 h 0\n",
            "",
            "2\n0 m 0 1\n",
            "2\n0 mat1 0 1 0 0 0 0 0 0 1\n",
            "2\n0 fs 0 1\n",
            "2\n0 rz 0\n",
            "2\n0\n",
        ):
            with self.assertRaises(ValueError, msg=content):
                self._parse(content)

    def test_a_cut_of_its_fsim_gives_the_circuits_value(self):
        import numpy as np
        import Qdislib as qd

        c = self._parse(
            "4\n0 hz_1_2 0\n0 x_1_2 1\n0 y_1_2 2\n0 rz 3 0.4\n0 x_1_2 3\n"
            "1 fs 0 1 0.6 0.3\n2 fs 1 2 0.5 0.2\n3 fs 2 3 0.4 0.7\n"
            "4 rz 1 0.9\n4 hz_1_2 2\n5 y_1_2 1\n5 x_1_2 2\n"
        )
        probabilities = np.abs(np.array(c().state())) ** 2
        # qubit 0 is qibo's most significant bit
        signs = np.array(
            [(-1) ** (((i >> 3) & 1) + (i & 1)) for i in range(16)], dtype=float
        )
        exact = float(probabilities @ signs)
        self.assertGreater(abs(exact), 1e-3)  # not a value any mistake would give too
        for backend in ("qibo", "qiskit"):
            value = qd.circuit_cutting(
                c, {"hadamard": ["FSIM_7"]}, "ZIIZ", backend=backend, shots=None
            )
            self.assertAlmostEqual(value, exact, places=8, msg=backend)


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class FindCutModesTest(BaseTimedTestCase):
    """find_cut should return sensible cuts in each mode and validate input."""

    def _circuit(self):
        import numpy as np
        from qibo import models, gates

        c = models.Circuit(4)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.RY(2, np.pi / 5))
        c.add(gates.CZ(2, 3))
        c.add(gates.RX(1, np.pi / 4))
        return c

    def test_gate_cut_mode(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        cut = find_cut(self._circuit(), gate_cut=True, wire_cut=False)
        self.assertTrue(cut)
        self.assertTrue(all(isinstance(g, str) for g in cut))

    def test_wire_cut_mode(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        cut = find_cut(self._circuit(), gate_cut=False, wire_cut=True)
        self.assertTrue(cut)
        self.assertTrue(all(isinstance(w, tuple) and len(w) == 2 for w in cut))

    def test_both_modes_returns_a_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        self.assertTrue(find_cut(self._circuit()))

    def test_requires_a_mode(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        with self.assertRaises(ValueError):
            find_cut(self._circuit(), gate_cut=False, wire_cut=False)

    def test_seed_is_reproducible(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        a = find_cut(self._circuit(), gate_cut=True, wire_cut=False, seed=123)
        b = find_cut(self._circuit(), gate_cut=True, wire_cut=False, seed=123)
        self.assertEqual(a, b)

    def test_unknown_implementation_raises(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        with self.assertRaises(ValueError):
            find_cut(self._circuit(), implementation="does-not-exist")

    def test_cut_cost(self):
        import Qdislib as qd

        self.assertEqual(
            qd.cut_cost(["CZ_3", "CZ_5"]),
            {"num_cuts": 2, "kind": "gate", "base": 6, "terms": 36},
        )
        self.assertEqual(qd.cut_cost([("H_1", "CZ_2")])["terms"], 8)
        self.assertEqual(qd.cut_cost([])["terms"], 1)

    def test_forbidden_gates_are_not_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        circ = self._circuit()
        base = find_cut(circ, gate_cut=True, wire_cut=False, seed=0)
        self.assertTrue(base)
        # Forbidding the gates it would pick must yield a cut that avoids them.
        avoided = find_cut(
            circ, gate_cut=True, wire_cut=False, seed=0, forbidden_gates=base
        )
        self.assertFalse(set(avoided) & set(base))

    def test_top_n_ranks_candidates(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        ranked = find_cut(self._circuit(), seed=0, top_n=3)
        self.assertIsInstance(ranked, list)
        self.assertLessEqual(len(ranked), 3)
        self.assertTrue(all(set(c) == {"cut", "kind", "score"} for c in ranked))
        scores = [c["score"] for c in ranked]
        self.assertEqual(scores, sorted(scores))  # best score first

    def test_forbidden_qubits_are_not_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.utils.backends import iter_gates

        circ = self._circuit()
        gate_qubits = {
            f"{name}_{i}": set(qs)
            for i, (name, qs) in enumerate(iter_gates(circ), start=1)
        }
        cut = find_cut(
            circ, gate_cut=True, wire_cut=False, seed=0, forbidden_qubits=[0]
        )
        for gate in cut:  # no cut gate may act on the forbidden qubit
            self.assertNotIn(0, gate_qubits.get(gate, set()))

    def test_wire_cut_counts_shared_qubits_as_multiple_wires(self):
        """A wire cut between two gates sharing both qubits weighs two cuts, not one.

        Regression: the gate-graph edge collapses the two parallel wires into a
        single edge, so cut-size scoring must use the per-edge shared-qubit weight
        -- otherwise a two-wire cut is costed as one and mis-ranked as cheaper.
        """
        from qibo import models, gates
        from Qdislib.core.find_cut.find_cut import (
            _circuit_to_gate_graph,
            _gate_edge_weights,
            _rank_candidates,
        )

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.CZ(1, 2))  # shares BOTH qubits 1 and 2 with the previous gate

        weights = _gate_edge_weights(_circuit_to_gate_graph(c))
        self.assertEqual(weights[("CZ_3", "CZ_4")], 2)  # two shared wires
        self.assertEqual(weights[("CZ_2", "CZ_3")], 1)  # one shared wire

        # Same partition shape, so only the cut-size term differs: the two-wire cut
        # must score strictly higher (worse) than the one-wire cut.
        one_wire = {"m": ([{0}, {1, 2}], [("CZ_2", "CZ_3")])}
        two_wire = {"m": ([{0}, {1, 2}], [("CZ_3", "CZ_4")])}
        s_one = _rank_candidates(one_wire, edge_weights=weights)[0][1]
        s_two = _rank_candidates(two_wire, edge_weights=weights)[0][1]
        self.assertLess(s_one, s_two)


class IdentityCutGuardTest(BaseTimedTestCase):
    """The identity-on-a-cut-qubit guard: CPU/GPU decouple the sign; QPU cannot."""

    def _cut_dag_identity_feeds_multiqubit(self):
        """A minimal cut DAG where a cut qubit carries an ``I`` observable AND feeds
        a multi-qubit gate after the cut -- the only case the guard rejects."""
        import networkx as nx

        dag = nx.DiGraph()
        dag.add_node("SUBS_0", gate="SUBS", qubits=(1,))
        dag.add_node("CZ_9", gate="CZ", qubits=(1, 2))
        dag.add_node("obsI", gate="Observable I", qubits=(1,))
        dag.add_edge("SUBS_0", "CZ_9")
        dag.add_edge("CZ_9", "obsI")
        return dag

    def test_gpu_path_is_allowed(self):
        # Regression: GPU routes to Aer's GPU simulator and runs the SAME decoupled
        # joint-counts reconstruction as CPU, so the guard must NOT raise for GPU.
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _assert_identity_observable_supported,
        )

        _assert_identity_observable_supported(  # must not raise
            self._cut_dag_identity_feeds_multiqubit(), software="qiskit", gpu=True
        )

    def test_qpu_path_still_raises(self):
        from Qdislib import InvalidCutError
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _assert_identity_observable_supported,
        )

        with self.assertRaises(InvalidCutError):
            _assert_identity_observable_supported(
                self._cut_dag_identity_feeds_multiqubit(),
                software="qiskit",
                qpu=True,
            )
