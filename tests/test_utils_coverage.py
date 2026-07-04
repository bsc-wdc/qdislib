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
    """parse_qsim must be deterministic and map gates correctly."""

    def test_parse_qsim_deterministic(self):
        from Qdislib.utils.circuit import parse_qsim

        content = "3\n0 rz 0 0.5\n1 x_1_2 1\n2 hz_1_2 0\n3 fs 1 2\n"
        with tempfile.NamedTemporaryFile("w", suffix=".qsim", delete=False) as handle:
            handle.write(content)
            path = handle.name
        try:
            c1 = parse_qsim(path)
            c2 = parse_qsim(path)
            self.assertEqual(c1.nqubits, 3)
            # Deterministic: same gate sequence on both parses (no random angles).
            self.assertEqual(
                [g.__class__.__name__ for g in c1.queue],
                [g.__class__.__name__ for g in c2.queue],
            )
        finally:
            os.remove(path)


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
