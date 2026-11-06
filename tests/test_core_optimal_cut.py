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

"""Utils Circuit Unit tests."""

from tests import BaseTimedTestCase


class MetisCallTest(BaseTimedTestCase):
    def test_no_deprecated_arguments(self):
        import warnings

        import networkx as nx

        from Qdislib.core.find_cut.find_cut import _metis_partition

        graph = nx.MultiGraph()
        for q in range(4):
            graph.add_node(q, gate=q)
        for a, b, name in ((0, 1, "CZ_1"), (1, 2, "CZ_2"), (2, 3, "CZ_3")):
            graph.add_edge(a, b, gate=name)
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            clusters, cut = _metis_partition(graph, 2, None, None)
        self.assertEqual(len(clusters), 2)
        self.assertEqual(len(cut), 1)


def _chain_for_kinds():
    from qibo import gates, models

    circuit = models.Circuit(3)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.CZ(1, 2))
    return circuit


class CircuitTest(BaseTimedTestCase):

    def test_optimal_cut_wire_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo

        qibo.set_backend("numpy")

        circuit = models.Circuit(4)
        circuit.add(gates.H(0))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1, 2))
        circuit.add(gates.RY(1, np.pi / 5))
        circuit.add(gates.RZ(2, np.pi / 3))
        circuit.add(gates.RZ(3, np.pi / 2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 3, gate_cut=False)

        reconstrucion = wire_cutting(circuit, cut)
        analytical = analytical_solution(circuit, "ZZZZ")

        if len(cut) == 1 and abs(reconstrucion - analytical) < 0.2:
            self.assertTrue(True)

        else:
            self.assertTrue(False)

    def test_optimal_cut_gate_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo

        qibo.set_backend("numpy")

        circuit = models.Circuit(4)
        circuit.add(gates.H(0))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1, 2))
        circuit.add(gates.RY(1, np.pi / 5))
        circuit.add(gates.RZ(2, np.pi / 3))
        circuit.add(gates.RZ(3, np.pi / 2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 2, wire_cut=False)

        reconstrucion = gate_cutting(circuit=circuit, cuts=cut)
        analytical = analytical_solution(circuit, "ZZZZ")

        if cut == ["CZ_5"] and abs(reconstrucion - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_no_kind_of_cut_enabled_is_refused_in_every_mode(self):
        from Qdislib.core.find_cut.find_cut import find_cut

        circuit = _chain_for_kinds()
        for options in ({}, {"top_n": 3}):
            with self.assertRaisesRegex(ValueError, "Enable at least one"):
                find_cut(circuit, gate_cut=False, wire_cut=False, **options)

    def test_zero_budgets_are_limits_not_none(self):
        """bsc-wdc/qdislib#15: max_cuts=0 allows no cut; max_qubits=0 is rejected."""
        from Qdislib.core.find_cut.find_cut import find_cut
        from qibo import models, gates

        def build():
            circuit = models.Circuit(3)
            for q in range(3):
                circuit.add(gates.RY(q, 0.3))
            circuit.add(gates.CZ(0, 1))
            circuit.add(gates.RY(1, 0.4))
            circuit.add(gates.CZ(1, 2))
            return circuit

        self.assertTrue(find_cut(build(), max_qubits=2, seed=0))
        for kwargs in ({}, {"gate_cut": False}, {"wire_cut": False}):
            self.assertEqual(
                find_cut(build(), max_qubits=2, max_cuts=0, seed=0, **kwargs), []
            )
        # Mixed mode returns an over-budget cut with a warning, but 0 is a hard stop.
        self.assertEqual(
            find_cut(build(), max_qubits=2, max_cuts=0, seed=0, mixed=True),
            {"gate": [], "hadamard": [], "wire": []},
        )
        # q1 and q2 are two hops apart here, so CZ(1, 2) is a SparseCut candidate.
        star = [(0, 1), (0, 2)]
        self.assertTrue(find_cut(build(), coupling_map=star))
        self.assertEqual(find_cut(build(), max_cuts=0, coupling_map=star), [])
        with self.assertRaises(ValueError):
            find_cut(build(), max_qubits=0)
        with self.assertRaises(ValueError):
            find_cut(build(), max_cuts=-1)

    def test_barriers_are_never_cut(self):
        """bsc-wdc/qdislib#16: a Qiskit barrier is not a cut point in any mode."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.core.find_cut.find_cut import find_cut
        from qiskit import QuantumCircuit

        def build():
            circuit = QuantumCircuit(3)
            circuit.h(0)
            circuit.barrier()
            circuit.ry(0.4, 1)
            circuit.cz(0, 1)
            circuit.cz(1, 2)
            circuit.ry(0.2, 2)
            return circuit

        wire = find_cut(build(), gate_cut=False, max_qubits=2, seed=0)
        gate = find_cut(build(), wire_cut=False, max_qubits=2, seed=0)
        mixed = find_cut(build(), max_qubits=2, seed=0, mixed=True)
        self.assertTrue(wire)
        self.assertNotIn("BARRIER", repr((wire, gate, mixed)))
        # The cut must be usable: the barrier's index is kept, so names line up.
        wire_cutting(build(), wire, observable="ZZZ", backend="qiskit", shots=1000)

    def test_wire_mode_respects_max_qubits_and_forbidden_qubits(self):
        """bsc-wdc/qdislib#10: wire cuts keep pieces within max_qubits and off
        forbidden qubits."""
        from Qdislib.core.cutting_algorithms.wire_cutting import _expand_wire_cuts
        from Qdislib.core.find_cut.find_cut import _wire_cut_fits, find_cut
        from Qdislib.utils.backends import convert_to_dag
        from qibo import models, gates

        chain = models.Circuit(6)
        for q in range(5):
            chain.add(gates.CZ(q, q + 1))
        cut = find_cut(chain, gate_cut=False, wire_cut=True, max_qubits=3, seed=0)
        self.assertTrue(cut)
        self.assertTrue(_wire_cut_fits(convert_to_dag(chain), cut, 3, None, None))

        layered = models.Circuit(4)
        for _ in range(2):
            for q in range(4):
                layered.add(gates.H(q))
            for q in range(3):
                layered.add(gates.CZ(q, q + 1))
        cut = find_cut(
            layered, gate_cut=False, wire_cut=True, forbidden_qubits=[1], seed=0
        )
        self.assertTrue(cut)
        cut_qubits = {q for _, _, q in _expand_wire_cuts(convert_to_dag(layered), cut)}
        self.assertNotIn(1, cut_qubits)


class MetisWeightsTest(BaseTimedTestCase):
    """METIS sees parallel gates: a pair joined by 4 gates costs 4 to cut, not 1."""

    def test_cuts_the_light_links_of_a_ring(self):
        import networkx as nx

        from Qdislib.core.find_cut.find_cut import _metis_partition

        graph = nx.MultiGraph()
        graph.add_nodes_from(range(4))
        gate = iter(range(1, 100))
        for a, b, count in ((0, 1, 4), (1, 2, 1), (2, 3, 4), (3, 0, 1)):
            for _ in range(count):
                graph.add_edge(a, b, gate=f"CZ_{next(gate)}")
        clusters, cut = _metis_partition(graph, 2, None, None)
        self.assertEqual(sorted(map(sorted, clusters)), [[0, 1], [2, 3]])
        self.assertEqual(len(cut), 2)
        # the budget counts the gates cut
        self.assertEqual(_metis_partition(graph, 2, None, 1), ([], []))


class FindCutSkipsWhatCannotBeCutTest(BaseTimedTestCase):
    """A candidate pair of gates not next to each other on any qubit is no wire cut: the
    search leaves it out, and does not stop with its error."""

    def test_a_candidate_that_is_no_wire_cut_is_skipped(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        qc = QuantumCircuit(8)
        qc.ry(0.6, range(8))
        for q in range(0, 7, 2):
            qc.cz(q, q + 1)
        for q in range(1, 7, 2):
            qc.cz(q, q + 1)
        qc.ry(0.5, range(8))
        for q in range(0, 7, 2):
            qc.cz(q, q + 1)
        for q in range(1, 7, 2):
            qc.cz(q, q + 1)
        for seed in (0, 1, 2):
            self.assertEqual(
                qd.find_cut(qc, max_qubits=2, seed=seed, gate_cut=False), []
            )


class MetisMissingTest(BaseTimedTestCase):
    """pymetis has no wheel on Linux arm64 and Windows: find_cut runs without it."""

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(6)
        for q in range(6):
            qc.ry(0.3 + 0.1 * q, q)
        for q in range(5):
            qc.cz(q, q + 1)
        for q in range(6):
            qc.ry(0.2 * q, q)
        return qc

    def test_find_cut_works_says_so_and_the_cut_is_right(self):
        import sys
        from unittest import mock

        import Qdislib as qd
        from Qdislib.utils import logging_config

        qc = self._circuit()
        exact = qd.analytical_solution(qc, "Z" * 6)
        logging_config._WARNED.clear()
        with mock.patch.dict(sys.modules, {"pymetis": None}):
            with self.assertWarnsRegex(Warning, "pymetis"):
                cut = qd.find_cut(qc, max_qubits=4, wire_cut=False, seed=0)
            wire = qd.find_cut(qc, max_qubits=4, gate_cut=False, seed=0)
        self.assertTrue(cut)
        self.assertTrue(wire)
        value = qd.gate_cutting(qc, cut, "Z" * 6, shots=None)
        self.assertAlmostEqual(value, exact, places=8)
        value = qd.wire_cutting(qc, wire, "Z" * 6, shots=None)
        self.assertAlmostEqual(value, exact, places=8)

    def test_a_platform_without_a_wheel_does_not_require_it(self):
        import os

        from packaging.requirements import Requirement

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(root, "requirements.txt"), encoding="utf-8") as handle:
            lines = [l.strip() for l in handle if l.strip().startswith("pymetis")]
        marker = Requirement(lines[0]).marker
        env = {"sys_platform": "linux", "platform_machine": "x86_64"}
        self.assertTrue(marker.evaluate(env))
        env = {"sys_platform": "darwin", "platform_machine": "arm64"}
        self.assertTrue(marker.evaluate(env))
        env = {"sys_platform": "linux", "platform_machine": "aarch64"}
        self.assertFalse(marker.evaluate(env))
        env = {"sys_platform": "win32", "platform_machine": "AMD64"}
        self.assertFalse(marker.evaluate(env))


class FitsAlreadyTest(BaseTimedTestCase):
    """A connected circuit within max_qubits needs no cut."""

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.3, 1)
        qc.cz(1, 2)
        qc.ry(0.5, 2)
        return qc

    def test_no_cut_when_it_fits(self):
        import Qdislib as qd

        qc = self._circuit()
        for max_qubits in (3, 10):
            for kinds in ({}, {"gate_cut": False}, {"wire_cut": False}):
                with self.subTest(max_qubits=max_qubits, kinds=kinds):
                    self.assertEqual(qd.find_cut(qc, max_qubits, **kinds), [])

    def test_ranked_mode_says_so(self):
        import Qdislib as qd

        self.assertEqual(
            qd.find_cut(self._circuit(), max_qubits=3, top_n=3),
            [{"cut": [], "kind": "none", "score": 0.0}],
        )

    def test_estimate_runs_it_uncut(self):
        import Qdislib as qd

        qc = self._circuit()
        value, metrics = qd.estimate(
            qc, "ZZZ", max_qubits=3, shots=None, return_metrics=True
        )
        self.assertEqual(metrics["cut"], [])
        self.assertAlmostEqual(value, qd.analytical_solution(qc, "ZZZ"))

    def test_a_circuit_that_does_not_fit_is_still_cut(self):
        import Qdislib as qd

        self.assertTrue(qd.find_cut(self._circuit(), max_qubits=2, seed=0))

    def test_without_max_qubits_it_still_looks_for_one(self):
        import Qdislib as qd

        qc = self._circuit()
        self.assertTrue(qd.find_cut(qc, seed=0))
        self.assertTrue(qd.find_cut(qc, gate_cut=False, wire_cut=True, seed=0))

    def test_neither_kind_is_still_refused(self):
        import Qdislib as qd

        with self.assertRaises(ValueError):
            qd.find_cut(self._circuit(), 3, gate_cut=False, wire_cut=False)

    def test_routing_gates_are_still_cut_on_a_coupling_map(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(4)
        qc.cz(0, 1)
        qc.cz(2, 3)
        qc.cz(0, 3)
        line = [(0, 1), (1, 2), (2, 3)]
        self.assertEqual(qd.find_cut(qc, coupling_map=line, seed=0), ["CZ_3"])
        self.assertEqual(
            qd.find_cut(qc, max_qubits=2, coupling_map=line, seed=0), ["CZ_3"]
        )


class RankedCandidatesTest(BaseTimedTestCase):
    """find_cut(top_n=): each candidate once, in the ranking's order."""

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(6)
        for q in range(6):
            qc.h(q)
        for q in range(5):
            qc.cz(q, q + 1)
        for q in range(6):
            qc.ry(0.1 * q, q)
        qc.cz(0, 5)
        return qc

    def _ranked(self, **options):
        import Qdislib as qd

        return qd.find_cut(self._circuit(), max_qubits=3, seed=0, **options)

    def test_no_duplicates(self):
        for kinds in ({}, {"wire_cut": False}, {"gate_cut": False}):
            with self.subTest(kinds=kinds):
                ranked = self._ranked(top_n=100, **kinds)
                self.assertTrue(ranked)
                seen = [(c["kind"], tuple(sorted(map(str, c["cut"])))) for c in ranked]
                self.assertEqual(len(seen), len(set(seen)), ranked)

    def test_best_first_and_as_many_as_asked(self):
        ranked = self._ranked(top_n=100)
        scores = [c["score"] for c in ranked]
        self.assertEqual(scores, sorted(scores))
        self.assertEqual(self._ranked(top_n=2), ranked[:2])

    def test_the_winner_is_still_the_winner(self):
        ranked = self._ranked(top_n=1, wire_cut=False)
        self.assertEqual(ranked[0]["cut"], self._ranked(wire_cut=False))
