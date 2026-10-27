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

"""Optimized (fragment-separated, 3**n + 4**n) wire cutting."""

import importlib.util
import unittest

import numpy as np

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None


def _example_circuit():
    """A 4-qubit circuit cuttable at two wires that both feed the same CZ."""
    from qibo import models, gates

    circuit = models.Circuit(4)
    circuit.add(gates.H(0))
    circuit.add(gates.H(1))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.CZ(1, 2))
    circuit.add(gates.H(1))
    circuit.add(gates.RZ(2, np.pi / 3))
    circuit.add(gates.CZ(1, 2))
    circuit.add(gates.RY(1, np.pi / 5))
    circuit.add(gates.RZ(3, np.pi / 2))
    circuit.add(gates.H(3))
    circuit.add(gates.CZ(2, 3))
    circuit.add(gates.H(3))
    return circuit


# Two wire cuts that both feed CZ_7 -> one measure fragment + one prepare fragment.
_TWO_CUTS = [("H_5", "CZ_7"), ("RZ_6", "CZ_7")]


def _chain_circuit():
    """A 3-qubit linear chain, cuttable at two DIFFERENT wires in series.

    q0 --CZ_2-- q1 --CZ_4-- q2. Cutting both wires puts the middle fragment
    (CZ_2's q1 side + RY_3) in the position of owning both a prepare marker
    (from cut 1) and a measure marker (from cut 2) -- a genuine chain/mixed
    topology, unlike ``_example_circuit``'s two cuts converging on one gate.
    """
    from qibo import models, gates

    circuit = models.Circuit(3)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.RY(1, np.pi / 7))
    circuit.add(gates.CZ(1, 2))
    circuit.add(gates.RZ(2, np.pi / 5))
    return circuit


# Two DIFFERENT wires cut in series -> the middle fragment (CZ_2/RY_3) owns
# both a prepare marker (cut 1) and a measure marker (cut 2): a mixed fragment.
_CHAIN_CUTS = [("H_1", "CZ_2"), ("RY_3", "CZ_4")]


def _fragment_runs(call):
    """The fragment circuits ``call()`` runs: one engine task each
    (:func:`~Qdislib.core.cutting_algorithms.mixed_cutting._run_fragment`)."""
    import Qdislib.core.cutting_algorithms.mixed_cutting as mc

    original = mc._run_fragment
    runs = []

    def counting(*args, **kwargs):
        runs.append(args)
        return original(*args, **kwargs)

    mc._run_fragment = counting
    try:
        call()
    finally:
        mc._run_fragment = original
    return runs


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class OptimizedWireCuttingTest(BaseTimedTestCase):
    def test_two_cuts_match_standard_and_analytical(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = _example_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        standard = wire_cutting(circuit, _TWO_CUTS, shots=10000, optimized=False)
        optimized = wire_cutting(circuit, _TWO_CUTS, shots=10000, optimized=True)

        self.assertLess(abs(optimized - analytical), 0.2)
        self.assertLess(abs(optimized - standard), 0.2)

    def test_one_cut_matches_analytical(self):
        from qibo import models, gates
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.utils.circuit import analytical_solution

        circuit = models.Circuit(4)
        circuit.add(gates.RY(0, np.pi / 6))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.RY(1, np.pi / 7))
        circuit.add(gates.CZ(1, 2))
        circuit.add(gates.RY(2, np.pi / 6))
        circuit.add(gates.CZ(2, 3))
        cut = find_cut(circuit, gate_cut=False)[:1]

        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        optimized = wire_cutting(circuit, cut, shots=10000, optimized=True)
        self.assertLess(abs(optimized - analytical), 0.2)

    def test_uses_3n_plus_4n_subcircuits(self):
        """Two cuts split into one measure + one prepare fragment: 3**2 + 4**2 = 25
        measured per wire, 5 + 4**2 = 21 measured jointly (the default)."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

        for joint, expected in ((False, 3**2 + 4**2), (True, 2**2 + 1 + 4**2)):
            runs = _fragment_runs(
                lambda: wire_cutting(
                    _example_circuit(),
                    _TWO_CUTS,
                    shots=4000,
                    optimized=True,
                    joint=joint,
                )
            )
            self.assertEqual(len(runs), expected, joint)

    def test_chain_topology_matches_standard_and_analytical(self):
        """Two cuts in series (a mixed fragment) still reconstruct correctly."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = _chain_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        standard = wire_cutting(circuit, _CHAIN_CUTS, shots=10000, optimized=False)
        optimized = wire_cutting(circuit, _CHAIN_CUTS, shots=10000, optimized=True)

        self.assertLess(abs(optimized - analytical), 0.2)
        self.assertLess(abs(optimized - standard), 0.2)

    def test_chain_topology_non_trivial_observable_matches_analytical(self):
        """A non-degenerate expectation value (not -1/0/1, unlike a pure
        GHZ/stabilizer circuit's Pauli expectations) with an X component
        directly on a cut qubit -- exercises the basis-rotation gate
        change_basis_dag inserts, which needs the same downstream relabeling
        as an Observable-I marker (see test_core_wire_cutting.py's
        identity-observable regression tests)."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = _chain_circuit()
        obs = "XXZ"  # X on both cut qubits (1 and 2 in _CHAIN_CUTS)
        analytical = analytical_solution(circuit, obs)
        self.assertGreater(abs(analytical), 0.05)  # sanity: not a degenerate ~0 case

        standard = wire_cutting(
            circuit, _CHAIN_CUTS, observable=obs, shots=20000, optimized=False
        )
        optimized = wire_cutting(
            circuit, _CHAIN_CUTS, observable=obs, shots=20000, optimized=True
        )

        self.assertLess(abs(standard - analytical), 0.1)
        self.assertLess(abs(optimized - analytical), 0.1)

    def test_optimized_runs_through_circuit_cutting(self):
        """optimized=True is the engine's optimized plan, not a separate path."""
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc

        original = mc._prepare_cut
        optimized = []

        def wrapped(*args, **kwargs):
            optimized.append(
                args[4]
            )  # (circuit, cuts, observable, hadamard, optimized)
            return original(*args, **kwargs)

        mc._prepare_cut = wrapped
        try:
            mc.circuit_cutting  # noqa: B018 -- the module the wire cut routes to
            from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

            wire_cutting(_chain_circuit(), _CHAIN_CUTS, shots=2000, optimized=True)
        finally:
            mc._prepare_cut = original
        self.assertEqual(optimized, [True])

    def test_weighted_sum_runs_each_term_optimized(self):
        """A Pauli-sum observable takes the optimized path for every term."""
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = _example_circuit()
        hamiltonian = [(0.5, "ZZZZ"), (-1.0, "ZIZI")]
        exact = sum(
            coeff * analytical_solution(circuit, pauli) for coeff, pauli in hamiltonian
        )

        original = mc._prepare_cut
        calls = []

        def wrapped(*args, **kwargs):
            calls.append(args[4])
            return original(*args, **kwargs)

        mc._prepare_cut = wrapped
        try:
            value = wire_cutting(
                circuit,
                _TWO_CUTS,
                observable=hamiltonian,
                shots=10000,
                optimized=True,
            )
        finally:
            mc._prepare_cut = original

        self.assertEqual(calls, [True] * len(hamiltonian))
        self.assertLess(abs(value - exact), 0.2)

    def test_chain_topology_uses_fewer_than_8n_subcircuits(self):
        """Measure-only + mixed + prepare-only fragments, not 8**2 = 64 circuits:
        3 + 3*4 + 4 = 19 in general. This chain is real (H, CZ, RY), so both wires
        are golden cut points -- their Y terms are zero and left out -- and it is
        2 + 2*3 + 3 = 11."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

        runs = _fragment_runs(
            lambda: wire_cutting(
                _chain_circuit(), _CHAIN_CUTS, shots=2000, optimized=True
            )
        )
        self.assertEqual(len(runs), 2 + 2 * 3 + 3)

    def test_three_step_matches_standard_and_analytical(self):
        import Qdislib as qd
        from Qdislib.utils.circuit import analytical_solution

        circuit = _example_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        # Two cuts -> one measure + one prepare fragment: 3**2 + 4**2 = 25 subcircuits
        # measured per wire (a Z-basis probe yields both <Z> and <I>), 5 + 16 = 21
        # measured jointly (the default).
        for joint, expected in ((False, 3**2 + 4**2), (True, 2**2 + 1 + 4**2)):
            subs = qd.wire_cutting_subcircuits(
                circuit, _TWO_CUTS, software="qibo", optimized=True, joint=joint
            )
            self.assertEqual(len(subs), expected)
            values = [qd.expectation_value(s, shots=10000) for s in subs]
            self.assertLess(abs(subs.reconstruct(values) - analytical), 0.2)

    def test_three_step_chain_topology_matches_standard_and_analytical(self):
        """The three-step API also handles the mixed-fragment chain topology."""
        import Qdislib as qd
        from Qdislib.utils.circuit import analytical_solution

        circuit = _chain_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        subs = qd.wire_cutting_subcircuits(
            circuit, _CHAIN_CUTS, software="qibo", optimized=True
        )
        self.assertIsNotNone(subs)
        # Measure-only + mixed + prepare-only, not 8**2 = 64: 3 + 3*4 + 4 = 19 in
        # general, 2 + 2*3 + 3 = 11 here (both wires are golden cut points), the
        # same circuits as the all-in-one call.
        self.assertEqual(len(subs), 2 + 2 * 3 + 3)
        values = [qd.expectation_value(s, shots=10000) for s in subs]
        self.assertLess(abs(subs.reconstruct(values) - analytical), 0.2)

    def test_sampling_matches_exact_distribution(self):
        import numpy as np
        import Qdislib as qd

        circuit = _example_circuit()
        n = circuit.nqubits
        probs = np.abs(np.asarray(circuit().state(numpy=True))) ** 2
        exact = {
            "".join(str((x >> (n - 1 - q)) & 1) for q in range(n)): float(probs[x])
            for x in range(2**n)
        }
        distribution = qd.wire_cutting_sampling(
            circuit, _TWO_CUTS, shots=60000, backend="qibo", seed=1, optimized=True
        )
        tvd = 0.5 * sum(
            abs(distribution.get(key, 0.0) - exact.get(key, 0.0))
            for key in set(distribution) | set(exact)
        )
        self.assertLess(tvd, 0.08)

    def test_sampling_chain_topology_matches_exact_distribution(self):
        """The sampling API also handles the mixed-fragment chain topology."""
        import numpy as np
        import Qdislib as qd

        circuit = _chain_circuit()
        n = circuit.nqubits
        probs = np.abs(np.asarray(circuit().state(numpy=True))) ** 2
        exact = {
            "".join(str((x >> (n - 1 - q)) & 1) for q in range(n)): float(probs[x])
            for x in range(2**n)
        }
        distribution = qd.wire_cutting_sampling(
            circuit, _CHAIN_CUTS, shots=60000, backend="qibo", seed=1, optimized=True
        )
        tvd = 0.5 * sum(
            abs(distribution.get(key, 0.0) - exact.get(key, 0.0))
            for key in set(distribution) | set(exact)
        )
        self.assertLess(tvd, 0.08)

    def test_sampling_chain_topology_does_not_fall_back(self):
        """The sampling path handles the chain directly, not via the fallback."""
        import Qdislib.core.cutting_algorithms.wire_cutting as wc

        original = wc._wire_sampling_optimized
        called = {"returned_none": False}

        def wrapped(*args, **kwargs):
            result = original(*args, **kwargs)
            called["returned_none"] = result is None
            return result

        wc._wire_sampling_optimized = wrapped
        try:
            wc.wire_cutting_sampling(
                _chain_circuit(),
                _CHAIN_CUTS,
                shots=2000,
                backend="qibo",
                seed=1,
                optimized=True,
            )
        finally:
            wc._wire_sampling_optimized = original

        self.assertFalse(called["returned_none"])


if __name__ == "__main__":
    unittest.main()
