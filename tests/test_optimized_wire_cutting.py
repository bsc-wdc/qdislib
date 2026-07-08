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
        """Two cuts split into one measure + one prepare fragment: 3**2 + 4**2 = 25."""
        import Qdislib.core.cutting_algorithms.wire_cutting as wc

        calls = {"n": 0}
        original = wc._sample_wire_fragment

        def counting(*args, **kwargs):
            calls["n"] += 1
            return original(*args, **kwargs)

        wc._sample_wire_fragment = counting
        try:
            wc.wire_cutting(_example_circuit(), _TWO_CUTS, shots=4000, optimized=True)
        finally:
            wc._sample_wire_fragment = original

        self.assertEqual(calls["n"], 3**2 + 4**2)

    def test_falls_back_to_standard_when_not_separable(self):
        """A topology the optimizer can't separate returns the standard result."""
        import Qdislib.core.cutting_algorithms.wire_cutting as wc
        from Qdislib.utils.circuit import analytical_solution

        circuit = _example_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)

        # Force the "not separable" path; wire_cutting must still reconstruct.
        original = wc._wire_cutting_optimized
        wc._wire_cutting_optimized = lambda *args, **kwargs: None
        try:
            value = wc.wire_cutting(circuit, _TWO_CUTS, shots=10000, optimized=True)
        finally:
            wc._wire_cutting_optimized = original

        self.assertLess(abs(value - analytical), 0.2)

    def test_three_step_matches_standard_and_analytical(self):
        import Qdislib as qd
        from Qdislib.utils.circuit import analytical_solution

        circuit = _example_circuit()
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        subs = qd.wire_cutting_subcircuits(
            circuit, _TWO_CUTS, software="qibo", optimized=True
        )
        # Two cuts -> one measure + one prepare fragment: 3**2 + 4**2 = 25 subcircuits
        # (a Z-basis measure probe yields both <Z> and <I>, so the measure side is 3**n).
        self.assertEqual(len(subs), 3**2 + 4**2)
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
            circuit, _TWO_CUTS, shots=60000, software="qibo", seed=1, optimized=True
        )
        tvd = 0.5 * sum(
            abs(distribution.get(key, 0.0) - exact.get(key, 0.0))
            for key in set(distribution) | set(exact)
        )
        self.assertLess(tvd, 0.08)


if __name__ == "__main__":
    unittest.main()
