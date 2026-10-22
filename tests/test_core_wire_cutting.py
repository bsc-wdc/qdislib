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

import importlib.util
import unittest

from tests import BaseTimedTestCase

CUDAQ_AVAILABLE = importlib.util.find_spec("cudaq") is not None


class CircuitTest(BaseTimedTestCase):

    def test_wire_cutting(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cz(0, 1)
        qc.ry(0.8, 0)

        from qibo import models, gates
        import numpy as np

        import qibo

        qibo.set_backend("numpy")

        circ = models.Circuit(2)
        circ.add(gates.H(0))
        circ.add(gates.CZ(0, 1))
        circ.add(gates.RY(0, 0.8))

        reconstruction = wire_cutting(qc, [("H_1", "CZ_2")])
        analytical = analytical_solution(circ, "ZZ")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_wire_cutting_more_components(self):
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
        circuit.add(gates.RY(1, np.pi / 5))
        circuit.add(gates.RZ(2, np.pi / 3))
        circuit.add(gates.RZ(3, np.pi / 2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2, 3))
        circuit.add(gates.H(3))

        print(circuit.draw())

        reconstruction = wire_cutting(circuit, [("CZ_3", "H_4")])

        analytical = analytical_solution(circuit, "ZZZZ")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_wire_cutting_qibo(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution
        from Qdislib.core.cutting_algorithms.wire_cutting import _expec_value_qibo
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            wire_cutting_subcircuits,
            wire_cutting_subcircuit_reconstruction,
        )
        from Qdislib.core.find_cut.find_cut import find_cut

        from qibo import models, gates
        import numpy as np

        import qibo

        nqubits = 10
        circuit = models.Circuit(nqubits)

        circuit.add(gates.H(0))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 6))
        circuit.add(gates.RZ(8, np.pi / 3))

        circuit.add(gates.RY(3, np.pi / 5))
        circuit.add(gates.RX(4, np.pi / 5))
        circuit.add(gates.CZ(0, 2))
        circuit.add(gates.CZ(5, 9))

        circuit.add(gates.CZ(3, 5))
        circuit.add(gates.CZ(3, 4))
        circuit.add(gates.CZ(6, 7))
        circuit.add(gates.RY(7, np.pi / 5))
        circuit.add(gates.RZ(1, np.pi / 5))

        circuit.add(gates.CZ(1, 5))
        circuit.add(gates.RX(6, np.pi / 5))
        circuit.add(gates.CZ(7, 8))

        circuit.add(gates.H(9))

        cut = find_cut(circuit, gate_cut=False)
        print(cut)
        subcircuits = wire_cutting_subcircuits(circuit, cut, software="qibo")

        results = []
        for subcircuit in subcircuits:
            print(subcircuit)
            # Execute individually each subcircuit
            result = _expec_value_qibo(subcircuit)
            results.append(result)

        # results = compss_wait_on(results)

        reconstruction = wire_cutting_subcircuit_reconstruction(results, number_cuts=1)

        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    @staticmethod
    def _two_cuts_circuit():
        """4-qubit circuit cuttable at two wires that both feed the same CZ,
        needing a full downstream chain (H_5/RZ_6/CZ_7/RY_8) on the cut qubits
        plus qubit3's own CZ_11 -- the shape that exposed the identity-
        observable regressions below (a shorter/simpler circuit doesn't)."""
        from qibo import models, gates
        import numpy as np

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

    def test_wire_cutting_identity_observable_on_regular_qubit(self):
        """Regression: evaluate_samples() used to drop the Observable-I marker
        list dag_to_circuit_* returns, so the optimized (raw-sampling) wire-cut
        path folded an identity-observed *regular* qubit into the Z-parity as
        if it were really observed with Z. Fixed by
        Backend.marginalize_identity_qubits; the standard (expectation-value)
        path was never affected, so this mainly guards the optimized path."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = self._two_cuts_circuit()
        cuts = [("H_5", "CZ_7"), ("RZ_6", "CZ_7")]
        obs = "IZZZ"  # qubit 0: not a cut qubit
        analytical = analytical_solution(circuit, obs)
        self.assertLess(abs(analytical), 0.05)  # sanity: this IS a degenerate ~0 case

        standard = wire_cutting(
            circuit, cuts, observable=obs, shots=20000, optimized=False
        )
        optimized = wire_cutting(
            circuit, cuts, observable=obs, shots=20000, optimized=True
        )
        self.assertLess(abs(standard - analytical), 0.15)
        self.assertLess(abs(optimized - analytical), 0.15)

    def test_wire_cutting_identity_observable_on_cut_qubit(self):
        """Regression: the wire-cut marker relabeling that moves the "prepared"
        side of a cut onto a fresh qubit line only walked one hop of successors,
        so a node further downstream on that qubit -- e.g. an Observable-I
        terminal marker change_basis_dag inserted after more gates on the same
        line -- kept its stale, pre-cut qubit index. It ended up correctly
        placed in the right fragment but pointing at the wrong (or a phantom)
        qubit position within it, corrupting the reconstruction. Fixed by
        making _relabel_qubit_downstream walk the whole downstream chain.
        Affects both the standard and optimized paths (same marker-insertion
        code); only reproduces with a real downstream chain past the cut
        target, which is why this needs _two_cuts_circuit rather than a
        minimal one."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = self._two_cuts_circuit()
        cuts = [("H_5", "CZ_7"), ("RZ_6", "CZ_7")]
        obs = "ZIZZ"  # qubit 1: IS a cut qubit
        analytical = analytical_solution(circuit, obs)
        self.assertLess(abs(analytical), 0.05)  # sanity: this IS a degenerate ~0 case

        standard = wire_cutting(
            circuit, cuts, observable=obs, shots=20000, optimized=False
        )
        optimized = wire_cutting(
            circuit, cuts, observable=obs, shots=20000, optimized=True
        )
        self.assertLess(abs(standard - analytical), 0.15)
        self.assertLess(abs(optimized - analytical), 0.15)

    def test_wire_cutting_non_partitioning_cut_raises(self):
        """A cut whose wire isn't the sole connectivity between its two sides
        (here: CZ_1 and CZ_2 also stay linked through CZ_3, via the OTHER two
        qubits) doesn't separate the circuit at all -- wire_cutting should
        refuse it rather than silently pay the cut's shot-count cost for none
        of its benefit."""
        from qibo import models, gates
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting

        circuit = models.Circuit(3)
        circuit.add(gates.CZ(0, 1))  # CZ_1
        circuit.add(gates.CZ(1, 2))  # CZ_2
        circuit.add(gates.CZ(0, 2))  # CZ_3 -- keeps CZ_1/CZ_2 connected regardless

        with self.assertRaises(ValueError):
            wire_cutting(circuit, [("CZ_1", "CZ_2")], shots=1000, optimized=False)

        with self.assertRaises(ValueError):
            wire_cutting(circuit, [("CZ_1", "CZ_2")], shots=1000, optimized=True)

    def test_wire_cutting_empty_cut_runs_unpartitioned(self):
        """An empty cut list is exempt from the partition check: it means 'no
        cut,' so the circuit just runs as-is."""
        from qibo import models, gates
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = models.Circuit(2)
        circuit.add(gates.H(0))
        circuit.add(gates.CZ(0, 1))

        result = wire_cutting(circuit, [], shots=1024)
        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        self.assertLess(abs(result - analytical), 0.2)

    def test_wire_cutting_qiskit(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution
        from Qdislib.core.cutting_algorithms.gate_cutting import _expec_value_qiskit
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            wire_cutting_subcircuits,
            wire_cutting_subcircuit_reconstruction,
        )
        from Qdislib.core.find_cut.find_cut import find_cut

        from qibo import models, gates
        import numpy as np

        import qibo

        nqubits = 10
        circuit = models.Circuit(nqubits)

        circuit.add(gates.H(0))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 6))
        circuit.add(gates.RZ(8, np.pi / 3))

        circuit.add(gates.RY(3, np.pi / 5))
        circuit.add(gates.RX(4, np.pi / 5))
        circuit.add(gates.CZ(0, 2))
        circuit.add(gates.CZ(5, 9))

        circuit.add(gates.CZ(3, 5))
        circuit.add(gates.CZ(3, 4))
        circuit.add(gates.CZ(6, 7))
        circuit.add(gates.RY(7, np.pi / 5))
        circuit.add(gates.RZ(1, np.pi / 5))

        circuit.add(gates.CZ(1, 5))
        circuit.add(gates.RX(6, np.pi / 5))
        circuit.add(gates.CZ(7, 8))

        circuit.add(gates.H(9))

        cut = find_cut(circuit, gate_cut=False)
        print(cut)
        subcircuits = wire_cutting_subcircuits(circuit, cut, software="qiskit")

        results = []
        for subcircuit in subcircuits:
            # Execute individually each subcircuit
            result = _expec_value_qiskit(subcircuit)
            results.append(result)

        reconstruction = wire_cutting_subcircuit_reconstruction(results, number_cuts=1)

        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_wire_cutting_default_software_follows_circuit(self):
        """bsc-wdc/qdislib#14: a Qiskit circuit with cx runs without naming software."""
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution
        from qibo import models, gates
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        qc.ry(0.3, 1)

        circ = models.Circuit(2)
        circ.add(gates.H(0))
        circ.add(gates.CNOT(0, 1))
        circ.add(gates.RY(1, 0.3))

        reconstruction = wire_cutting(
            qc, [("CX_2", "RY_3")], observable="ZZ", shots=20000
        )
        self.assertLess(abs(reconstruction - analytical_solution(circ, "ZZ")), 0.05)

    @unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq not installed")
    def test_wire_cutting_cudaq(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution
        from Qdislib.core.cutting_algorithms.wire_cutting import _expec_value_cudaq
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            wire_cutting_subcircuits,
            wire_cutting_subcircuit_reconstruction,
        )
        from Qdislib.utils.runtime import compss_wait_on
        from Qdislib.core.find_cut.find_cut import find_cut

        from qibo import models, gates
        import numpy as np

        import qibo

        nqubits = 10
        circuit = models.Circuit(nqubits)

        circuit.add(gates.H(0))
        circuit.add(gates.CZ(0, 1))
        circuit.add(gates.CZ(2, 6))
        circuit.add(gates.RZ(8, np.pi / 3))

        circuit.add(gates.RY(3, np.pi / 5))
        circuit.add(gates.RX(4, np.pi / 5))
        circuit.add(gates.CZ(0, 2))
        circuit.add(gates.CZ(5, 9))

        circuit.add(gates.CZ(3, 5))
        circuit.add(gates.CZ(3, 4))
        circuit.add(gates.CZ(6, 7))
        circuit.add(gates.RY(7, np.pi / 5))
        circuit.add(gates.RZ(1, np.pi / 5))

        circuit.add(gates.CZ(1, 5))
        circuit.add(gates.RX(6, np.pi / 5))
        circuit.add(gates.CZ(7, 8))

        circuit.add(gates.H(9))

        cut = find_cut(circuit, gate_cut=False)
        print(cut)
        subcircuits = wire_cutting_subcircuits(circuit, cut, software="cudaq")

        # subcircuits = compss_wait_on(subcircuits)

        results = []
        for subcircuit in subcircuits:
            # Execute individually each subcircuit
            result = _expec_value_cudaq(subcircuit)
            results.append(result)

        results = compss_wait_on(results)
        print(results)

        reconstruction = wire_cutting_subcircuit_reconstruction(results, number_cuts=1)

        analytical = analytical_solution(circuit, "Z" * circuit.nqubits)
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
