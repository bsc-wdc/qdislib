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


class CircuitTest(BaseTimedTestCase):

    def test_gate_cutting(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
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
        circ.add(gates.CZ(0,1))
        circ.add(gates.RY(0,0.8))

        reconstruction = gate_cutting(qc, ["CZ_2"])
        analytical = analytical_solution(circ,"ZZ")
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_gate_cutting_more_components(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo
        qibo.set_backend("numpy")
        
        circuit = models.Circuit(4)
        circuit.add(gates.H(0))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0,1))
        circuit.add(gates.H(1))
        circuit.add(gates.RY(1,np.pi/5))
        circuit.add(gates.RZ(2,np.pi/3))
        circuit.add(gates.RZ(3,np.pi/2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2,3))
        circuit.add(gates.H(3))

        reconstruction = gate_cutting(circuit, ["CZ_3"])

        analytical = analytical_solution(circuit, "ZZZZ")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)