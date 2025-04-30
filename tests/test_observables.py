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
from qibo import models, gates

def entire_circuit():
    circuit = models.Circuit(2)
    circuit.add(gates.H(0))
    circuit.add(gates.RY(1, 0.8))
    circuit.add(gates.H(1))
    circuit.add(gates.CZ(0,1))
    circuit.add(gates.H(1))
    circuit.add(gates.X(0))
    return circuit


class CircuitTest(BaseTimedTestCase):

    def test_observables_gate_cutting(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = entire_circuit()

        cut = ["CZ_4"]

        analytical = analytical_solution(circuit,"YY")
        circuit = entire_circuit()
        reconstruction = gate_cutting(circuit,cut,observables="YY")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_observables_gate_cutting2(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution
        
        circuit = entire_circuit()
        
        cut = ["CZ_4"]

        analytical = analytical_solution(circuit,"XX")
        circuit = entire_circuit()
        reconstruction = gate_cutting(circuit,cut,observables="XX")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_observables_gate_cutting3(self):
        from Qdislib.core.cutting_algorithms.gate_cutting import gate_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = entire_circuit()

        cut = ["CZ_4"]

        analytical = analytical_solution(circuit,"YX")
        circuit = entire_circuit()
        reconstruction = gate_cutting(circuit,cut,observables="YX")

        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)


    def test_observables_wire_cutting(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = entire_circuit()

        cut = [("H_3","CZ_4")]

        analytical = analytical_solution(circuit,"YY")
        circuit = entire_circuit()
        reconstruction = wire_cutting(circuit,cut,observables="YY")
        
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)


    def test_observables_wire_cutting2(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = entire_circuit()

        cut = [("H_3","CZ_4")]

        analytical = analytical_solution(circuit,"XX")
        circuit = entire_circuit()
        reconstruction = wire_cutting(circuit,cut,observables="XX")
        
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
    
    def test_observables_wire_cutting3(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import wire_cutting
        from Qdislib.utils.circuit import analytical_solution

        circuit = entire_circuit()

        cut = [("H_3","CZ_4")]

        analytical = analytical_solution(circuit,"XY")

        circuit = entire_circuit()
        reconstruction = wire_cutting(circuit,cut,observables="XY")
        
        if abs(reconstruction - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)


