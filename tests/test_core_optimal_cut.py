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

    def test_optimal_cut_no_cut(self):
        from Qdislib.core.find_cut.find_cut import find_cut
        from Qdislib.utils.circuit import analytical_solution

        from qibo import models, gates
        import numpy as np

        import qibo
        qibo.set_backend("numpy")
        
        circ = models.Circuit(2)
        circ.add(gates.H(0))
        circ.add(gates.CZ(0,1))
        circ.add(gates.RY(0,0.8))

        cut = find_cut(circ, max_qubits=2, gate_cut=False)
        if cut == None:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

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
        circuit.add(gates.CZ(0,1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1,2))
        circuit.add(gates.RY(1,np.pi/5))
        circuit.add(gates.RZ(2,np.pi/3))
        circuit.add(gates.RZ(3,np.pi/2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2,3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 3, gate_cut=False)

        reconstrucion = wire_cutting(circuit, cut)
        analytical = analytical_solution(circuit,"ZZZZ")

        if cut == [('CZ_5', 'RZ_7')] and abs(reconstrucion - analytical) < 0.2:
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
        circuit.add(gates.CZ(0,1))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(1,2))
        circuit.add(gates.RY(1,np.pi/5))
        circuit.add(gates.RZ(2,np.pi/3))
        circuit.add(gates.RZ(3,np.pi/2))
        circuit.add(gates.H(3))
        circuit.add(gates.CZ(2,3))
        circuit.add(gates.H(3))

        cut = find_cut(circuit, 2, wire_cut=False)

        reconstrucion = gate_cutting(dag=circuit, gates_cut=cut)
        analytical = analytical_solution(circuit,"ZZZZ")

        if cut == ['CZ_5'] and abs(reconstrucion - analytical) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)