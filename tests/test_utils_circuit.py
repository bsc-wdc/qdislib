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

    def test_random_circuit(self):
        from Qdislib.utils.circuit import _random_circuit

        import qibo

        circuit = _random_circuit(5, 5, 3, None)

        if type(circuit) == qibo.models.Circuit:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_draw_to_circuit(self):
        from Qdislib.utils.circuit import draw_to_circuit

        circuit_draw = """
        q0: ─H─S─X───SDG───H───
        q1: ─X─S───o─X─T─o─S───
        q2: ───────Z─T───|───o─
        q3: ─────────────|─X─|─
        q4: ──RX─────RX──Z───Z─
        q0: ─H─S─X───SDG───H───
        q1: ─X─S───o─X─T─o─S───
        q2: ───────Z─T───|───o─
        q3: ─────────────|─X─|─
        q4: ──RX─────RX──Z───Z─
        """

        param = {4: [("RX",0.6),("RX",0.3)], 9: [("RX",0.2),("RX",0.8)]}

        circuit = draw_to_circuit(circuit_draw, param)

        # Do something to check that the draw_to_circuit function works
        self.assertTrue(True)

    def test_analytical_solution(self):
        from Qdislib.utils.circuit import analytical_solution

        import qibo
        from qibo import models, gates
        qibo.set_backend("numpy")

        circuit = models.Circuit(2)
        circuit.add(gates.H(0))
        circuit.add(gates.RY(1, 0.8))
        circuit.add(gates.H(1))
        circuit.add(gates.CZ(0,1))
        circuit.add(gates.H(1))
        circuit.add(gates.X(0))

        solution = analytical_solution(circuit,"ZZ")

        if abs(solution + 0.69) < 0.2:
            self.assertTrue(True)
        else:
            self.assertTrue(False)