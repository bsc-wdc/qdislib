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

"""What a CUDA-Q kernel's source says: its registers, measurements and adjoints."""

import importlib.util
import unittest
import warnings

from tests import BaseTimedTestCase

CUDAQ_AVAILABLE = importlib.util.find_spec("cudaq") is not None

if CUDAQ_AVAILABLE:
    import cudaq

    # Kernels are compiled when defined, from their own source: so here, one by one.

    @cudaq.kernel
    def two_registers():
        a = cudaq.qvector(1)
        b = cudaq.qvector(1)
        x(b[0])

    @cudaq.kernel
    def wider_registers():
        a = cudaq.qvector(2)
        b = cudaq.qvector(1)
        c = cudaq.qubit()
        x(b[0])
        x.ctrl(b[0], a[1])

    @cudaq.kernel
    def outside_its_register():
        a = cudaq.qvector(1)
        b = cudaq.qvector(2)
        x(a[1])


@unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq not installed")
class CudaqRegistersTest(BaseTimedTestCase):
    def test_qubits_of_each_register_are_told_apart(self):
        import Qdislib as qd
        from Qdislib.utils.backends import convert_to_dag

        dag = convert_to_dag(two_registers)
        self.assertEqual([data["qubits"] for _, data in dag.nodes(data=True)], [(1,)])
        self.assertAlmostEqual(qd.gate_cutting(two_registers, [], "ZI", shots=None), 1)
        self.assertAlmostEqual(qd.gate_cutting(two_registers, [], "IZ", shots=None), -1)
        self.assertEqual(
            qd.gate_cutting(two_registers, [], "IZ", shots=100, seed=1), -1.0
        )
        self.assertEqual(
            qd.gate_cutting(two_registers, [], "IZ", shots=None),
            qd.analytical_solution(two_registers, "IZ"),
        )

    def test_a_register_starts_after_the_ones_before_it(self):
        from Qdislib.utils.backends import convert_to_dag

        dag = convert_to_dag(wider_registers)
        self.assertEqual(
            [data["qubits"] for _, data in dag.nodes(data=True)], [(2,), (2, 1)]
        )

    def test_a_qubit_outside_its_register_is_refused(self):
        from Qdislib.utils.backends import convert_to_dag

        with self.assertRaisesRegex(ValueError, "Cannot read this CUDA-Q kernel"):
            convert_to_dag(outside_its_register)


if __name__ == "__main__":
    unittest.main()
