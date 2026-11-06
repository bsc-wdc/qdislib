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

    @cudaq.kernel
    def plain():
        q = cudaq.qvector(3)
        ry(0.4, q[0])
        x.ctrl(q[0], q[1])
        ry(0.3, q[1])
        z.ctrl(q[1], q[2])
        ry(0.5, q[2])

    @cudaq.kernel
    def final_each():
        q = cudaq.qvector(3)
        ry(0.4, q[0])
        x.ctrl(q[0], q[1])
        ry(0.3, q[1])
        z.ctrl(q[1], q[2])
        ry(0.5, q[2])
        mz(q[0])
        mz(q[1])
        mz(q[2])

    @cudaq.kernel
    def final_all():
        q = cudaq.qvector(3)
        ry(0.4, q[0])
        x.ctrl(q[0], q[1])
        ry(0.3, q[1])
        z.ctrl(q[1], q[2])
        ry(0.5, q[2])
        mz(q)

    @cudaq.kernel
    def final_assigned():
        q = cudaq.qvector(3)
        ry(0.4, q[0])
        x.ctrl(q[0], q[1])
        ry(0.3, q[1])
        z.ctrl(q[1], q[2])
        ry(0.5, q[2])
        m = mz(q)

    @cudaq.kernel
    def final_part():
        q = cudaq.qvector(2)
        x(q[1])
        mz(q[1])

    @cudaq.kernel
    def measured_before_another_qubit():
        q = cudaq.qvector(2)
        x(q[0])
        mz(q[0])
        x(q[1])

    @cudaq.kernel
    def final_two_registers():
        a = cudaq.qvector(1)
        b = cudaq.qvector(1)
        x(b[0])
        mz(a)
        mz(b)

    @cudaq.kernel
    def mid_circuit():
        q = cudaq.qvector(2)
        h(q[0])
        m = mz(q[0])
        h(q[0])
        x.ctrl(q[0], q[1])

    @cudaq.kernel
    def gate_after_whole_register():
        q = cudaq.qvector(2)
        h(q[0])
        mz(q)
        h(q[1])

    @cudaq.kernel
    def x_basis():
        q = cudaq.qvector(1)
        h(q[0])
        mx(q[0])

    @cudaq.kernel
    def adjoint_t():
        q = cudaq.qvector(1)
        h(q[0])
        t.adj(q[0])
        h(q[0])

    @cudaq.kernel
    def adjoint_s():
        q = cudaq.qvector(1)
        h(q[0])
        s.adj(q[0])
        h(q[0])

    @cudaq.kernel
    def adjoints_cut():
        q = cudaq.qvector(2)
        h(q[0])
        t.adj(q[0])
        h(q[0])
        s.adj(q[1])
        h(q[1])
        z.ctrl(q[0], q[1])
        h(q[1])

    @cudaq.kernel
    def adjoints_of_paulis():
        q = cudaq.qvector(1)
        x.adj(q[0])
        h.adj(q[0])
        z.adj(q[0])
        y.adj(q[0])
        h(q[0])


@unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq not installed")
class CudaqAdjointsTest(BaseTimedTestCase):
    def test_every_adjoint_a_kernel_can_write(self):
        import Qdislib as qd

        for kernel, observable, cut in (
            (adjoint_t, "Z", []),
            (adjoint_s, "Z", []),
            (adjoints_cut, "ZZ", []),
            (adjoints_cut, "ZZ", ["CZ_6"]),
        ):
            with self.subTest(kernel=kernel.name, cut=cut):
                exact = qd.analytical_solution(kernel, observable)
                self.assertAlmostEqual(
                    qd.gate_cutting(kernel, cut, observable, shots=None), exact
                )
                self.assertAlmostEqual(
                    qd.gate_cutting(kernel, cut, observable, shots=40000, seed=1),
                    exact,
                    delta=0.05,
                )

    def test_the_self_adjoint_ones(self):
        import Qdislib as qd

        # read from the source: CUDA-Q itself does not compile h.adj in a kernel
        self.assertAlmostEqual(
            qd.gate_cutting(adjoints_of_paulis, [], "Z", shots=None), -1.0
        )
        self.assertEqual(
            qd.gate_cutting(adjoints_of_paulis, [], "Z", shots=100, seed=1), -1.0
        )

    def test_they_build_as_the_adjoint(self):
        import numpy as np

        from Qdislib.utils.backends import convert_to_dag
        from Qdislib.utils.graph_cudaq import dag_to_circuit_cudaq_kernel

        kernel, _ = dag_to_circuit_cudaq_kernel(convert_to_dag(adjoint_t), 1)
        state = np.asarray(cudaq.get_state(kernel), dtype=complex)
        want = np.asarray(cudaq.get_state(adjoint_t), dtype=complex)
        self.assertTrue(np.allclose(state, want))
        self.assertAlmostEqual(abs(state[0]) ** 2 - abs(state[1]) ** 2, 2**-0.5)


@unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq not installed")
class CudaqMeasurementsTest(BaseTimedTestCase):
    def test_final_measurements_do_not_change_the_answer(self):
        import Qdislib as qd

        exact = qd.analytical_solution(plain, "ZZZ")
        self.assertAlmostEqual(exact, 0.83838664, places=6)
        for kernel in (final_each, final_all, final_assigned):
            with self.subTest(kernel=kernel.name):
                for _ in range(8):  # get_state collapsed: +-1 at random
                    self.assertAlmostEqual(
                        qd.analytical_solution(kernel, "ZZZ"), exact, places=6
                    )
                self.assertAlmostEqual(
                    qd.gate_cutting(kernel, ["CZ_4"], "ZZZ", shots=None),
                    exact,
                    places=6,
                )
                for run in (
                    lambda: qd.gate_cutting(
                        kernel, ["CZ_4"], "ZZZ", shots=8000, seed=1
                    ),
                    lambda: qd.wire_cutting(
                        kernel, [("RY_3", "CZ_4")], "ZZZ", shots=8000, seed=1
                    ),
                    lambda: qd.gate_cutting_hadamard(
                        kernel, ["CZ_4"], "ZZZ", shots=8000, seed=1
                    ),
                ):
                    self.assertAlmostEqual(run(), exact, delta=0.15)

    def test_they_are_not_a_gate_to_cut(self):
        from Qdislib.utils.backends import convert_to_dag, iter_gates

        dag = convert_to_dag(final_each)
        self.assertEqual(
            [data["gate"] for _, data in dag.nodes(data=True)],
            ["ry", "cx", "ry", "cz", "ry"],
        )
        self.assertEqual(
            [qubits for _, qubits in iter_gates(final_each)][-3:], [(), (), ()]
        )

    def test_a_native_run_reads_every_qubit(self):
        import Qdislib as qd

        # sampled as it is, only the measured qubit makes it into the keys
        self.assertEqual(qd.expectation_value([final_part, "ZI"], shots=100), 1.0)
        self.assertEqual(qd.expectation_value([final_part, "IZ"], shots=100), -1.0)
        self.assertAlmostEqual(qd.analytical_solution(final_part, "IZ"), -1.0)
        self.assertEqual(
            qd.expectation_value([final_two_registers, "IZ"], shots=100), -1.0
        )

    def test_a_measured_qubit_followed_by_gates_on_others_only(self):
        import Qdislib as qd
        from Qdislib.utils.backends import convert_to_dag

        dag = convert_to_dag(measured_before_another_qubit)
        # the numbering keeps the measurement's place, as in the other SDKs
        self.assertEqual(list(dag.nodes), ["X_1", "X_3"])
        self.assertAlmostEqual(
            qd.analytical_solution(measured_before_another_qubit, "ZZ"), 1.0
        )
        self.assertEqual(
            qd.gate_cutting(measured_before_another_qubit, [], "ZZ", shots=100), 1.0
        )

    def test_a_gate_after_a_measurement_is_refused(self):
        import Qdislib as qd

        for kernel in (mid_circuit, gate_after_whole_register):
            for run in (
                lambda: qd.analytical_solution(kernel, "ZI"),
                lambda: qd.gate_cutting(kernel, [], "ZI", shots=None),
                lambda: qd.gate_cutting(kernel, [], "ZI", shots=100),
                lambda: qd.expectation_value([kernel, "ZI"], shots=100),
            ):
                with self.subTest(kernel=kernel.name):
                    with self.assertRaisesRegex(ValueError, "after a measurement"):
                        run()

    def test_another_basis_is_refused(self):
        import Qdislib as qd

        for run in (
            lambda: qd.analytical_solution(x_basis, "Z"),
            lambda: qd.gate_cutting(x_basis, [], "Z", shots=100),
            lambda: qd.expectation_value([x_basis, "Z"], shots=100),
        ):
            with self.assertRaisesRegex(ValueError, "the Z basis"):
                run()

    def test_a_cut_kernel_keeps_its_own_measurements(self):
        import Qdislib as qd

        exact = qd.analytical_solution(plain, "ZZZ")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for observable in ("ZZZ", "ZIZ"):
                exact = qd.analytical_solution(plain, observable)
                self.assertAlmostEqual(
                    qd.gate_cutting(plain, ["CZ_4"], observable, shots=8000, seed=1),
                    exact,
                    delta=0.15,
                )
                self.assertAlmostEqual(
                    qd.gate_cutting(
                        plain,
                        ["CZ_4"],
                        observable,
                        shots=8000,
                        seed=1,
                        method="automatic",
                    ),
                    exact,
                    delta=0.15,
                )


@unittest.skipUnless(CUDAQ_AVAILABLE, "cudaq not installed")
class CudaqRegistersTest(BaseTimedTestCase):
    def test_a_bare_kernel_runs_with_shots(self):
        import cudaq

        import Qdislib as qd

        @cudaq.kernel
        def bell():
            q = cudaq.qvector(2)
            h(q[0])
            x.ctrl(q[0], q[1])

        self.assertAlmostEqual(qd.expectation_value(bell, shots=None), 1.0)
        self.assertEqual(qd.expectation_value(bell, shots=100), 1.0)

    def test_a_register_with_an_initial_state_is_refused(self):
        # qvector([0.0, 1.0]) starts in |1>: reading only its gates answers for |0>
        import cudaq

        import Qdislib as qd

        @cudaq.kernel
        def from_a_state():
            q = cudaq.qvector([0.0, 1.0])
            h(q[0])
            h(q[0])

        with self.assertRaisesRegex(ValueError, "initial state"):
            qd.gate_cutting(from_a_state, [], "Z", shots=None)

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
