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

"""Free parameters: a circuit cut once with symbolic angles, then run for many
sets of values (``subs.bind(values)``, ``parameter_values=``)."""

import importlib.util
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

VALUES = [[0.4, 1.1, -0.7], [0.9, -0.2, 1.5], {"a": 2.0, "b": 0.1, "c": 0.3}]
CUTS = {
    "gate": {"gate": ["CZ_4"]},
    "hadamard": {"hadamard": ["CZ_4"]},
    "wire": {"wire": [("CZ_4", "RY_5")]},
    "mixed": {"gate": ["CZ_4"], "wire": [("CX_3", "CZ_4")]},
}
HAMILTONIAN = [(0.5, "ZZZZ"), (0.7, "XIZI"), (-0.3, "IIZZ"), (0.2, "IIII")]


def _circuit():
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter

    a, b, c = Parameter("a"), Parameter("b"), Parameter("c")
    qc = QuantumCircuit(4)
    qc.ry(a, 0)
    qc.rx(b, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(c, 2)
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    return qc


def _bound(qc, values):
    """``qc`` with ``values`` (a list in parameter order, or a name mapping)."""
    params = {p.name: p for p in qc.parameters}
    if isinstance(values, dict):
        return qc.assign_parameters({params[k]: v for k, v in values.items()})
    return qc.assign_parameters(values)


@unittest.skipUnless(QISKIT, "qiskit required")
class ParametersTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        warnings.simplefilter("ignore")
        import Qdislib as qd

        self.qd, self.qc = qd, _circuit()

    def test_bind_is_exact_for_every_cut_type(self):
        qd = self.qd
        for name, cuts in CUTS.items():
            for observable in ("ZIZZ", HAMILTONIAN):
                subs = qd.circuit_cutting_subcircuits(self.qc, cuts, observable)
                self.assertEqual(subs.parameters, ("a", "b", "c"))
                for values in VALUES:
                    bound = subs.bind(values)
                    value = bound.reconstruct(
                        [qd.expectation_value(s, shots=None) for s in bound]
                    )
                    exact = qd.analytical_solution(_bound(self.qc, values), observable)
                    self.assertAlmostEqual(value, exact, places=9, msg=name)

    def test_parameter_values_match_cutting_each_bound_circuit(self):
        from Qdislib.core.cutting_algorithms._reconstruction import next_seed

        qd = self.qd
        for name, cuts in CUTS.items():
            for options in ({}, {"circuits_per_task": 5}, {"error": 30}):
                values = qd.circuit_cutting(
                    self.qc,
                    cuts,
                    HAMILTONIAN,
                    shots=3000,
                    seed=11,
                    parameter_values=VALUES,
                    **options,
                )
                self.assertEqual(len(values), len(VALUES))
                for i, v in enumerate(VALUES):
                    reference = qd.circuit_cutting(
                        _bound(self.qc, v),
                        cuts,
                        HAMILTONIAN,
                        shots=3000,
                        seed=next_seed(11, "parameter_values", i),
                        **options,
                    )
                    self.assertAlmostEqual(values[i], reference, places=12, msg=name)
                    if "error" in options:
                        self.assertAlmostEqual(
                            values[i].error, reference.error, places=12
                        )
            distributions = qd.circuit_cutting_sampling(
                self.qc,
                cuts,
                shots=3000,
                seed=5,
                parameter_values=VALUES,
                circuits_per_task=4,
            )
            for i, v in enumerate(VALUES):
                reference = qd.circuit_cutting_sampling(
                    _bound(self.qc, v),
                    cuts,
                    shots=3000,
                    seed=next_seed(5, "parameter_values", i),
                    circuits_per_task=4,
                )
                for bits in set(distributions[i]) | set(reference):
                    self.assertAlmostEqual(
                        distributions[i].get(bits, 0), reference.get(bits, 0), places=12
                    )

    def test_every_cut_type_takes_parameter_values(self):
        qd = self.qd
        for fn, cuts in (
            (qd.gate_cutting, ["CZ_4"]),
            (qd.gate_cutting_hadamard, ["CZ_4"]),
            (qd.wire_cutting, [("CZ_4", "RY_5")]),
        ):
            values = fn(self.qc, cuts, "ZIZZ", shots=None, parameter_values=VALUES)
            for value, v in zip(values, VALUES):
                exact = qd.analytical_solution(_bound(self.qc, v), "ZIZZ")
                self.assertAlmostEqual(value, exact, places=9, msg=fn.__name__)
        for fn, cuts in (
            (qd.gate_cutting_sampling, ["CZ_4"]),
            (qd.gate_cutting_hadamard_sampling, ["CZ_4"]),
            (qd.wire_cutting_sampling, [("CZ_4", "RY_5")]),
        ):
            distributions = fn(self.qc, cuts, shots=None, parameter_values=VALUES)
            self.assertEqual(len(distributions), len(VALUES), fn.__name__)

    def test_the_cache_keeps_bindings_apart(self):
        qd = self.qd
        cuts = CUTS["gate"]
        exact = [qd.analytical_solution(_bound(self.qc, v), "ZIII") for v in VALUES[:2]]
        with qd.use_cache("memory"):
            first = qd.circuit_cutting(
                self.qc, cuts, "ZIII", shots=None, parameter_values=VALUES[:2]
            )
            again = qd.circuit_cutting(
                self.qc, cuts, "ZIII", shots=None, parameter_values=VALUES[:2]
            )
        for values in (first, again):
            for value, reference in zip(values, exact):
                self.assertAlmostEqual(value, reference, places=9)

    def test_values_are_bound_once_per_fragment_and_set(self):
        import Qdislib.core.cutting_algorithms.parameters as parameters

        calls = []
        original = parameters.bind_dag

        def counting(dag, values):
            calls.append(1)
            return original(dag, values)

        parameters.bind_dag = counting
        try:
            for per_task in (None, 3):
                calls.clear()
                self.qd.circuit_cutting(
                    self.qc,
                    CUTS["gate"],
                    "ZIZZ",
                    shots=100,
                    seed=1,
                    parameter_values=VALUES,
                    circuits_per_task=per_task,
                )
                # 2 fragments x 3 sets, not once per circuit (12 of them per set)
                self.assertEqual(len(calls), 2 * len(VALUES), per_task)
        finally:
            parameters.bind_dag = original

    def test_errors(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit import Parameter

        qd = self.qd
        subs = qd.circuit_cutting_subcircuits(self.qc, CUTS["gate"], "ZZZZ")
        with self.assertRaisesRegex(ValueError, "free parameters"):
            qd.expectation_value(subs[0], shots=100)
        with self.assertRaisesRegex(ValueError, "expected 3 values"):
            subs.bind([1.0])
        with self.assertRaisesRegex(ValueError, "missing"):
            subs.bind({"a": 1.0})
        with self.assertRaisesRegex(ValueError, "list of value sets"):
            qd.circuit_cutting(
                self.qc, CUTS["gate"], "ZZZZ", parameter_values={"a": 1, "b": 2, "c": 3}
            )
        q = QuantumCircuit(3)
        q.ry(0.2, 0)
        q.cx(0, 1)
        q.rzz(Parameter("t"), 1, 2)
        for kind in ("gate", "hadamard"):
            with self.assertRaisesRegex(ValueError, "Cut gate"):
                qd.circuit_cutting_subcircuits(q, {kind: ["RZZ_3"]}, "ZZZ")
        # an unbound list still gives its circuits, symbolic, to bind yourself
        circuit = subs.circuits("qiskit")[0]
        self.assertTrue(circuit.parameters)

    def test_symbolic_ry_keeps_golden_wire_cuts(self):
        from qiskit.circuit import Parameter

        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan
        from tests.test_cut_costs import _chain, _chain_cuts, _exact

        numeric = _chain()
        symbolic = numeric.copy_empty_like()
        angles = []
        for ins in numeric.data:
            op = ins.operation
            if op.name == "ry":
                p = Parameter(f"t{len(angles):02d}")
                angles.append(float(op.params[0]))
                symbolic.ry(p, ins.qubits[0])
            else:
                symbolic.append(ins)
        cuts = _chain_cuts(numeric)
        plan = _CutPlan(symbolic, cuts, "ZZZZZZ", False, False, True)
        self.assertEqual([len(d) for d in plan.wire_digits], [6, 6])  # still golden
        value = self.qd.circuit_cutting(
            symbolic, cuts, "ZZZZZZ", shots=None, parameter_values=[angles]
        )[0]
        self.assertAlmostEqual(value, _exact(numeric, "ZZZZZZ"), places=9)


if __name__ == "__main__":
    unittest.main()
