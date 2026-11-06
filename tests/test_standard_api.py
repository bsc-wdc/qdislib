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

"""The calling convention every cut type shares (2.0.0): same signatures and
defaults, ``backend`` for where circuits run, ``software`` for the SDK the
three-step subcircuits are written in, one kind of subcircuit object, and the
1.x names still accepted with a DeprecationWarning."""

import importlib.util
import inspect
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None
QIBO = importlib.util.find_spec("qibo") is not None

FAMILIES = {
    "gate_cutting": {"optimized": False, "num_samples": None},
    "wire_cutting": {"optimized": False, "joint": True, "num_samples": None},
    "gate_cutting_hadamard": {"joint": True, "num_samples": None},
    "circuit_cutting": {
        "hadamard": False,
        "optimized": False,
        "joint": True,
        "num_samples": None,
    },
}


def _circuit():
    """Blocks {q0, q1} and {q2, q3} joined only by CZ_8 (gate cut), or by the wire
    of q1 between RY_7 and CZ_8 (wire cut)."""
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(4)
    for q in range(4):
        qc.h(q)  # H_1 .. H_4
    qc.cz(0, 1)  # CZ_5
    qc.cz(2, 3)  # CZ_6
    qc.ry(0.3, 1)  # RY_7
    qc.cz(1, 2)  # CZ_8
    qc.ry(0.4, 1)  # RY_9
    qc.ry(0.5, 2)  # RY_10
    qc.cz(2, 3)  # CZ_11
    qc.rx(0.2, 0)  # RX_12
    return qc


CUTS = {
    "gate_cutting": ["CZ_8"],
    "gate_cutting_hadamard": ["CZ_8"],
    "wire_cutting": [("RY_7", "CZ_8")],
    "circuit_cutting": {"gate": ["CZ_8"]},
}


def _exact(qc, observable):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(qc).expectation_value(SparsePauliOp(observable[::-1])).real
    )


def _exact_distribution(qc):
    from qiskit.quantum_info import Statevector

    return {k[::-1]: p for k, p in Statevector(qc).probabilities_dict().items()}


def _tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


class SignatureTest(BaseTimedTestCase):
    def _params(self, name):
        import Qdislib as qd

        return inspect.signature(getattr(qd, name)).parameters

    def test_same_signature_for_every_cut_type(self):
        for family, options in FAMILIES.items():
            forms = {
                family: ["circuit", "cuts", "observable"],
                f"{family}_sampling": ["circuit", "cuts"],
                f"{family}_subcircuits": ["circuit", "cuts", "observable"],
                f"{family}_sampling_subcircuits": ["circuit", "cuts"],
            }
            for name, positional in forms.items():
                with self.subTest(name=name):
                    params = self._params(name)
                    self.assertEqual(
                        [
                            p
                            for p, v in params.items()
                            if v.kind is v.POSITIONAL_OR_KEYWORD
                        ],
                        positional,
                    )
                    if "observable" in params:
                        self.assertIsNone(params["observable"].default)
                    runs = not name.endswith("subcircuits")
                    expected = (
                        {"backend": None, "shots": 1024, "seed": None, "sync": True}
                        if runs
                        else {"software": None}
                    )
                    for option, default in {**expected, **options}.items():
                        if option == "num_samples" or option in params:
                            self.assertIn(option, params)
                            self.assertEqual(params[option].default, default, option)
                    self.assertNotIn("software" if runs else "backend", params)
                    for legacy in ("qpu", "qpu_dict", "observables", "method"):
                        self.assertNotIn(legacy, params)


@unittest.skipUnless(QISKIT and QIBO, "needs qiskit and qibo")
class ExactEverywhereTest(BaseTimedTestCase):
    def test_shots_none_is_exact_for_every_family(self):
        import Qdislib as qd

        qc = _circuit()
        exact = _exact(qc, "ZZZZ")
        exact_dist = _exact_distribution(qc)
        for family in FAMILIES:
            with self.subTest(family=family):
                value = getattr(qd, family)(qc, CUTS[family], shots=None)
                self.assertAlmostEqual(value, exact, places=9)
                dist = getattr(qd, f"{family}_sampling")(qc, CUTS[family], shots=None)
                self.assertLess(_tvd(dist, exact_dist), 1e-9)

    def test_three_step_exact_for_every_family(self):
        import Qdislib as qd

        qc = _circuit()
        for family in FAMILIES:
            with self.subTest(family=family):
                subs = getattr(qd, f"{family}_subcircuits")(qc, CUTS[family], "XZZY")
                values = [qd.expectation_value(s, shots=None) for s in subs]
                self.assertAlmostEqual(
                    subs.reconstruct(values), _exact(qc, "XZZY"), places=9
                )
                subs = getattr(qd, f"{family}_sampling_subcircuits")(qc, CUTS[family])
                sampled = [qd.sample_subcircuit(s, shots=None) for s in subs]
                self.assertLess(
                    _tvd(subs.reconstruct(sampled), _exact_distribution(qc)), 1e-9
                )


@unittest.skipUnless(QISKIT and QIBO, "needs qiskit and qibo")
class SubcircuitObjectTest(BaseTimedTestCase):
    def test_every_item_is_a_subcircuit_in_the_chosen_language(self):
        import qibo
        from qiskit import QuantumCircuit

        import Qdislib as qd

        qc = _circuit()
        native = {"qiskit": QuantumCircuit, "qibo": qibo.models.Circuit}
        for family in FAMILIES:
            for form in ("_subcircuits", "_sampling_subcircuits"):
                for software, cls in native.items():
                    with self.subTest(family=family, form=form, software=software):
                        subs = getattr(qd, family + form)(
                            qc, CUTS[family], software=software
                        )
                        for s in subs:
                            self.assertIsInstance(s, qd.Subcircuit)
                            self.assertEqual(s.software, software)
                            self.assertIsInstance(s.circuit, cls)
                            self.assertGreater(s.num_qubits, 0)
                            if form == "_sampling_subcircuits":
                                self.assertIsNone(s.observable)

    def test_default_language_is_the_circuits(self):
        from qiskit import QuantumCircuit

        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(_circuit(), ["CZ_8"])
        self.assertIsInstance(subs[0].circuit, QuantumCircuit)

    def test_either_run_helper_runs_a_subcircuit(self):
        """A fragment circuit returns its correlators, which serve the expectation
        value and the distribution, so either helper runs it."""
        import Qdislib as qd

        qc = _circuit()
        subs = qd.gate_cutting_subcircuits(qc, ["CZ_8"], "ZZZZ")
        values = [qd.sample_subcircuit(s, shots=None) for s in subs]
        self.assertAlmostEqual(subs.reconstruct(values), _exact(qc, "ZZZZ"), places=9)

    def test_route_each_subcircuit_by_size(self):
        """The three-step loop picks the backend per subcircuit."""
        import Qdislib as qd

        qc = _circuit()
        subs = qd.gate_cutting_subcircuits(qc, ["CZ_8"], "ZZZZ")
        values = [
            qd.expectation_value(
                s, shots=None, backend="qiskit" if s.num_qubits <= 2 else "qibo"
            )
            for s in subs
        ]
        self.assertAlmostEqual(subs.reconstruct(values), _exact(qc, "ZZZZ"), places=9)


@unittest.skipUnless(QISKIT, "qiskit not installed")
class BackendAndSeedTest(BaseTimedTestCase):
    def test_seed_makes_runs_reproducible(self):
        import Qdislib as qd

        qc = _circuit()
        for family in FAMILIES:
            with self.subTest(family=family):
                run = getattr(qd, family)
                a = run(qc, CUTS[family], backend="qiskit", shots=500, seed=7)
                b = run(qc, CUTS[family], backend="qiskit", shots=500, seed=7)
                c = run(qc, CUTS[family], backend="qiskit", shots=500, seed=8)
                self.assertEqual(a, b)
                self.assertNotEqual(a, c)

    def test_backend_object(self):
        import Qdislib as qd

        qc = _circuit()
        aer = qd.build_backend("AER_BACKEND", qd.BackendConfig())
        value = qd.gate_cutting(qc, ["CZ_8"], backend=aer, shots=20000, seed=1)
        self.assertLess(abs(value - _exact(qc, "ZZZZ")), 0.05)
        subs = qd.wire_cutting_subcircuits(qc, [("RY_7", "CZ_8")])
        values = [qd.expectation_value(s, shots=20000, backend=aer) for s in subs]
        self.assertLess(abs(subs.reconstruct(values) - _exact(qc, "ZZZZ")), 0.08)

    def test_unknown_backend(self):
        import Qdislib as qd

        with self.assertRaises(qd.UnsupportedBackendError):
            qd.gate_cutting(_circuit(), ["CZ_8"], backend="no-such-backend")


@unittest.skipUnless(QISKIT, "qiskit not installed")
class DeprecatedNamesTest(BaseTimedTestCase):
    def test_old_keywords_still_work_and_warn(self):
        import Qdislib as qd

        qc = _circuit()
        with self.assertWarnsRegex(
            DeprecationWarning, "observables= is now observable="
        ):
            value = qd.gate_cutting(qc, ["CZ_8"], observables="ZZZZ", shots=None)
        self.assertAlmostEqual(value, _exact(qc, "ZZZZ"), places=9)
        with self.assertWarnsRegex(DeprecationWarning, "software= is now backend="):
            qd.wire_cutting(qc, [("RY_7", "CZ_8")], software="qiskit", shots=100)
        with self.assertWarnsRegex(DeprecationWarning, "Qibo platform"):
            qd.wire_cutting(qc, [("RY_7", "CZ_8")], backend="numpy", shots=None)
        with self.assertWarnsRegex(DeprecationWarning, "replaced by backend="):
            qd.gate_cutting(qc, ["CZ_8"], qpu=False, qpu_dict={}, shots=100)
        with self.assertWarns(DeprecationWarning):
            qd.gate_cutting_subcircuits(qc, gates_cut=["CZ_8"])

    def test_old_positional_order_is_refused(self):
        import Qdislib as qd

        qc = _circuit()
        with self.assertRaisesRegex(TypeError, "keyword"):
            qd.wire_cutting(qc, [("RY_7", "CZ_8")], True, "ZZZZ")
        with self.assertRaisesRegex(TypeError, "keyword-only"):
            qd.gate_cutting(qc, ["CZ_8"], "ZZZZ", 1024)

    def test_deprecated_functions(self):
        import Qdislib as qd

        qc = _circuit()
        with self.assertWarns(DeprecationWarning):
            qd.gate_wire_cutting(qc, ["CZ_8"], [], "ZZZZ", software="qiskit", shots=10)
        with self.assertWarns(DeprecationWarning):
            qd.cut_cost(["CZ_8"])
        subs = qd.gate_cutting_subcircuits(qc, ["CZ_8"])
        values = [qd.expectation_value(s, shots=None) for s in subs]
        # the deprecated positional reconstruction cannot read 2.0 lists, and says
        # what to use instead
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(TypeError, "subs.reconstruct"):
                qd.gate_cutting_subcircuit_reconstruction(values, 1, 2)
        with self.assertRaisesRegex(AttributeError, "subs.reconstruct"):
            subs.number_components
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            self.assertAlmostEqual(
                subs.reconstruct(values), _exact(qc, "Z" * qc.num_qubits), places=9
            )


@unittest.skipUnless(QISKIT, "qiskit not installed")
class InvalidInputTest(BaseTimedTestCase):
    """What cannot be right is refused, by name, before anything runs."""

    def _error(self, message, *args, **kwargs):
        import Qdislib as qd

        with self.assertRaisesRegex(qd.InvalidCutError, message):
            qd.circuit_cutting(_circuit(), *args, shots=None, **kwargs)

    def test_a_cut_of_no_gate_of_the_circuit(self):
        self._error("No gate .*CZ_99", ["CZ_99"])
        self._error("No gate .*CZ_99", {"hadamard": ["CZ_99"]})
        self._error("No gate .*CZ_99", [("RY_7", "CZ_99")])

    def test_a_cut_given_twice(self):
        self._error("CZ_8.* cut twice", ["CZ_8", "CZ_8"])
        self._error("CZ_8.* cut twice", {"gate": ["CZ_8"], "hadamard": ["CZ_8"]})
        self._error("cut twice", [("RY_7", "CZ_8"), ("RY_7", "CZ_8")])

    def test_a_wire_between_gates_that_share_no_qubit_cuts_nothing(self):
        import Qdislib as qd

        qc = _circuit()
        names = {f"{n}_{i}": q for i, (n, q) in enumerate(qd.iter_gates(qc), 1)}
        source, target = next(
            (a, b)
            for a, qa in names.items()
            for b, qb in names.items()
            if qa and qb and not set(qa) & set(qb)
        )
        self._error("do not join", [(source, target)])

    def test_an_observable_with_a_letter_per_qubit(self):
        import Qdislib as qd

        qc = _circuit()
        width = qc.num_qubits
        cut = ["CZ_8"]
        # a shorter one read as if the rest were I, a longer one as idle qubits
        for wrong in ("Z" * (width - 1), "Z" * (width + 1), "X" * (width + 1)):
            for call in (qd.circuit_cutting, qd.circuit_cutting_subcircuits):
                with self.assertRaisesRegex(ValueError, "one per qubit"):
                    call(qc, cut, wrong)
        with self.assertRaisesRegex(ValueError, "one per qubit"):
            qd.gate_cutting(qc, cut, [(0.5, "Z" * (width - 1))])
        self.assertAlmostEqual(
            qd.circuit_cutting(qc, cut, "Z" * width, shots=None),
            _exact(qc, "Z" * width),
            places=9,
        )

    def test_shots(self):
        import Qdislib as qd

        for shots in (0, -5):
            with self.assertRaisesRegex(ValueError, "shots must be"):
                qd.circuit_cutting(_circuit(), ["CZ_8"], shots=shots)
            with self.assertRaisesRegex(ValueError, "shots must be"):
                qd.circuit_cutting_sampling(_circuit(), ["CZ_8"], shots=shots)

    def test_no_circuit(self):
        import Qdislib as qd

        with self.assertRaisesRegex(TypeError, "circuit"):
            qd.circuit_cutting("not a circuit", [], "ZZ")
        with self.assertRaisesRegex(TypeError, "ircuit"):
            qd.circuit_cutting(object(), [])


@unittest.skipUnless(QISKIT, "qiskit not installed")
class DegenerateCircuitTest(BaseTimedTestCase):
    def _empty(self):
        from qiskit import QuantumCircuit

        return QuantumCircuit(2)

    def test_no_gates(self):
        import Qdislib as qd

        for shots in (None, 100):
            self.assertAlmostEqual(
                qd.circuit_cutting(self._empty(), [], "ZZ", shots=shots), 1.0
            )
            self.assertEqual(
                dict(qd.circuit_cutting_sampling(self._empty(), [], shots=shots)),
                {"00": 1.0, "01": 0.0, "10": 0.0, "11": 0.0},
            )
        subs = qd.circuit_cutting_subcircuits(self._empty(), [], "ZI")
        self.assertAlmostEqual(
            subs.reconstruct([qd.expectation_value(s, shots=None) for s in subs]), 1.0
        )

    def test_an_observable_of_identities_with_free_parameters(self):
        from qiskit import QuantumCircuit
        from qiskit.circuit import Parameter

        import Qdislib as qd

        theta = Parameter("t")
        qc = QuantumCircuit(2)
        qc.ry(theta, 0)
        qc.cz(0, 1)
        values = qd.circuit_cutting(
            qc, ["CZ_2"], "II", shots=None, parameter_values=[[0.1], [0.2]]
        )
        self.assertEqual(list(values), [1.0, 1.0])


@unittest.skipUnless(QISKIT, "qiskit not installed")
class ReconstructAnotherObservableTest(BaseTimedTestCase):
    def test_lower_case_reads_as_upper_case(self):
        # the basis check took "iz..." and the mask read each i as a Z
        import Qdislib as qd

        qc = _circuit()
        width = qc.num_qubits
        subs = qd.circuit_cutting_subcircuits(qc, ["CZ_8"], "Z" * width)
        results = [
            qd.sample_subcircuit(s, shots=500, seed=n) for n, s in enumerate(subs)
        ]
        target = "Z" + "I" * (width - 1)
        self.assertEqual(
            subs.reconstruct(results, target), subs.reconstruct(results, target.lower())
        )


class OverheadTest(BaseTimedTestCase):
    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_one_cost_function_for_every_kind(self):
        import Qdislib as qd

        qc = _circuit()
        gate = qd.cutting_overhead(qc, {"gate": ["CZ_8"]})
        self.assertEqual((gate["subcircuits"], gate["gamma"]), (12, 3))
        had = qd.cutting_overhead(qc, {"hadamard": ["CZ_8"]})
        self.assertEqual(had["gamma"], 3)
        wire = qd.cutting_overhead(qc, {"wire": [("RY_7", "CZ_8")]}, optimized=True)
        self.assertEqual((wire["subcircuits"], wire["gamma"]), (7, 4))


if __name__ == "__main__":
    unittest.main()
