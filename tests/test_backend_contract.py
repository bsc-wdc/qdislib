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

"""A backend needs one method: ``evaluate_samples``. Cutting derives the rest."""

import importlib.util
import pickle
import types
import unittest
from unittest import mock

from tests import BaseTimedTestCase

AER = importlib.util.find_spec("qiskit_aer") is not None

from Qdislib.core.backend.base_backend import (  # noqa: E402
    BackendConfig,
    BaseBackend,
    sample_with_qiskit,
)


class Minimal(BaseBackend):
    """The whole of a new backend, as the user guide shows it (module level: it has
    to pickle)."""

    def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
        from qiskit_aer import AerSimulator

        return sample_with_qiskit(
            dag,
            num_qubits,
            lambda circuit: AerSimulator()
            .run(circuit, shots=shots, seed_simulator=seed)
            .result()
            .get_counts(),
        )


def _circuit():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(-0.7, 2)  # qubits 1 and 2 keep evolving: their cut points are deferred
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    return qc


def _terminal():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(2)
    qc.ry(0.4, 0)
    qc.rx(1.1, 1)
    qc.cz(0, 1)  # nothing after the cut: read from the final bits
    return qc


def _exact(qc):
    from qiskit.quantum_info import Statevector

    return {k[::-1]: v for k, v in Statevector(qc).probabilities_dict().items()}


@unittest.skipUnless(AER, "qiskit-aer required")
class MinimalBackendTest(BaseTimedTestCase):
    def setUp(self):
        super().setUp()
        import Qdislib as qd

        self.qd = qd
        self.device = Minimal("MINIMAL", BackendConfig(num_qubits=8))

    def assertClose(self, got, want, delta=0.05):
        for key in set(got) | set(want):
            self.assertAlmostEqual(
                got.get(key, 0.0), want.get(key, 0.0), delta=delta, msg=key
            )

    def test_distribution_for_every_cut_type(self):
        qc, want = _circuit(), _exact(_circuit())
        cuts = {
            "gate": {"gate": ["CZ_4"]},
            "hadamard": {"hadamard": ["CZ_4"]},
            "wire": {"wire": [("CZ_4", "RY_5")]},
        }
        for name, cut in cuts.items():
            with self.subTest(cut=name):
                got = self.qd.circuit_cutting_sampling(
                    qc, cut, backend=self.device, shots=40000, seed=3
                )
                self.assertClose(got, want)

    def test_a_cut_point_with_gates_after_it_is_deferred(self):
        # the gate cut's own rz follows its measurement: one more qubit, read at the end
        qc = _terminal()
        widths = []
        evaluate = self.device.evaluate_samples

        def spy(dag, num_qubits, *a, **k):
            widths.append(num_qubits)
            return evaluate(dag, num_qubits, *a, **k)

        self.device.evaluate_samples = spy
        got = self.qd.circuit_cutting_sampling(
            qc, {"gate": ["CZ_3"]}, backend=self.device, shots=40000, seed=1
        )
        self.assertClose(got, _exact(qc))
        self.assertEqual(set(widths), {1, 2})

    def test_a_terminal_cut_point_is_read_from_the_final_bits(self):
        import networkx as nx

        dag = nx.DiGraph()
        dag.add_node("RY_1", gate="ry", qubits=[0], parameters=[0.3])
        dag.add_node("CX_2", gate="cx", qubits=[0, 1], parameters=[])
        dag.add_node("MEASURE_3", gate="measure", qubits=[1], parameters=[])
        dag.add_edge("RY_1", "CX_2")
        dag.add_edge("CX_2", "MEASURE_3")
        widths = []

        class Fixed(BaseBackend):
            def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
                widths.append(num_qubits)
                return {"00": 3, "11": 5, "10": 2}

        joint = Fixed("FIXED", BackendConfig(num_qubits=2)).evaluate_joint_counts(
            dag, 2, [1], 10
        )
        self.assertEqual(widths, [2])  # no extra qubit
        self.assertEqual(
            joint, {((0, 0), (0,)): 3, ((1, 1), (1,)): 5, ((1, 0), (0,)): 2}
        )

    def test_expectation_value(self):
        from qiskit.quantum_info import SparsePauliOp, Statevector

        qc = _circuit()
        want = Statevector(qc).expectation_value(SparsePauliOp("ZZZZ")).real
        got = self.qd.gate_cutting(
            qc, ["CZ_4"], "ZZZZ", backend=self.device, shots=40000
        )
        self.assertAlmostEqual(got, want, delta=0.06)

    def test_it_can_go_to_a_task(self):
        copy = pickle.loads(pickle.dumps(self.device))
        got = self.qd.circuit_cutting_sampling(
            _terminal(), {"hadamard": ["CZ_3"]}, backend=copy, shots=20000, seed=2
        )
        self.assertClose(got, _exact(_terminal()))


class ContractTest(BaseTimedTestCase):
    def test_nothing_is_abstract(self):
        BaseBackend("NOTHING", BackendConfig(num_qubits=1))

    def test_what_is_missing_says_so(self):
        device = BaseBackend("NOTHING", BackendConfig(num_qubits=1))
        with self.assertRaisesRegex(NotImplementedError, "NOTHING.*evaluate_samples"):
            device.evaluate_samples(None, 1)
        with self.assertRaisesRegex(NotImplementedError, "NOTHING.*run_estimator"):
            device.run_estimator("OPENQASM", False, "Z", 1)

    def test_transpile_defaults_to_the_circuit(self):
        device = BaseBackend("NOTHING", BackendConfig(num_qubits=1))
        circuit = object()
        self.assertIs(device.transpile(circuit), circuit)


@unittest.skipUnless(AER, "qiskit-aer required")
class IBMOptionsTest(BaseTimedTestCase):
    def test_options_reach_the_sampler(self):
        from qiskit import QuantumCircuit
        from qiskit.providers.fake_provider import GenericBackendV2
        from Qdislib.core.backend.qpu_backends import IBMBackend
        from Qdislib.core.cutting_algorithms import gate_cutting as gc
        from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

        device = GenericBackendV2(2, basis_gates=["cz", "rz", "sx", "x"], seed=1)
        options = {"twirling": {"enable_gates": True}}
        ibm = IBMBackend(
            BackendConfig(
                num_qubits=2,
                extra={"batch": device, "backend": device, "options": options},
            )
        )
        qc = QuantumCircuit(2)
        qc.x(0)
        dag = circuit_qiskit_to_dag(qc)
        stop = RuntimeError("stop before running")
        for call in (
            lambda: ibm.evaluate_joint_counts(dag, 2, (), 10),
            lambda: ibm.evaluate_samples(dag, 2, 10),
            lambda: ibm.evaluate_expectation(dag, 2, shots=10),
        ):
            with mock.patch.object(gc, "Sampler") as sampler:
                sampler.return_value.run.side_effect = stop
                with self.assertRaises(RuntimeError):
                    call()
                self.assertEqual(sampler.call_args.kwargs["options"], options)


@unittest.skipUnless(AER, "qiskit-aer required")
class AerOptionsTest(BaseTimedTestCase):
    """What a backend object is configured with applies to the circuits a cut runs."""

    def _noisy(self):
        from qiskit_aer.noise import NoiseModel, depolarizing_error

        import Qdislib as qd

        noise = NoiseModel()
        noise.add_all_qubit_quantum_error(
            depolarizing_error(0.3, 1), ["x", "h", "rz", "z", "s", "sdg", "u"]
        )
        return qd.build_backend("AER_BACKEND", qd.BackendConfig(noise_model=noise))

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.x(0)
        qc.cz(0, 1)
        qc.x(1)
        return qc

    def test_a_cut_runs_with_the_noise_model(self):
        import Qdislib as qd

        qc = self._circuit()
        kept = {"gate": ["CZ_2"]}
        clean = qd.circuit_cutting(qc, kept, "ZZ", shots=4000, seed=1)
        self.assertAlmostEqual(clean, 1.0, places=6)
        noisy = qd.circuit_cutting(
            qc, kept, "ZZ", backend=self._noisy(), shots=4000, seed=1
        )
        self.assertLess(noisy, 0.8)

    def test_a_hadamard_cut_runs_with_the_noise_model(self):
        import Qdislib as qd

        noisy = qd.circuit_cutting(
            self._circuit(),
            {"hadamard": ["CZ_2"]},
            "ZZ",
            backend=self._noisy(),
            shots=4000,
            seed=1,
        )
        self.assertLess(noisy, 0.8)

    def test_evaluate_expectation_runs_with_the_noise_model(self):
        from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

        dag = circuit_qiskit_to_dag(self._circuit())
        self.assertLess(
            self._noisy().evaluate_expectation(dag, 2, (), 4000, seed=1), 0.8
        )

    def test_gpu_runs_on_the_backends_simulator(self):
        # the cutting code's own simulator is a CPU one: a run asked for on a GPU
        # went there, silently
        from qiskit_aer import AerSimulator

        import Qdislib as qd

        qc = _circuit()
        exact = qd.analytical_solution(qc, "ZZZZ")
        backend = qd.build_backend("AER_BACKEND", qd.BackendConfig(use_gpu=True))
        backend._simulator_cache = mock.Mock(wraps=AerSimulator())
        value = qd.circuit_cutting(
            qc, ["CZ_4"], "ZZZZ", backend=backend, shots=40000, seed=3
        )
        self.assertTrue(backend._simulator_cache.run.called)
        self.assertLess(abs(value - exact), 0.1)

    def test_precision_keeps_the_value_and_a_cut_point_read_at_the_end(self):
        import Qdislib as qd

        qc = _circuit()
        exact = qd.analytical_solution(qc, "ZIZI")
        single = qd.build_backend("AER_BACKEND", qd.BackendConfig(precision="single"))
        value = qd.circuit_cutting(
            qc, ["CZ_4"], "ZIZI", backend=single, shots=40000, seed=3
        )
        self.assertLess(abs(value - exact), 0.1)

    def test_goes_to_a_task_after_running(self):
        import Qdislib as qd
        from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

        backend = self._noisy()
        dag = circuit_qiskit_to_dag(self._circuit())
        before = backend.evaluate_samples(dag.copy(), 2, 100, seed=1)
        again = pickle.loads(pickle.dumps(backend))
        self.assertEqual(again.evaluate_samples(dag.copy(), 2, 100, seed=1), before)
        self.assertIsInstance(
            again, type(qd.build_backend("AER_BACKEND", qd.BackendConfig()))
        )


@unittest.skipUnless(importlib.util.find_spec("qibo") is not None, "qibo required")
class QiboOptionsTest(BaseTimedTestCase):
    def test_a_cut_runs_on_the_platform(self):
        from qibo import gates, models

        import Qdislib as qd
        from Qdislib.utils.graph_qibo import circuit_qibo_to_dag

        qc = models.Circuit(1)
        qc.add(gates.X(0))
        backend = qd.build_backend("QIBO_BACKEND", qd.BackendConfig(platform="mine"))
        with mock.patch("qibo.set_backend") as set_backend:
            joint = backend.evaluate_joint_counts(circuit_qibo_to_dag(qc), 1, (), 10)
        set_backend.assert_called_with("mine")
        self.assertEqual(joint, {((1,), ()): 10})


@unittest.skipUnless(importlib.util.find_spec("cudaq") is not None, "cudaq required")
class CudaqOptionsTest(BaseTimedTestCase):
    def test_a_cut_runs_on_the_target(self):
        # it went to the "nvidia" target, whatever the backend was given
        import cudaq
        import networkx as nx

        import Qdislib as qd

        dag = nx.DiGraph()
        dag.add_node("X_1", gate="x", qubits=(0,), parameters=())
        dag.add_node("M_2", gate="measure", qubits=(0,), parameters=())
        dag.add_edge("X_1", "M_2", color="blue")
        backend = qd.build_backend("CUDAQ_BACKEND", qd.BackendConfig(target="qpp-cpu"))
        for prep in ((), (0,)):
            with mock.patch.object(
                cudaq, "set_target", wraps=cudaq.set_target
            ) as chosen:
                joint = backend.evaluate_joint_counts(dag.copy(), 1, prep, 10, seed=1)
            self.assertEqual(
                {call.args[0] for call in chosen.call_args_list}, {"qpp-cpu"}
            )
            self.assertEqual(sum(joint.values()), 10)


if __name__ == "__main__":
    unittest.main()
