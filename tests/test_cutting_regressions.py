#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Regression tests for two cutting-correctness fixes.

1. ``wire_cutting`` forwards the user's ``shots``/``backend`` into subcircuit
   evaluation on the standard (non-optimized) path (they were dropped, so every
   subcircuit ran at the hardcoded default).
2. A gate/wire cut that isolates a qubit into a component whose only marker is an
   Observable-I node (stripped before conversion) left an EMPTY component, which
   crashed ``_update_qubits_serie`` with ``ValueError: max() arg is empty``. An
   empty component contributes <I> = 1 and must reconstruct without crashing.
"""

import importlib.util
import unittest
import warnings

import numpy as np

from tests import BaseTimedTestCase

QIBO = importlib.util.find_spec("qibo") is not None
QISKIT = importlib.util.find_spec("qiskit") is not None


def _z_oracle(qibo_circuit, obs):
    """Independent statevector oracle for a Z/I Pauli string (qubit-0-first)."""
    psi = qibo_circuit().state(numpy=True)
    z = np.diag([1.0, -1.0]).astype(complex)
    eye = np.eye(2)
    op = np.array([[1]], dtype=complex)
    for ch in obs:
        op = np.kron(op, z if ch == "Z" else eye)
    return float(np.real(psi.conj() @ op @ psi))


@unittest.skipUnless(QIBO, "qibo not installed")
class WireCuttingForwardsShotsTest(BaseTimedTestCase):
    """The standard wire-cutting path must honour the user's ``shots`` argument."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _spy_shots(self, circuit, cut, user_shots):
        import Qdislib.core.cutting_algorithms.mixed_cutting as mc
        import Qdislib.core.cutting_algorithms.wire_cutting as wc

        seen = []
        original = mc._run_fragment

        def spy(template, spec, runner, shots, *args, **kwargs):
            seen.append(shots)
            return original(template, spec, runner, shots, *args, **kwargs)

        mc._run_fragment = spy
        try:
            wc.wire_cutting(circuit, cut, shots=user_shots, backend="qibo")
        finally:
            mc._run_fragment = original
        return seen

    def test_single_component_forwards_shots(self):
        from qibo import models, gates

        c = models.Circuit(3)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))  # CZ_2
        c.add(gates.RY(1, 0.5))  # RY_3
        c.add(gates.CZ(1, 2))
        seen = self._spy_shots(c, [("CZ_2", "RY_3")], user_shots=7777)
        self.assertTrue(seen, "no subcircuits were evaluated")
        # Every subcircuit must run at the user's shot count, not the old hardcoded
        # _generate_wire_cutting default (10000).
        self.assertEqual(set(seen), {7777})

    def test_multi_component_forwards_shots(self):
        from qibo import models, gates

        c = models.Circuit(4)
        c.add(gates.H(0))
        c.add(gates.CZ(0, 1))  # CZ_2  (component A: q0, q1)
        c.add(gates.RY(1, 0.4))  # RY_3
        c.add(gates.H(2))
        c.add(gates.CZ(2, 3))  # component B: q2, q3
        seen = self._spy_shots(c, [("CZ_2", "RY_3")], user_shots=7777)
        self.assertTrue(seen, "no subcircuits were evaluated")
        self.assertEqual(set(seen), {7777})


@unittest.skipUnless(QIBO, "qibo not installed")
class IsolatedQubitAfterCutTest(BaseTimedTestCase):
    """A cut that isolates a qubit must not crash when that side is observed with a
    single-qubit Z (previously an empty component hit ``max() arg is empty``)."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _circuit(self):
        from qibo import models, gates

        # Cutting CZ(0, 2) isolates qubit 1 (carries RY_2, no other edge) as its own
        # component; observing "ZII"/"IZI" leaves the other component's observable
        # all-identity -> an empty subcircuit after Observable-I stripping.
        c = models.Circuit(3)
        c.add(gates.RY(0, 0.9))  # RY_1
        c.add(gates.RY(1, 1.1))  # RY_2 -> isolated component after the cut
        c.add(gates.CZ(0, 2))  # CZ_3 -> gate cut
        return c

    def test_single_z_on_isolated_side_reconstructs(self):
        import Qdislib as qd

        for obs in ("ZII", "IZI", "IIZ", "ZIZ", "ZZZ"):
            oracle = _z_oracle(self._circuit(), obs)
            got = qd.gate_cutting(
                self._circuit(), ["CZ_3"], observable=obs, backend="qibo", shots=8192
            )
            self.assertLess(abs(got - oracle), 0.15, f"qibo {obs}")

    @unittest.skipUnless(QISKIT, "qiskit not installed")
    def test_single_z_on_isolated_side_reconstructs_qiskit(self):
        import Qdislib as qd

        for obs in ("ZII", "IZI"):  # the two that used to crash
            oracle = _z_oracle(self._circuit(), obs)
            got = qd.gate_cutting(
                self._circuit(),
                ["CZ_3"],
                observable=obs,
                backend="qiskit",
                shots=8192,
            )
            self.assertLess(abs(got - oracle), 0.15, f"qiskit {obs}")

    def test_update_qubits_serie_empty_component(self):
        # Unit check: an empty component reports one idle |0> qubit (contributes 1.0)
        # instead of raising on max()/min() of an empty set.
        import networkx
        from Qdislib.utils.graph_qibo import _update_qubits_serie

        s, highest, smallest = _update_qubits_serie(networkx.DiGraph())
        self.assertEqual((highest, smallest), (1, 0))


def _exact(circuit, pauli):
    """<pauli> of a Qiskit circuit, pauli qubit-0-first."""
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(pauli[::-1])).real
    )


def _two_blocks():
    """Qubits 0-1 and qubits 2-3 never interact: two disconnected pieces."""
    from qiskit import QuantumCircuit

    circuit = QuantumCircuit(4)
    circuit.ry(0.5, 0)
    circuit.ry(0.7, 1)
    circuit.cz(0, 1)
    circuit.rx(0.4, 1)
    circuit.ry(1.1, 2)
    circuit.ry(0.3, 3)
    circuit.cz(2, 3)
    circuit.rx(0.9, 2)
    return circuit


def _three_step(subcircuits, shots=60000):
    import Qdislib as qd

    return subcircuits.reconstruct(
        [qd.expectation_value(s, shots=shots) for s in subcircuits]
    )


@unittest.skipUnless(QIBO and QISKIT, "qibo and qiskit required")
class SubcircuitsGroupingTest(BaseTimedTestCase):
    """bsc-wdc/qdislib#9: a piece that is not cut must not shift the grouping."""

    def test_gate_cut_with_an_uncut_piece(self):
        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(_two_blocks(), ["CZ_3"], observable="ZZZZ")
        self.assertLess(abs(_three_step(subs) - _exact(_two_blocks(), "ZZZZ")), 0.05)

    def test_gate_cut_on_both_pieces(self):
        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(
            _two_blocks(), ["CZ_3", "CZ_7"], observable="ZZZZ"
        )
        self.assertLess(abs(_three_step(subs) - _exact(_two_blocks(), "ZZZZ")), 0.05)

    def test_wire_cut_with_an_uncut_piece(self):
        import Qdislib as qd

        subs = qd.wire_cutting_subcircuits(
            _two_blocks(), [("RY_1", "CZ_3")], observable="ZZZZ"
        )
        self.assertLess(abs(_three_step(subs) - _exact(_two_blocks(), "ZZZZ")), 0.05)

    def test_cache_with_an_uncut_piece(self):
        import Qdislib as qd

        value = qd.gate_cutting(
            _two_blocks(),
            ["CZ_3"],
            observable="ZZZZ",
            shots=60000,
            cache=qd.InMemoryCircuitCache(),
        )
        self.assertLess(abs(value - _exact(_two_blocks(), "ZZZZ")), 0.05)


@unittest.skipUnless(QIBO and QISKIT, "qibo and qiskit required")
class SubcircuitsPauliSumTest(BaseTimedTestCase):
    """bsc-wdc/qdislib#6: the three-step workflow honours a weighted Pauli sum."""

    H = [(0.5, "ZZZZ"), (-1.0, "IZIZ"), (0.3, "XIZI")]

    def _exact_h(self):
        return sum(c * _exact(_two_blocks(), p) for c, p in self.H)

    def test_gate_cutting_subcircuits(self):
        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(_two_blocks(), ["CZ_3"], observable=self.H)
        self.assertLess(abs(_three_step(subs) - self._exact_h()), 0.05)

    def test_wire_cutting_subcircuits(self):
        import Qdislib as qd

        subs = qd.wire_cutting_subcircuits(
            _two_blocks(), [("RY_1", "CZ_3")], observable=self.H
        )
        self.assertLess(abs(_three_step(subs) - self._exact_h()), 0.05)

    def test_reconstruct_checks_the_value_count(self):
        import Qdislib as qd

        subs = qd.gate_cutting_subcircuits(_two_blocks(), ["CZ_3"], observable=self.H)
        with self.assertRaises(ValueError):
            subs.reconstruct([0.0] * (len(subs) - 1))

    def test_hadamard_takes_a_pauli_sum(self):
        import Qdislib as qd

        value = qd.gate_cutting_hadamard(
            _two_blocks(), ["CZ_3"], observable=self.H, shots=None
        )
        self.assertAlmostEqual(value, self._exact_h(), places=8)
        with self.assertRaises(ValueError):  # not a Pauli letter
            qd.gate_cutting_hadamard(
                _two_blocks(), ["CZ_3"], observable="AZZZ", shots=None
            )


def _serial_wire():
    """Qubit 0 runs H, RY, RX, CZ, RY: two wire cuts in series fit on it."""
    from qiskit import QuantumCircuit

    circuit = QuantumCircuit(2)
    circuit.h(0)
    circuit.ry(0.7, 0)
    circuit.rx(0.4, 0)
    circuit.h(1)
    circuit.cz(0, 1)
    circuit.ry(0.3, 0)
    return circuit


@unittest.skipUnless(QIBO and QISKIT, "qibo and qiskit required")
class SerialWireCutsTest(BaseTimedTestCase):
    """bsc-wdc/qdislib#7: two cuts one after the other on the same wire."""

    CUTS = [("H_1", "RY_2"), ("RY_2", "RX_3")]

    def test_expectation_paths(self):
        import Qdislib as qd

        exact = _exact(_serial_wire(), "ZZ")
        for kwargs in ({}, {"optimized": True}):
            value = qd.wire_cutting(
                _serial_wire(),
                self.CUTS,
                observable="ZZ",
                backend="qiskit",
                shots=40000,
                **kwargs,
            )
            self.assertLess(abs(value - exact), 0.05, kwargs)
        for cuts in (self.CUTS[::-1],):
            value = qd.wire_cutting(
                _serial_wire(), cuts, observable="ZZ", backend="qiskit", shots=40000
            )
            self.assertLess(abs(value - exact), 0.05, cuts)
        subs = qd.wire_cutting_subcircuits(_serial_wire(), self.CUTS, observable="ZZ")
        self.assertLess(abs(_three_step(subs, shots=40000) - exact), 0.05)

    def test_sampling(self):
        import Qdislib as qd
        from qiskit.quantum_info import Statevector

        # Qiskit keys are little-endian; Qdislib's are qubit-0-first.
        exact = {
            key[::-1]: p
            for key, p in Statevector(_serial_wire()).probabilities_dict().items()
        }
        counts = qd.wire_cutting_sampling(
            _serial_wire(), self.CUTS, shots=40000, backend="qiskit"
        )
        total = sum(counts.values())
        for key, p in exact.items():
            self.assertLess(abs(counts.get(key, 0) / total - p), 0.03, key)


@unittest.skipUnless(QIBO and QISKIT, "qibo and qiskit required")
class FinalMeasurementTest(BaseTimedTestCase):
    """bsc-wdc/qdislib#4: a circuit's final measurements change nothing."""

    @staticmethod
    def _qiskit(measure):
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.ry(0.7, 0)
        circuit.ry(1.2, 1)
        circuit.cz(0, 1)
        circuit.rx(0.5, 0)
        circuit.ry(0.3, 1)
        if measure:
            circuit.measure_all()
        return circuit

    def test_qiskit_measure_all(self):
        import Qdislib as qd

        exact = _exact(self._qiskit(False), "ZZ")
        for cut in (["CZ_3"], []):
            value = qd.gate_cutting(
                self._qiskit(True), cut, observable="ZZ", shots=40000
            )
            self.assertLess(abs(value - exact), 0.05, cut)

    def test_qibo_final_m(self):
        import Qdislib as qd
        from qibo import models, gates

        def build(measure):
            circuit = models.Circuit(2)
            circuit.add(gates.RY(0, 0.7))
            circuit.add(gates.RY(1, 1.2))
            circuit.add(gates.CZ(0, 1))
            circuit.add(gates.RX(0, 0.5))
            if measure:
                circuit.add(gates.M(0, 1))
            return circuit

        exact = qd.analytical_solution(build(False), "ZZ")
        value = qd.gate_cutting(
            build(True), ["CZ_3"], observable="ZZ", backend="qibo", shots=40000
        )
        self.assertLess(abs(value - exact), 0.05)
        # M(0, 1) must not look like a two-qubit gate to find_cut.
        self.assertEqual(qd.find_cut(build(True), wire_cut=False, seed=0), ["CZ_3"])

    def test_estimate(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(3)
        circuit.h(0)
        circuit.cz(0, 1)
        circuit.cz(1, 2)
        circuit.measure_all()
        value = qd.estimate(circuit, "ZZZ", max_qubits=2, backend="qiskit")
        self.assertLess(abs(value - 0.0), 0.1)

    def test_mid_circuit_measurement_is_kept(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2, 1)
        circuit.h(0)
        circuit.measure(0, 0)
        circuit.h(0)
        circuit.cz(0, 1)
        gates = [
            data["gate"] for _, data in qd.convert_to_dag(circuit).nodes(data=True)
        ]
        self.assertIn("measure", gates)


@unittest.skipUnless(QIBO and QISKIT, "qibo and qiskit required")
class RepeatedGateCutsTest(BaseTimedTestCase):
    """bsc-wdc/qdislib#3: several gate cuts on one qubit keep every sign."""

    @staticmethod
    def _two_cz():
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.ry(0.7, 0)
        circuit.rx(1.1, 1)
        circuit.cz(0, 1)
        circuit.ry(0.9, 0)
        circuit.rx(0.6, 1)
        circuit.cz(0, 1)
        circuit.ry(0.4, 0)
        circuit.rx(0.8, 1)
        return circuit

    def test_two_cz_cuts_on_one_pair(self):
        import Qdislib as qd

        for observable in ("ZZ", "IZ"):
            value = qd.gate_cutting(
                self._two_cz(),
                ["CZ_3", "CZ_6"],
                observable=observable,
                backend="qiskit",
                shots=40000,
            )
            exact = _exact(self._two_cz(), observable)
            self.assertLess(abs(value - exact), 0.08, observable)

    def test_general_gate_decomposed_into_several_cz(self):
        import Qdislib as qd
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.ry(0.7, 0)
        circuit.rx(1.1, 1)
        circuit.rxx(0.5, 0, 1)  # two CZs once decomposed
        circuit.ry(0.4, 0)
        circuit.rx(0.8, 1)
        value = qd.gate_cutting(
            circuit, ["RXX_3"], observable="ZZ", backend="qiskit", shots=40000
        )
        self.assertLess(abs(value - _exact(circuit, "ZZ")), 0.08)

    def test_fold_cut_outcomes(self):
        from Qdislib.core.cutting_algorithms.wire_cutting import _fold_cut_outcomes

        for outcomes in ([0], [1], [0, 0], [0, 1], [1, 1], [1, 0, 1]):
            product = 1
            for m in outcomes:
                product *= -((-1) ** m)
            self.assertEqual(-((-1) ** _fold_cut_outcomes(outcomes)), product)


@unittest.skipUnless(QISKIT, "qiskit not installed")
class DeprecatedQpuPathTest(BaseTimedTestCase):
    def test_a_pauli_sum_keeps_the_qpu_config(self):
        # each term connected to IBM from the environment instead
        from unittest import mock

        from qiskit import QuantumCircuit

        import Qdislib as qd
        from Qdislib.core.cutting_algorithms import gate_cutting as gc

        qc = QuantumCircuit(2)
        qc.cz(0, 1)
        config = {"token": "t", "channel": "ibm_cloud"}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with (
                mock.patch.object(
                    gc, "_check_ibm_qc", return_value=(None, None)
                ) as connect,
                mock.patch.object(gc, "_execute_gate_cutting", return_value=[1.0]),
            ):
                qd.gate_cutting(
                    qc,
                    ["CZ_1"],
                    [(0.5, "ZZ"), (0.5, "ZI")],
                    qpu=True,
                    qpu_dict={"IBM_Quantum": 5},
                    qpu_config=config,
                )
        self.assertEqual(connect.call_count, 2)
        for call in connect.call_args_list:
            self.assertEqual(call.args[0], config)


class GeneralGateAfterWireCutTest(BaseTimedTestCase):
    """A general gate cut that follows a wire cut on the same qubit
    (bsc-wdc/qdislib#17): the CY decomposes into a CZ on the wire the cut moved."""

    def _circuit(self):
        from qiskit import QuantumCircuit

        qc = QuantumCircuit(3)
        qc.ry(0.6, 0)
        qc.ry(0.4, 1)
        qc.ry(0.5, 2)
        qc.cz(0, 1)
        qc.cy(1, 2)
        qc.ry(0.3, 2)
        qc.rx(0.9, 1)
        return qc

    def _exact(self, qc):
        from qiskit.quantum_info import SparsePauliOp, Statevector

        return float(Statevector(qc).expectation_value(SparsePauliOp("ZZZ")).real)

    def test_each_cut_and_both_give_the_circuits_value(self):
        import Qdislib as qd

        qc = self._circuit()
        exact = self._exact(qc)
        self.assertAlmostEqual(exact, 0.1386, places=3)  # as the report has it
        wire = [("RY_2", "CZ_4")]
        for backend in ("qiskit", "qibo"):
            for cuts in (
                {"gate": ["CY_5"], "wire": wire},
                {"wire": wire},
                {"gate": ["CY_5"]},
                {"gate": ["CZ_4"], "wire": wire},
            ):
                value = qd.circuit_cutting(qc, cuts, "ZZZ", backend=backend, shots=None)
                self.assertAlmostEqual(value, exact, places=8, msg=f"{backend} {cuts}")

    def test_with_shots_and_the_deprecated_call(self):
        import Qdislib as qd

        qc = self._circuit()
        with self.assertWarns(DeprecationWarning):
            value = qd.gate_wire_cutting(
                qc,
                ["CY_5"],
                [("RY_2", "CZ_4")],
                observable="ZZZ",
                software="qiskit",
                shots=20000,
            )
        self.assertAlmostEqual(value, self._exact(qc), delta=0.1)


if __name__ == "__main__":
    unittest.main()
