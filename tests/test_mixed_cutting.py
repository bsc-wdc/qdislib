#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Mixed wire + gate cutting on one circuit (circuit knitting)."""

import importlib.util
import unittest

from tests import BaseTimedTestCase

QIBO_AVAILABLE = importlib.util.find_spec("qibo") is not None


class ConfigSignTest(BaseTimedTestCase):
    """The combined quasiprobability sign (no qibo needed)."""

    def test_matches_signed_product_convention(self):
        from Qdislib.core.cutting_algorithms.mixed_cutting import _config_sign
        from Qdislib.core.cutting_algorithms._reconstruction import signed_product

        # _config_sign must agree with the standalone signed_product sign rule
        # (gate cut: base 6, flip {3,5}; wire cut: base 8, flip {3,5,7}).
        for index in range(40):
            self.assertEqual(
                _config_sign(index, 6, {3, 5}), signed_product([1.0], index, 6, {3, 5})
            )
            self.assertEqual(
                _config_sign(index, 8, {3, 5, 7}),
                signed_product([1.0], index, 8, {3, 5, 7}),
            )


@unittest.skipUnless(QIBO_AVAILABLE, "qibo not installed")
class MixedCuttingTest(BaseTimedTestCase):
    """The product QPD reconstructs the same expectation as the uncut circuit."""

    def setUp(self):
        super().setUp()
        import qibo

        qibo.set_backend("numpy")

    def _circuit(self):
        from qibo import gates, models

        c = models.Circuit(3)
        c.add(gates.H(0))  # H_1
        c.add(gates.CZ(0, 1))  # CZ_2  -> gate cut
        c.add(gates.RY(1, 0.5))  # RY_3
        c.add(gates.CZ(1, 2))  # CZ_4
        c.add(gates.RY(2, 0.3))  # RY_5
        c.add(gates.RX(0, 0.4))  # RX_6
        return c

    def test_gate_and_wire_cut_matches_analytical(self):
        # Cut one CZ (gate cut) AND one wire on the same circuit; the tensored QPD
        # must reconstruct the analytical expectation (within 1024-shot noise).
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.mixed_cutting import gate_wire_cutting

        circuit = self._circuit()
        analytical = qd.analytical_solution(circuit, "IZZ")  # non-degenerate
        got = gate_wire_cutting(
            self._circuit(),
            gates_cut=["CZ_2"],
            wire_cut=[("RY_3", "CZ_4")],
            observable="IZZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.2)

    def test_reduces_to_gate_cut_only(self):
        # With no wire cuts it must equal plain gate cutting (both shot-based, 1024).
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.mixed_cutting import gate_wire_cutting

        circuit = self._circuit()
        analytical = qd.analytical_solution(circuit, "ZIZ")
        got = gate_wire_cutting(
            self._circuit(),
            gates_cut=["CZ_2"],
            wire_cut=[],
            observable="ZIZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.2)

    def test_reduces_to_wire_cut_only(self):
        # With no gate cuts it must equal plain wire cutting.
        import Qdislib as qd
        from Qdislib.core.cutting_algorithms.mixed_cutting import gate_wire_cutting

        circuit = self._circuit()
        analytical = qd.analytical_solution(circuit, "IZZ")
        got = gate_wire_cutting(
            self._circuit(),
            gates_cut=[],
            wire_cut=[("RY_3", "CZ_4")],
            observable="IZZ",
            software="qibo",
            shots=1024,
        )
        self.assertLess(abs(analytical - got), 0.2)

    def test_exposed_in_public_api(self):
        import Qdislib as qd

        self.assertTrue(hasattr(qd, "gate_wire_cutting"))


if __name__ == "__main__":
    unittest.main()
