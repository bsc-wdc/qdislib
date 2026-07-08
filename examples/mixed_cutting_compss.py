#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Mixed wire + gate cutting under PyCOMPSs -- a task-graph showcase.

Cuts one CZ (gate cut) AND one wire on the same circuit, so the fragments over the
product QPD (6 gate configs x 8 wire configs) all fan out as independent COMPSs tasks::

    runcompss -g --pythonpath=$(pwd) --python_interpreter=python3 \\
        examples/mixed_cutting_compss.py

Exits non-zero unless the reconstruction matches the analytical reference.
"""

import os
import sys

import qibo
from qibo import gates, models

import Qdislib as qd

TOLERANCE = 0.2
SOFTWARE = os.environ.get("QDISLIB_SOFTWARE", "qibo")
SHOTS = int(os.environ.get("QDISLIB_SHOTS", "1024"))


def build_circuit():
    c = models.Circuit(3)
    c.add(gates.H(0))  # H_1
    c.add(gates.CZ(0, 1))  # CZ_2  -> gate cut
    c.add(gates.RY(1, 0.5))  # RY_3
    c.add(gates.CZ(1, 2))  # CZ_4
    c.add(gates.RY(2, 0.3))  # RY_5
    return c


def main():
    qibo.set_backend("numpy")
    observable = "IZZ"
    analytical = qd.analytical_solution(build_circuit(), observable)
    print(f"Backend: {SOFTWARE}  shots={SHOTS}")
    print(f"Analytical <{observable}> = {analytical:+.4f}")

    got = qd.gate_wire_cutting(
        build_circuit(),
        gates_cut=["CZ_2"],
        wire_cut=[("RY_3", "CZ_4")],
        observable=observable,
        software=SOFTWARE,
        shots=SHOTS,
    )
    print(
        f"Mixed gate+wire cutting -> {got:+.4f} (|diff| = {abs(got - analytical):.4f})"
    )

    ok = abs(got - analytical) < TOLERANCE
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
