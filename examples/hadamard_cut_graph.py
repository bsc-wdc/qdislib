#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Render the Hadamard-cut PyCOMPSs task graph as Graphviz DOT.

This introspects the *actual* fan-out that :func:`gate_cutting_hadamard` produces --
one ``_evaluate_hadamard_fragment`` task per (configuration, fragment), all feeding a
single ``_reduce_hadamard`` fan-in task -- and writes it to a ``.dot`` file, WITHOUT
needing the PyCOMPSs runtime. It is a faithful preview of the dependency graph
``runcompss --graph`` emits (see ``examples/compss_app.py`` for the real run)::

    python3 examples/hadamard_cut_graph.py                 # default 1-cut example
    python3 examples/hadamard_cut_graph.py --cuts 2        # richer 2-cut fan-out
    dot -Tpng hadamard_cut_task_graph.dot -o graph.png     # render (needs graphviz)
"""

import argparse
import sys

import qibo
from qibo import gates, models

from Qdislib.core.cutting_algorithms.hadamard_cut import _hadamard_fragment_specs


def _one_cut():
    c = models.Circuit(2)
    c.add(gates.H(0))
    c.add(gates.CZ(0, 1))
    c.add(gates.RY(0, 0.8))
    return c, ["CZ_2"], "ZZ"


def _two_cuts():
    c = models.Circuit(3)
    c.add(gates.RY(0, 0.7))
    c.add(gates.RY(1, 0.9))
    c.add(gates.RY(2, 1.1))
    c.add(gates.CZ(0, 1))
    c.add(gates.CZ(1, 2))
    return c, ["CZ_4", "CZ_5"], "ZZZ"


def build_dot(specs):
    """Graphviz DOT for the fan-out -> fan-in task graph of ``specs`` (from
    :func:`_hadamard_fragment_specs`)."""
    lines = [
        "digraph hadamard_cut {",
        "  rankdir=TB;",
        "  labelloc=t;",
        '  label="gate_cutting_hadamard task graph";',
        "  node [style=filled, fontname=Helvetica];",
        '  circuit [label="circuit + cut", shape=note, fillcolor="#e8e0ff"];',
        '  reduce [label="_reduce_hadamard\\n(weighted sum of products)", '
        'shape=box, fillcolor="#ffd9b3"];',
        '  result [label="<observable>", shape=doublecircle, fillcolor="#d9ffd9"];',
    ]
    frag_id = 0
    for cfg_idx, (weight, fragments) in enumerate(specs):
        lines.append(f"  subgraph cluster_cfg{cfg_idx} {{")
        lines.append(f'    label="config {cfg_idx}  (weight {weight:+.3f})";')
        lines.append('    style=dashed; color="#9090b0";')
        for f_in_cfg, (frag, z_qubits) in enumerate(fragments):
            frag_id += 1
            lines.append(
                f'    t{frag_id} [label="evaluate_fragment\\n'
                f'{frag.nqubits}q, Z on {z_qubits}", '
                f'shape=ellipse, fillcolor="#cfe8ff"];'
            )
        lines.append("  }")
    lines.append(f"  // {frag_id} independent fragment tasks -> 1 reduction")
    for t in range(1, frag_id + 1):
        lines.append(f"  circuit -> t{t};")
        lines.append(f"  t{t} -> reduce;")
    lines.append("  reduce -> result;")
    lines.append("}")
    return "\n".join(lines), frag_id


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cuts", type=int, choices=(1, 2), default=1)
    parser.add_argument("--out", default="hadamard_cut_task_graph.dot")
    args = parser.parse_args(argv)

    qibo.set_backend("numpy")
    circuit, gates_cut, observable = _one_cut() if args.cuts == 1 else _two_cuts()
    specs = _hadamard_fragment_specs(circuit, gates_cut, observable)
    dot, n_tasks = build_dot(specs)
    with open(args.out, "w") as fh:
        fh.write(dot + "\n")

    print(f"cuts={len(gates_cut)}  observable={observable}")
    print(f"configurations = {len(specs)}")
    print(f"fragment tasks  = {n_tasks}  (fan-out) -> 1 reduction task (fan-in)")
    print(f"wrote {args.out}  (render: dot -Tpng {args.out} -o graph.png)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
