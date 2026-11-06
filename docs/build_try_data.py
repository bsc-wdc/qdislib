#!/usr/bin/env python3
"""Generate the data of the in-browser demo on the "Try it" page.

    python docs/build_try_data.py          # writes docs/source/_static/try_data.js

The page runs no quantum code: it replays these runs. Every number comes from Qdislib
through the same functions the poster and the notebook use: cutting_overhead and the
*_subcircuits lists for the costs, analytical_solution and shots=None for the exact
values, and gate_cutting, gate_cutting_hadamard and wire_cutting (with their _sampling
forms) on Qiskit Aer for the shot runs, several seeds per point. draw_cut makes the figures.

The methods are set up as the poster counts them, each of k cuts between two fragments:
joint=False for the Hadamard cut and for the optimized wire cut, and a circuit with
complex gates, so that no wire cut is a golden one. The page also shows the fewer
circuits the defaults of 2.0 run, joint=True.
"""

import base64
import io
import json
import multiprocessing
import pathlib
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

import Qdislib as qd

TARGET = pathlib.Path(__file__).resolve().parent / "source" / "_static" / "try_data.js"
SHOTS = [256, 1024, 4096, 16384]
CUTS = [1, 2, 3]
SEEDS = 5
TABLE_CUTS = range(1, 11)
METHODS = {
    "a": {"label": "Gate cut, mid-circuit measures", "formula": "6^k x 2"},
    "b": {"label": "Gate cut, Hadamard test", "formula": "6^k x 2"},
    "c": {"label": "Wire cut", "formula": "8^k x 2"},
    "d": {"label": "Optimized wire cut", "formula": "3^k + 4^k"},
}


def circuit_with_cuts(k):
    """Two blocks of k + 1 qubits, joined by k CZ gates; complex gates, so nothing is golden.

    Returns the circuit, the k gates between the blocks, and the wire of each: the
    edge between that qubit's last gate in its block and the gate that leaves it.
    """
    qc = QuantumCircuit(2 * k + 2)
    for q in range(2 * k + 2):
        qc.rx(0.4 + 0.15 * q, q)
    last = {}
    for block in (0, k + 1):
        for q in range(block, block + k):
            qc.cz(q, q + 1)
            last[q] = last[q + 1] = len(qc.data)
    gates, wires = [], []
    for q in range(k):
        qc.cz(q, k + 1 + q)
        gates.append(f"CZ_{len(qc.data)}")
        wires.append((f"CZ_{last[q]}", f"CZ_{len(qc.data)}"))
    for q in range(2 * k + 2):
        qc.rx(0.3, q)
    return qc, gates, wires


def setup(method, k):
    """(cuts for cutting_overhead, options, expectation fn, sampling fn, list fn)."""
    qc, gates, wires = circuit_with_cuts(k)
    if method == "a":
        return qc, {"gate": gates}, {}, gates, qd.gate_cutting, qd.gate_cutting_sampling, qd.gate_cutting_subcircuits
    if method == "b":
        return qc, {"hadamard": gates}, {"joint": False}, gates, qd.gate_cutting_hadamard, qd.gate_cutting_hadamard_sampling, qd.gate_cutting_hadamard_subcircuits
    if method == "c":
        return qc, {"wire": wires}, {}, wires, qd.wire_cutting, qd.wire_cutting_sampling, qd.wire_cutting_subcircuits
    return qc, {"wire": wires}, {"optimized": True, "joint": False}, wires, qd.wire_cutting, qd.wire_cutting_sampling, qd.wire_cutting_subcircuits


def exact_distribution(qc):
    return {b[::-1]: float(p) for b, p in Statevector(qc).probabilities_dict().items()}


def tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def figure(qc, cuts, k):
    """The drawing as SVG, sharp at any size and on any screen. Its text stays text:
    glyph outlines would make the page's data ten times larger."""
    with plt.rc_context({"svg.fonttype": "none", "font.family": "sans-serif"}):
        ax = qd.draw_cut(qc, cuts, show=False)
        buffer = io.BytesIO()
        ax.figure.savefig(buffer, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(ax.figure)
    return "data:image/svg+xml;base64," + base64.b64encode(buffer.getvalue()).decode()


def one_run(task):
    method, k, shots, seed = task
    qc, _, options, cuts, expectation, sampling, _ = setup(method, k)
    observable = "Z" * qc.num_qubits
    value = float(expectation(qc, cuts, observable, shots=shots, seed=seed, **options))
    distribution = sampling(qc, cuts, shots=shots, seed=seed, **options)
    exact = exact_distribution(qc)
    keys = sorted(exact, key=lambda b: -exact[b])[:8]
    return task, value, tvd(exact, distribution), [round(float(distribution.get(b, 0.0)), 4) for b in keys]


def method_entry(method, k):
    qc, spec, options, cuts, expectation, sampling, subcircuits = setup(method, k)
    observable = "Z" * qc.num_qubits
    cost = qd.cutting_overhead(qc, spec, **options)
    subs = subcircuits(qc, cuts, observable, **options)
    assert len(subs) == cost["subcircuits"], (method, k, len(subs), cost)
    widths = {}
    for s in subs:
        widths[str(s).split("fragment ")[1].split(",")[0]] = s.num_qubits
    exact = exact_distribution(qc)
    keys = sorted(exact, key=lambda b: -exact[b])[:8]
    default = {
        "b": lambda: qd.cutting_overhead(qc, spec),
        "d": lambda: qd.cutting_overhead(qc, spec, optimized=True),
    }.get(method)
    return {
        "qubits": qc.num_qubits,
        "subcircuits": cost["subcircuits"],
        "default_subcircuits": default()["subcircuits"] if default else None,
        "default_gamma": round(default()["gamma"], 3) if default else None,
        "gamma": round(cost["gamma"], 3),
        "fragments": sorted(widths.values()),
        "uncut": round(qd.analytical_solution(qc, observable), 6),
        "reconstructed": float(expectation(qc, cuts, observable, shots=None, **options)),
        "uncut_full": float(qd.analytical_solution(qc, observable)),
        "tvd_exact": tvd(exact, sampling(qc, cuts, shots=None, **options)),
        "keys": keys,
        "exact_top": [round(exact[b], 4) for b in keys],
        "figure": figure(qc, spec, k),
    }


def table():
    rows = {m: [] for m in METHODS}
    for k in TABLE_CUTS:
        for m in METHODS:
            qc, spec, options, *_ = setup(m, k)
            rows[m].append(qd.cutting_overhead(qc, spec, **options)["subcircuits"])
    return rows


def redraw_figures():
    """Only the drawings again, keeping the stored runs: python build_try_data.py --figures"""
    text = TARGET.read_text()
    header, payload = text.split("window.QD_TRY = ", 1)
    data = json.loads(payload.rstrip().rstrip(";"))
    for m in METHODS:
        for k in CUTS:
            qc, spec, *_ = setup(m, k)
            data["runs"][m][str(k)]["figure"] = figure(qc, spec, k)
    TARGET.write_text(
        header + "window.QD_TRY = " + json.dumps(data, separators=(",", ":")) + ";\n"
    )
    print(f"wrote {TARGET.name}: {TARGET.stat().st_size // 1024} KB")


def main():
    start = time.time()
    data = {
        "version": qd.__version__,
        "shots": SHOTS,
        "seeds": SEEDS,
        "cuts": CUTS,
        "methods": METHODS,
        "table_cuts": list(TABLE_CUTS),
        "table": table(),
        "runs": {m: {str(k): method_entry(m, k) for k in CUTS} for m in METHODS},
    }
    tasks = [(m, k, s, seed + 1) for m in METHODS for k in CUTS for s in SHOTS for seed in range(SEEDS)]
    with multiprocessing.Pool() as pool:
        results = pool.imap_unordered(one_run, tasks, chunksize=2)
        for count, ((m, k, s, seed), value, distance, top) in enumerate(results, 1):
            entry = data["runs"][m][str(k)].setdefault("runs", {}).setdefault(str(s), [None] * SEEDS)
            entry[seed - 1] = [round(value, 5), round(distance, 5), top]
            if count % 40 == 0:
                print(f"{count}/{len(tasks)} runs, {time.time() - start:.0f} s", flush=True)
    TARGET.write_text(
        "// Generated by docs/build_try_data.py from real Qdislib runs. Do not edit.\n"
        "window.QD_TRY = " + json.dumps(data, separators=(",", ":")) + ";\n"
    )
    print(f"wrote {TARGET.name}: {TARGET.stat().st_size // 1024} KB, {time.time() - start:.0f} s")
    for m in METHODS:
        print(m, [(k, data["runs"][m][str(k)]["subcircuits"], data["runs"][m][str(k)]["default_subcircuits"],
                   data["runs"][m][str(k)]["gamma"], data["runs"][m][str(k)]["fragments"],
                   f"{abs(data['runs'][m][str(k)]['reconstructed'] - data['runs'][m][str(k)]['uncut_full']):.1e}",
                   f"{data['runs'][m][str(k)]['tvd_exact']:.1e}") for k in CUTS], data["table"][m][:4])


if __name__ == "__main__":
    if "--figures" in sys.argv:
        redraw_figures()
    else:
        main()
