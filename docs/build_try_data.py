#!/usr/bin/env python3
"""Generate the data of the in-browser demo on the "Try it" page.

    OMP_NUM_THREADS=1 python docs/build_try_data.py [processes]    # writes docs/source/_static/try_data.js

The page runs no quantum code: it replays these runs. Every number comes from Qdislib:
cutting_overhead and the *_subcircuits lists for the costs, analytical_solution and
shots=None for the exact values, and gate_cutting, gate_cutting_hadamard and wire_cutting
(with their _sampling forms) on Qiskit Aer for the shot runs, several seeds per point.
draw_cut makes the figures.

Each method runs as Qdislib does by default, with joint cuts on, and with joint=False where
that changes anything: the Hadamard cut and the optimized wire cut, from two cuts up. Those
are the forms the poster counts: 6^k x 2, 8^k x 2 and 3^k + 4^k. The circuit has complex
gates, as a real state would make a single wire a golden cut point, and the counts of the
form without joint cuts are then lower than those. The page gives that case as a note,
from cutting_overhead.

Shot runs are stored as they finish in docs/.try_runs.json, which is not committed: run
the script again to go on after a stop. Two processes take about 15 minutes.
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
SEEDS = 4
PROGRESS = pathlib.Path(__file__).resolve().parent / ".try_runs.json"
TABLE_CUTS = range(1, 11)
METHODS = {
    "a": {"label": "Gate cut, mid-circuit measures", "formula": {"off": "6^k x 2", "on": "6^k x 2"}},
    "b": {"label": "Gate cut, Hadamard test", "formula": {"off": "6^k x 2", "on": "2^k (2^(k+1) - 1) x 2"}},
    "c": {"label": "Wire cut", "formula": {"off": "8^k x 2", "on": "8^k x 2"}},
    "d": {"label": "Optimized wire cut", "formula": {"off": "3^k + 4^k", "on": "(2^k + 1) + 4^k"}},
}
#: The methods joint cuts change, and from how many cuts.
JOINT = {"b": 2, "d": 2}
STATES = {"on": True, "off": False}


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


def setup(method, k, joint=True):
    """(circuit, cuts for cutting_overhead, options, cuts, expectation, sampling, list)."""
    qc, gates, wires = circuit_with_cuts(k)
    if method == "a":
        return qc, {"gate": gates}, {}, gates, qd.gate_cutting, qd.gate_cutting_sampling, qd.gate_cutting_subcircuits
    if method == "b":
        return qc, {"hadamard": gates}, {"joint": joint}, gates, qd.gate_cutting_hadamard, qd.gate_cutting_hadamard_sampling, qd.gate_cutting_hadamard_subcircuits
    if method == "c":
        return qc, {"wire": wires}, {}, wires, qd.wire_cutting, qd.wire_cutting_sampling, qd.wire_cutting_subcircuits
    return qc, {"wire": wires}, {"optimized": True, "joint": joint}, wires, qd.wire_cutting, qd.wire_cutting_sampling, qd.wire_cutting_subcircuits


def states(method, k):
    """The forms that differ: "on" alone where joint cuts change nothing."""
    return ("on", "off") if method in JOINT and k >= JOINT[method] else ("on",)


def exact_distribution(qc):
    return {b[::-1]: float(p) for b, p in Statevector(qc).probabilities_dict().items()}


def tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def figure(qc, cuts, k):
    """The drawing as SVG, sharp at any size and on any screen. Its text stays text:
    glyph outlines would make the page's data ten times larger."""
    with plt.rc_context({"svg.fonttype": "none", "font.family": "sans-serif"}):
        ax = qd.draw_cut(qc, cuts, show=False)
        ax.set_title("")  # it gives the default's gamma: the page gives the state shown
        buffer = io.BytesIO()
        ax.figure.savefig(buffer, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(ax.figure)
    return "data:image/svg+xml;base64," + base64.b64encode(buffer.getvalue()).decode()


def one_run(task):
    method, k, state, shots, seed = task
    qc, _, options, cuts, expectation, sampling, _ = setup(method, k, STATES[state])
    observable = "Z" * qc.num_qubits
    value = float(expectation(qc, cuts, observable, shots=shots, seed=seed, **options))
    distribution = sampling(qc, cuts, shots=shots, seed=seed, **options)
    exact = exact_distribution(qc)
    keys = sorted(exact, key=lambda b: -exact[b])[:8]
    return task, value, tvd(exact, distribution), [round(float(distribution.get(b, 0.0)), 4) for b in keys]


def state_entry(method, k, state):
    qc, spec, options, cuts, expectation, sampling, subcircuits = setup(method, k, STATES[state])
    observable = "Z" * qc.num_qubits
    cost = qd.cutting_overhead(qc, spec, **options)
    subs = subcircuits(qc, cuts, observable, **options)
    assert len(subs) == cost["subcircuits"], (method, k, state, len(subs), cost)
    widths = {}
    for s in subs:
        widths[str(s).split("fragment ")[1].split(",")[0]] = s.num_qubits
    exact = exact_distribution(qc)
    return {
        "subcircuits": cost["subcircuits"],
        "gamma": round(cost["gamma"], 3),
        "fragments": sorted(widths.values()),
        "reconstructed": float(expectation(qc, cuts, observable, shots=None, **options)),
        "tvd_exact": tvd(exact, sampling(qc, cuts, shots=None, **options)),
        "runs": {},
    }


def method_entry(method, k):
    qc, spec, *_ = setup(method, k)
    observable = "Z" * qc.num_qubits
    exact = exact_distribution(qc)
    keys = sorted(exact, key=lambda b: -exact[b])[:8]
    return {
        "qubits": qc.num_qubits,
        "uncut": float(qd.analytical_solution(qc, observable)),
        "keys": keys,
        "exact_top": [round(exact[b], 4) for b in keys],
        "states": {state: state_entry(method, k, state) for state in states(method, k)},
    }


def table():
    rows = {state: {m: [] for m in METHODS} for state in STATES}
    for k in TABLE_CUTS:
        for m in METHODS:
            for state, joint in STATES.items():
                qc, spec, options, *_ = setup(m, k, joint)
                rows[state][m].append(qd.cutting_overhead(qc, spec, **options)["subcircuits"])
    return rows


def figures():
    """One drawing per cut: the optimized wire cut draws as the wire cut does."""
    out = {}
    for m in "abc":
        for k in CUTS:
            qc, spec, *_ = setup(m, k)
            out[f"{m}:{k}"] = figure(qc, spec, k)
    return out


def golden_note():
    """One wire leaving a real state: what the same cut costs without the rx gates."""
    qc = QuantumCircuit(4)
    qc.ry(0.5, range(4))
    qc.cz(0, 1)
    qc.cz(2, 3)
    wire = len(qc.data)
    qc.cz(0, 2)
    spec = {"wire": [(f"CZ_{wire - 1}", f"CZ_{wire + 1}")]}
    cost = qd.cutting_overhead(qc, spec)
    optimized = qd.cutting_overhead(qc, spec, optimized=True)
    return {"subcircuits": cost["subcircuits"], "optimized": optimized["subcircuits"], "gamma": cost["gamma"]}


def redraw_figures():
    """Only the drawings again, keeping the stored runs: python build_try_data.py --figures"""
    text = TARGET.read_text()
    header, payload = text.split("window.QD_TRY = ", 1)
    data = json.loads(payload.rstrip().rstrip(";"))
    data["figures"] = figures()
    TARGET.write_text(
        header + "window.QD_TRY = " + json.dumps(data, separators=(",", ":")) + ";\n"
    )
    print(f"wrote {TARGET.name}: {TARGET.stat().st_size // 1024} KB")


def main():
    start = time.time()
    processes = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    tasks = [
        (m, k, state, s, seed + 1)
        for m in METHODS for k in CUTS for state in states(m, k) for s in SHOTS for seed in range(SEEDS)
    ]
    done = json.loads(PROGRESS.read_text()) if PROGRESS.exists() else {}
    todo = [t for t in tasks if json.dumps(t) not in done]
    print(f"{len(todo)} of {len(tasks)} runs to do", flush=True)
    with multiprocessing.Pool(processes) as pool:
        for count, (task, value, distance, top) in enumerate(pool.imap_unordered(one_run, todo), 1):
            done[json.dumps(task)] = [round(value, 5), round(distance, 5), top]
            if count % 8 == 0:
                PROGRESS.write_text(json.dumps(done))
                print(f"{count}/{len(todo)} runs, {time.time() - start:.0f} s", flush=True)
    PROGRESS.write_text(json.dumps(done))
    data = {
        "version": qd.__version__,
        "shots": SHOTS,
        "seeds": SEEDS,
        "cuts": CUTS,
        "methods": METHODS,
        "joint_from": JOINT,
        "table_cuts": list(TABLE_CUTS),
        "table": table(),
        "golden": golden_note(),
        "figures": figures(),
        "runs": {m: {str(k): method_entry(m, k) for k in CUTS} for m in METHODS},
    }
    for task in tasks:
        m, k, state, s, seed = task
        entry = data["runs"][m][str(k)]["states"][state]["runs"].setdefault(str(s), [None] * SEEDS)
        entry[seed - 1] = done[json.dumps(task)]
    TARGET.write_text(
        "// Generated by docs/build_try_data.py from real Qdislib runs. Do not edit.\n"
        "window.QD_TRY = " + json.dumps(data, separators=(",", ":")) + ";\n"
    )
    print(f"wrote {TARGET.name}: {TARGET.stat().st_size // 1024} KB, {time.time() - start:.0f} s")
    for m in METHODS:
        for k in CUTS:
            for state, e in data["runs"][m][str(k)]["states"].items():
                print(m, k, state, e["subcircuits"], e["gamma"], e["fragments"],
                      f"{abs(e['reconstructed'] - data['runs'][m][str(k)]['uncut']):.1e} {e['tvd_exact']:.1e}")
    print({state: {m: data["table"][state][m][:4] for m in METHODS} for state in STATES}, data["golden"])


if __name__ == "__main__":
    if "--figures" in sys.argv:
        redraw_figures()
    else:
        main()
