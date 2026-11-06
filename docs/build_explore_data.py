#!/usr/bin/env python3
"""Generate the stored runs of the Explore page.

    OMP_NUM_THREADS=1 python docs/build_explore_data.py    # writes docs/source/_static/explore_data.js

Everything but the page's circuit playground replays these: real Qdislib, on Qiskit Aer. The
script runs one thing at a time. The SDK section needs CUDA-Q, PennyLane, QiliSDK and Qibo;
a section can be run on its own, keeping the others: build_explore_data.py sdk noise
"""

import json
import math
import pathlib
import random
import re
import subprocess
import sys
import time
import warnings

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

import Qdislib as qd

TARGET = pathlib.Path(__file__).resolve().parent / "source" / "_static" / "explore_data.js"
SDK_TARGET = TARGET.with_name("explore_sdk_data.js")
HERE = pathlib.Path(__file__).resolve().parent


def ops_of(qc):
    """A Qiskit circuit as the page's list of gates."""
    out = []
    for instruction in qc.data:
        op = {"name": instruction.operation.name, "qubits": [qc.find_bit(q).index for q in instruction.qubits]}
        if instruction.operation.params:
            op["theta"] = round(float(instruction.operation.params[0]), 4)
        out.append(op)
    return out


def engine_cuts(qc, spec):
    """{"gate", "hadamard", "wire"} as the page's cuts, by gate position."""
    index = lambda name: int(name.rsplit("_", 1)[1]) - 1
    cuts = []
    for kind in ("gate", "hadamard"):
        cuts += [{"kind": kind, "op": index(name)} for name in spec.get(kind) or []]
    ops = ops_of(qc)
    for source, target in spec.get("wire") or []:
        a, b = index(source), index(target)
        qubit = sorted(set(ops[a]["qubits"]) & set(ops[b]["qubits"]))[0]
        cuts.append({"kind": "wire", "qubit": qubit, "from": a, "to": b})
    return cuts


def exact_distribution(qc, qubits=None):
    out = {}
    for bits, p in Statevector(qc).probabilities_dict().items():
        bits = bits[::-1]
        key = bits if qubits is None else "".join(bits[q] for q in qubits)
        out[key] = out.get(key, 0.0) + float(p)
    return out


def tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


# --------------------------------------------------------------------- find_cut
def family_chain():
    qc = QuantumCircuit(10)
    for q in range(10):
        qc.ry(0.4 + 0.07 * q, q)
    for q in range(9):
        qc.cz(q, q + 1)
    return qc


def family_qaoa():
    qc = QuantumCircuit(8)
    qc.h(range(8))
    for q in range(8):
        qc.rzz(0.6, q, (q + 1) % 8)
    qc.rx(0.8, range(8))
    return qc


def family_blocks():
    qc = QuantumCircuit(10)
    qc.h(range(10))
    for start, stop in ((0, 5), (5, 10)):
        for q in range(start, stop - 1):
            qc.rzz(1.2, q, q + 1)
    qc.rzz(0.1, 4, 5)
    qc.rzz(0.1, 0, 9)
    qc.rx(0.4, range(10))
    return qc


def family_ladder():
    qc = QuantumCircuit(8)
    qc.ry(0.6, range(8))
    for q in range(0, 7, 2):
        qc.cz(q, q + 1)
    for q in range(1, 7, 2):
        qc.cz(q, q + 1)
    qc.ry(0.5, range(8))
    for q in range(0, 7, 2):
        qc.cz(q, q + 1)
    for q in range(1, 7, 2):
        qc.cz(q, q + 1)
    return qc


def family_ring():
    qc = QuantumCircuit(8)
    for q in range(8):
        qc.ry(0.5 + 0.05 * q, q)
    for q in range(8):
        qc.cz(q, (q + 1) % 8)
    qc.rx(0.3, range(8))
    return qc


def family_star():
    qc = QuantumCircuit(7)
    for q in range(7):
        qc.ry(0.4 + 0.08 * q, q)
    for q in range(1, 7):
        qc.cz(0, q)
    qc.rx(0.3, range(7))
    return qc


def family_random():
    rng = np.random.default_rng(12)
    qc = QuantumCircuit(8)
    for q in range(8):
        qc.ry(float(rng.uniform(0.2, 2.5)), q)
    for _ in range(14):
        a, b = (int(q) for q in rng.choice(8, 2, replace=False))
        if abs(a - b) > 3:
            b = a + 1 if a < 7 else a - 1
        qc.cz(a, b)
        qc.rx(float(rng.uniform(0.2, 2.5)), a)
    return qc


def find_cut_section():
    families = {
        "chain": ("A chain of CZ gates", family_chain),
        "ring": ("A ring of CZ gates", family_ring),
        "ladder": ("A ladder: two layers of CZ gates", family_ladder),
        "blocks": ("Two blocks, joined by weak gates", family_blocks),
        "qaoa": ("A QAOA ring of rzz(0.6) gates", family_qaoa),
        "star": ("A star: one qubit joined to all", family_star),
        "random": ("A random circuit, seed 12", family_random),
    }
    out = []
    for key, (label, make) in families.items():
        qc = make()
        entry = {"key": key, "label": label, "n": qc.num_qubits, "ops": ops_of(qc), "runs": {}}
        for kind, options in (
            ("gate", {"wire_cut": False}),
            ("wire", {"gate_cut": False}),
            ("mixed", {"mixed": True, "time_limit": 3}),
        ):
            for max_qubits in range(2, qc.num_qubits):
                found = qd.find_cut(qc, max_qubits=max_qubits, seed=0, **options)
                spec = found if isinstance(found, dict) else ({"gate": found} if kind == "gate" else {"wire": found})
                run = {"spec": {k: v for k, v in spec.items() if v}}
                if any(spec.values()):
                    cost = qd.cutting_overhead(qc, spec)
                    run.update(cost, cuts=engine_cuts(qc, spec))
                    if cost["configurations"] <= 10**7 and cost["fragments"] <= 8:
                        value = qd.circuit_cutting(qc, spec, "Z" * qc.num_qubits, shots=None)
                        run["difference"] = abs(float(value) - float(qd.analytical_solution(qc, "Z" * qc.num_qubits)))
                entry["runs"][f"{kind}:{max_qubits}"] = run
            print("find_cut", key, kind, flush=True)
        out.append(entry)
    return out


# ---------------------------------------------------------------------- planner
def between_blocks(k, gate="cz", angle=None, real=False):
    """Two connected blocks of k + 1 qubits with k gates between them, and the wires."""
    qc = QuantumCircuit(2 * k + 2)
    for q in range(2 * k + 2):
        (qc.ry if real else qc.rx)(0.4 + 0.1 * q, q)
    last = {}
    for block in (0, k + 1):
        for q in range(block, block + k):
            qc.cz(q, q + 1)
            last[q] = last[q + 1] = len(qc.data)
    names, wires = [], []
    for q in range(k):
        if gate == "rzz":
            qc.rzz(angle, q, k + 1 + q)
            names.append(f"RZZ_{len(qc.data)}")
        else:
            qc.cz(q, k + 1 + q)
            names.append(f"CZ_{len(qc.data)}")
        wires.append((f"CZ_{last[q]}", names[-1]))
    return qc, names, wires


def golden_chain(k):
    """k + 1 fragments in a line with a real state: every wire golden."""
    qc = QuantumCircuit(k + 2)
    qc.ry(0.5, range(k + 2))
    qc.cz(0, 1)
    previous, wires = len(qc.data), []
    for q in range(1, k + 1):
        qc.cz(q, q + 1)
        wires.append((f"CZ_{previous}", f"CZ_{len(qc.data)}"))
        previous = len(qc.data)
    return qc, wires


def planner_section():
    kinds = {
        "gate": ("Gate cut, CZ", "gate_cutting", lambda k: (between_blocks(k)[0], {"gate": between_blocks(k)[1]}), {}),
        "rzz": ("Gate cut, rzz(0.3)", "gate_cutting", lambda k: (between_blocks(k, "rzz", 0.3)[0], {"gate": between_blocks(k, "rzz", 0.3)[1]}), {}),
        "had": ("Hadamard cut, one ancilla per gate", "gate_cutting_hadamard, joint=False", lambda k: (between_blocks(k)[0], {"hadamard": between_blocks(k)[1]}), {"joint": False}),
        "joint": ("Hadamard cut, joint", "gate_cutting_hadamard", lambda k: (between_blocks(k)[0], {"hadamard": between_blocks(k)[1]}), {}),
        "wire": ("Wire cut", "wire_cutting", lambda k: (between_blocks(k)[0], {"wire": between_blocks(k)[2]}), {}),
        "wireopt": ("Wire cut, optimized", "wire_cutting, optimized=True", lambda k: (between_blocks(k)[0], {"wire": between_blocks(k)[2]}), {"optimized": True}),
        "golden": ("Wire cut, real state in a chain", "wire_cutting", lambda k: (golden_chain(k)[0], {"wire": golden_chain(k)[1]}), {}),
    }
    out = {"cuts": list(range(1, 9)), "epsilon": 0.01, "delta": 0.05, "kinds": {}}
    for key, (label, call, make, options) in kinds.items():
        rows = []
        for k in out["cuts"]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                qc, spec = make(k)
                cost = qd.cutting_overhead(qc, spec, **options)
                shots = qd.shots_for(qc, spec, out["epsilon"], out["delta"], **options)
                # the page scales shots with 1 / epsilon**2: check that it may
                other = qd.shots_for(qc, spec, 0.1, out["delta"], **options)
            assert abs(shots / other - 100) < 1e-3 * 100 + 1, (key, k, shots, other)
            rows.append({"subcircuits": cost["subcircuits"], "gamma": round(cost["gamma"], 4), "fragments": cost["fragments"], "shots": shots})
        out["kinds"][key] = {"label": label, "call": call, "rows": rows}
        print("planner", key, flush=True)
    return out


# --------------------------------------------------------- what lowers the cost
def lower_section():
    observable = lambda qc: "Z" * qc.num_qubits
    cuts = list(range(1, 7))
    largest = [0.0]

    def values(make, options, exact_call):
        rows = []
        for k in cuts:
            qc, spec = make(k)
            cost = qd.cutting_overhead(qc, spec, **options)
            rows.append({"subcircuits": cost["subcircuits"], "gamma": round(cost["gamma"], 4)})
            if k <= 3:
                value = exact_call(qc, spec)
                largest[0] = max(largest[0], abs(float(value) - float(qd.analytical_solution(qc, observable(qc)))))
        return rows

    def run(**options):
        return lambda qc, spec: qd.circuit_cutting(qc, spec, observable(qc), shots=None, **options)

    gate = lambda k: (between_blocks(k)[0], {"gate": between_blocks(k)[1]})
    hadamard = lambda k: (between_blocks(k)[0], {"hadamard": between_blocks(k)[1]})
    wire = lambda k: (between_blocks(k)[0], {"wire": between_blocks(k)[2]})
    golden = lambda k: (golden_chain(k)[0], {"wire": golden_chain(k)[1]})
    groups = [
        {"label": "Gate cuts", "note": "optimized=True runs each fragment once per setting of the cuts that touch it, and reuses it in every term: a sum over fragments. gamma does not change.", "rows": [
            {"label": "As they are", "values": values(gate, {}, run())},
            {"label": "optimized=True", "values": values(gate, {"optimized": True}, run(optimized=True))},
        ]},
        {"label": "Hadamard cuts", "note": "With joint=True, the default, the gates crossing one boundary share one ancilla per side: gamma 2^(k+1) - 1 instead of 3^k.", "rows": [
            {"label": "joint=False", "values": values(hadamard, {"joint": False}, run(joint=False))},
            {"label": "joint=True", "values": values(hadamard, {}, run())},
        ]},
        {"label": "Wire cuts", "note": "With optimized=True a fragment runs once per setting of its wires, and with joint=True, the default, the wires it sends are measured together: 2^k + 1 circuits instead of 3^k.", "rows": [
            {"label": "As they are", "values": values(wire, {}, run())},
            {"label": "optimized, joint=False", "values": values(wire, {"optimized": True, "joint": False}, run(optimized=True, joint=False))},
            {"label": "optimized, joint=True", "values": values(wire, {"optimized": True}, run(optimized=True))},
        ]},
        {"label": "Golden wire cuts", "note": "When the fragment that sends a wire alone has a real state, the wire's Y terms are zero and are left out: 6 terms instead of 8, gamma 3 instead of 4. It is exact and automatic.", "rows": [
            {"label": "A real state", "values": values(golden, {}, run())},
            {"label": "optimized=True", "values": values(golden, {"optimized": True}, run(optimized=True))},
        ]},
    ]
    angles = []
    for phi in (0.05, 0.1, 0.3, 0.6, 1.0, 1.5708):
        make = lambda k: (between_blocks(k, "rzz", phi)[0], between_blocks(k, "rzz", phi)[1])
        row = {"phi": phi, "gamma": [], "joint": []}
        for k in cuts:
            qc, names = make(k)
            row["gamma"].append(round(qd.cutting_overhead(qc, {"gate": names})["gamma"], 4))
            row["joint"].append(round(qd.cutting_overhead(qc, {"hadamard": names})["gamma"], 4))
            if k <= 2:
                value = qd.circuit_cutting(qc, {"gate": names}, observable(qc), shots=None)
                largest[0] = max(largest[0], abs(float(value) - float(qd.analytical_solution(qc, observable(qc)))))
        angles.append(row)
    print("lower: largest difference", largest[0], flush=True)
    return {"cuts": cuts, "groups": groups, "angles": angles, "largest_difference": largest[0]}


# ---------------------------------------------------------------- distributions
def dist_circuit():
    qc = QuantumCircuit(8)
    for q, a in enumerate([0.9, 0.5, 1.2, 0.7, 0.4, 1.1, 0.6, 0.8]):
        qc.ry(a, q)
    for q in range(7):
        qc.cz(q, q + 1)
    for q, a in enumerate([0.3, 0.8, 0.2, 0.6, 0.9, 0.4, 0.7, 0.5]):
        qc.ry(a, q)
    return qc


def big_circuit(n=40):
    rng = np.random.default_rng(3)
    qc = QuantumCircuit(n)
    for q in range(n):
        qc.ry(float(rng.uniform(0.2, 1.2)), q)
    for q in range(n - 1):
        qc.cz(q, q + 1)
    for q in range(n):
        qc.ry(float(rng.uniform(0.2, 1.2)), q)
    return qc


SMALL_SHOTS = [256, 1024, 4096]


def distribution_section():
    qc = dist_circuit()
    cut = {"gate": ["CZ_12"]}
    exact = exact_distribution(qc)
    keys = sorted(exact, key=lambda b: -exact[b])[:12]
    shots_levels = [256, 1024, 4096, 16384]
    full = {"keys": keys, "exact": [round(exact[k], 4) for k in keys], "shots": shots_levels, "runs": []}
    for shots in shots_levels:
        runs = []
        for seed in (1, 2, 3):
            d = qd.circuit_cutting_sampling(qc, cut, shots=shots, seed=seed)
            runs.append({"top": [round(float(d.get(k, 0.0)), 4) for k in keys], "tvd": round(tvd(exact, d), 4)})
        full["runs"].append(runs)
    print("distributions: full", flush=True)
    marginals = {}
    order = ["00", "01", "10", "11"]
    for a in range(8):
        for b in range(a + 1, 8):
            e = exact_distribution(qc, [a, b])
            runs = []
            for shots in SMALL_SHOTS:
                d = qd.circuit_cutting_sampling(qc, cut, shots=shots, seed=1, qubits=[a, b])
                runs.append([round(float(d.get(k, 0.0)), 4) for k in order])
            marginals[f"{a},{b}"] = {"runs": runs, "exact": [round(e.get(k, 0.0), 4) for k in order]}
    print("distributions: marginals", flush=True)
    big = big_circuit()
    big_cut = {"gate": ["CZ_53", "CZ_66"]}
    cost = qd.cutting_overhead(big, big_cut)
    start = time.time()
    d = qd.circuit_cutting_sampling(big, big_cut, shots=None)
    elapsed = time.time() - start
    exact_pairs = {q: d.marginal([q, q + 1]) for q in range(39)}
    pairs = {f"{q},{q + 1}": [round(float(exact_pairs[q].get(k, 0.0)), 4) for k in order] for q in range(39)}
    sampled = {}
    for shots in SMALL_SHOTS:
        run = qd.circuit_cutting_sampling(big, big_cut, shots=shots, seed=1)
        sampled[str(shots)] = {f"{q},{q + 1}": [round(float(run.marginal([q, q + 1]).get(k, 0.0)), 4) for k in order] for q in range(39)}
    head = QuantumCircuit(10)
    for instruction in big.data:
        indices = [big.find_bit(q).index for q in instruction.qubits]
        if max(indices) < 10:
            head.append(instruction.operation, indices)
    check = max(
        abs(float(exact_pairs[q].get(k, 0.0)) - exact_distribution(head, [q, q + 1]).get(k, 0.0))
        for q in range(2, 7) for k in ("00", "01", "10", "11")
    )
    factored = {
        "text": repr(d),
        "seconds": round(elapsed, 1),
        "gamma": cost["gamma"],
        "subcircuits": cost["subcircuits"],
        "most_likely": {k: round(v, 8) for k, v in d.most_likely(5).items()},
        "pairs": pairs,
        "sampled": sampled,
        "check": check,
    }
    print("distributions: factored", factored["text"], "check", check, flush=True)
    return {"cut": "CZ_12", "n": 8, "ops": ops_of(qc), "cuts": engine_cuts(qc, cut), "full": full, "small_shots": SMALL_SHOTS, "marginals": marginals, "factored": factored}


# ------------------------------------------------------------------- error bars
def error_section():
    qc = QuantumCircuit(6)
    for q, a in enumerate([0.8, 0.5, 0.6, 0.4, 0.3, 0.7]):
        qc.ry(a, q)
    for q in range(5):
        qc.cz(q, q + 1)
    for q, a in enumerate([0.2, 0.1, 0.7, 0.3, 0.9, 0.4]):
        qc.rx(a, q)
    observable = "Z" * 6
    exact = float(qd.analytical_solution(qc, observable))
    cases = {
        "gate": ("One gate cut", lambda **k: qd.gate_cutting(qc, ["CZ_9"], observable, **k)),
        "wire": ("One wire cut", lambda **k: qd.wire_cutting(qc, [("CZ_8", "CZ_9")], observable, **k)),
        "two": ("Two gate cuts", lambda **k: qd.gate_cutting(qc, ["CZ_8", "CZ_10"], observable, **k)),
    }
    shots_levels = [256, 1024, 4096, 16384]
    out = {"exact": round(exact, 5), "shots": shots_levels, "cases": {}}
    for key, (label, run) in cases.items():
        rows = []
        for shots in shots_levels:
            values = [round(float(run(shots=shots, seed=100 + i)), 4) for i in range(24)]
            reported = run(shots=shots, seed=100, error=True)
            rows.append({"values": values, "reported": round(reported.error, 5), "value": round(float(reported), 5), "total": reported.shots, "spread": round(float(np.std(values, ddof=1)), 5)})
        out["cases"][key] = {"label": label, "rows": rows}
        print("error", key, flush=True)
    return out


# ------------------------------------------------------------------------ noise
NOISE_LEVELS = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05]


def noise_model(p):
    """On device gates, so the backend compiles every circuit to them and nothing runs clean:
    depolarizing p on cz, p / 10 on sx and x, and a readout error of p / 2. rz has no error."""
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error

    noise = NoiseModel(basis_gates=["cz", "sx", "x", "rz"])
    if p:
        noise.add_all_qubit_quantum_error(depolarizing_error(p, 2), ["cz"])
        noise.add_all_qubit_quantum_error(depolarizing_error(p / 10, 1), ["sx", "x"])
        noise.add_all_qubit_readout_error(ReadoutError([[1 - p / 2, p / 2], [p / 2, 1 - p / 2]]))
    return noise


def noise_circuit():
    """Two blocks of 3 qubits joined by one CZ, with the block on the right all after it: the
    same place to cut for all three kinds, a gate cut of that CZ or the wire just after it."""
    qc = QuantumCircuit(6)
    for q, a in enumerate([0.3, 0.2, 0.25, 0.15, 0.3, 0.2]):
        qc.ry(a, q)
    qc.cz(0, 1)
    qc.cz(1, 2)
    qc.cz(2, 3)
    middle = len(qc.data)
    qc.ry(0.25, 3)
    after = len(qc.data)
    qc.cz(3, 4)
    qc.cz(4, 5)
    qc.ry(0.25, range(6))
    qc.cz(0, 1)
    qc.cz(1, 2)
    qc.cz(3, 4)
    qc.cz(4, 5)
    return qc, f"CZ_{middle}", (f"CZ_{middle}", f"RY_{after}")


def noise_section():
    from qiskit import transpile
    from qiskit_aer import AerSimulator

    qc, gate, wire = noise_circuit()
    observable = "Z" * 6
    exact = float(qd.analytical_solution(qc, observable))
    # all three kinds, noiseless and exact: the cut is the right one
    for value in (
        qd.gate_cutting(qc, [gate], observable, shots=None),
        qd.gate_cutting_hadamard(qc, [gate], observable, shots=None),
        qd.wire_cutting(qc, [wire], observable, shots=None),
    ):
        assert abs(float(value) - exact) < 1e-9
    shots = 8000
    rows = []
    for p in NOISE_LEVELS:
        noise = noise_model(p)
        backend = qd.build_backend("AER_BACKEND", qd.BackendConfig(extra={"noise_model": noise}))
        simulator = AerSimulator(noise_model=noise, seed_simulator=1)
        measured = qc.copy()
        measured.measure_all()
        counts = simulator.run(transpile(measured, backend=simulator, optimization_level=0), shots=shots).result().get_counts()
        uncut = sum((-1) ** key.count("1") * c for key, c in counts.items()) / shots
        row = {"p": p, "uncut": [round(uncut, 4), round(((1 - uncut**2) / shots) ** 0.5, 4)]}
        for key, run in (
            ("gate", lambda **k: qd.gate_cutting(qc, [gate], observable, **k)),
            ("hadamard", lambda **k: qd.gate_cutting_hadamard(qc, [gate], observable, **k)),
            ("wire", lambda **k: qd.wire_cutting(qc, [wire], observable, **k)),
        ):
            value = run(backend=backend, shots=shots, seed=1, error=True)
            row[key] = [round(float(value), 4), round(value.error, 4)]
        rows.append(row)
        print("noise", row, flush=True)
    return {"exact": round(exact, 4), "shots": shots, "rows": rows, "gate": gate, "wire": list(wire), "n": 6, "ops": ops_of(qc), "cuts": engine_cuts(qc, {"gate": [gate]})}


# -------------------------------------------------------------------------- SDKs
SDK_SNIPPETS = {
    "Qiskit": ("qdislib[qiskit]", '''
from qiskit import QuantumCircuit
import Qdislib as qd

circuit = QuantumCircuit(4)
circuit.ry(0.8, 0)
circuit.cz(0, 1)
circuit.ry(0.5, 1)
circuit.cz(1, 2)
circuit.ry(0.6, 2)
circuit.cz(2, 3)

cut = qd.find_cut(circuit, max_qubits=2, wire_cut=False)
value = qd.gate_cutting(circuit, cut, backend="qiskit", shots=8192, seed=1)
print(cut, value, qd.analytical_solution(circuit, "ZZZZ"))
'''),
    "Qibo": ("qdislib[qibo]", '''
from qibo import gates, models
import Qdislib as qd

circuit = models.Circuit(4)
circuit.add(gates.RY(0, 0.8))
circuit.add(gates.CZ(0, 1))
circuit.add(gates.RY(1, 0.5))
circuit.add(gates.CZ(1, 2))
circuit.add(gates.RY(2, 0.6))
circuit.add(gates.CZ(2, 3))

cut = qd.find_cut(circuit, max_qubits=2, wire_cut=False)
value = qd.gate_cutting(circuit, cut, backend="qibo", shots=8192, seed=1)
print(cut, value, qd.analytical_solution(circuit, "ZZZZ"))
'''),
    "CUDA-Q": ("qdislib[cudaq]", '''
import cudaq
import Qdislib as qd

@cudaq.kernel
def circuit():
    q = cudaq.qvector(4)
    ry(0.8, q[0])
    z.ctrl(q[0], q[1])
    ry(0.5, q[1])
    z.ctrl(q[1], q[2])
    ry(0.6, q[2])
    z.ctrl(q[2], q[3])

cut = qd.find_cut(circuit, max_qubits=2, wire_cut=False)
value = qd.gate_cutting(circuit, cut, backend="cudaq", shots=8192, seed=1)
print(cut, value, qd.analytical_solution(circuit, "ZZZZ"))
'''),
    "PennyLane": ("qdislib[pennylane]", '''
import pennylane as qml
import Qdislib as qd

with qml.tape.QuantumTape() as circuit:
    qml.RY(0.8, wires=0)
    qml.CZ(wires=[0, 1])
    qml.RY(0.5, wires=1)
    qml.CZ(wires=[1, 2])
    qml.RY(0.6, wires=2)
    qml.CZ(wires=[2, 3])

cut = qd.find_cut(circuit, max_qubits=2, wire_cut=False)
value = qd.gate_cutting(circuit, cut, backend="pennylane", shots=8192, seed=1)
print(cut, value, qd.analytical_solution(circuit, "ZZZZ"))
'''),
    "QiliSDK": ("qdislib[qilisdk]", '''
from qilisdk.digital import CZ, RY, Circuit
import Qdislib as qd

circuit = Circuit(4)
circuit.add(RY(0, theta=0.8))
circuit.add(CZ(0, 1))
circuit.add(RY(1, theta=0.5))
circuit.add(CZ(1, 2))
circuit.add(RY(2, theta=0.6))
circuit.add(CZ(2, 3))

cut = qd.find_cut(circuit, max_qubits=2, wire_cut=False)
value = qd.gate_cutting(circuit, cut, backend="qilisdk", shots=8192, seed=1)
print(cut, value, qd.analytical_solution(circuit, "ZZZZ"))
'''),
}


def sdk_section():
    """Each snippet runs as written, in its own process."""
    out = []
    for name, (extra, code) in SDK_SNIPPETS.items():
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code], capture_output=True, text=True, timeout=170,
        )
        lines = [line for line in result.stdout.splitlines() if line.startswith("[")]
        if result.returncode or not lines:
            raise RuntimeError(f"{name}: {result.stderr[-600:]}")
        match = re.match(r"(\[.*\]) (\S+) (\S+)$", lines[-1])
        out.append({
            "name": name,
            "install": f'pip install "{extra}"',
            "code": code.strip(),
            "cut": match.group(1),
            "value": round(float(match.group(2)), 4),
            "exact": round(float(match.group(3)), 4),
        })
        print("sdk", out[-1]["name"], out[-1]["cut"], out[-1]["value"], flush=True)
    return out


# ------------------------------------------------------------------------- pack
#: The layouts on the page, in order of size. To add a device, add a row:
#: "fake": a backend qiskit_ibm_runtime ships, a real device's coupling map with a stored
#: calibration, read here with no connection to the device;
#: "read": a calibration read from IBM, in docs/explore_ibm_layouts.json;
#: "layout": a coupling map published by a package, in docs/explore_devices.json (made by
#: docs/build_explore_devices.py), with no calibration: it gets equal placeholder errors;
#: "drawn": a chip's published topology, written out below, also with placeholder errors.
PACK_DEVICES = [
    {"key": "manila", "name": "IBM Manila", "family": "", "fake": "FakeManilaV2"},
    {"key": "adonis", "name": "IQM Adonis", "family": "", "layout": "adonis"},
    {"key": "guadalupe", "name": "IBM Guadalupe", "family": "", "fake": "FakeGuadalupeV2"},
    {"key": "apollo", "name": "IQM Apollo", "family": "", "layout": "apollo"},
    {"key": "ona_blue", "name": "MareNostrum-Ona Blue", "family": "", "drawn": "ona_blue"},
    {"key": "montreal", "name": "IBM Montreal", "family": "Falcon", "fake": "FakeMontrealV2"},
    {"key": "sycamore", "name": "Google Sycamore", "family": "", "layout": "sycamore"},
    {"key": "aphrodite", "name": "IQM Aphrodite", "family": "", "layout": "aphrodite"},
    {"key": "ona_red", "name": "MareNostrum-Ona Red", "family": "", "drawn": "ona_red"},
    {"key": "sherbrooke", "name": "IBM Sherbrooke", "family": "Eagle", "fake": "FakeSherbrooke"},
    {"key": "torino", "name": "IBM Torino", "family": "Heron", "fake": "FakeTorino"},
    {"key": "kobe", "name": "IBM Kobe", "family": "Heron", "read": "ibm_kobe", "shape_of": "FakeFez"},
]
#: For a layout without a calibration: equal for every qubit and coupler, so nothing is
#: told apart by quality.
PLACEHOLDER = {"cz": 0.005, "sx": 0.0005, "readout": 0.01}


def ona_red():
    """MareNostrum-Ona Red, as its webpage draws it: a grid of (row, column), without its
    eight qubits that have no couplers."""
    grid = [
        [0, 1, 2, 3, 8, 9, 10, 11],
        [4, 5, 6, 7, 12, 13, 14, 15],
        [16, 17, 18, 19, 24, 25, 26, 27],
        [20, 21, 22, 23, 28, 29, 30, 31],
        [32, 33, 34, 35, 40, 41, 42, 43],
        [36, 37, 38, 39, 44, 45, 46, 47],
        [48, 49, 50, 51],
        [52, 53, 54, 55],
    ]
    across, down = [], []
    for r, row in enumerate(grid):
        for c in range(len(row) - 1):
            if r == 0 and c == 3:  # 3-8 is not connected
                continue
            across.append((row[c], row[c + 1]))
    for r in range(len(grid) - 1):
        for c in range(len(grid[r + 1])):
            down.append((grid[r][c], grid[r + 1][c]))
    assert len(across) == 47 and len(down) == 48
    edges = sorted(tuple(sorted(e)) for e in across + down)
    assert len(edges) == 95
    points = {str(q): [c, r] for r, row in enumerate(grid) for c, q in enumerate(row)}
    assert len(points) == 56
    return {"count": 56, "edges": edges, "points": points}


def ona_blue():
    """MareNostrum-Ona Blue, as its webpage draws it, numbered row by row for this page only."""
    spots = [(0, 2), (0, 3)] + [(r, c) for r in (1, 2, 3) for c in range(5)] + [(4, 1), (4, 2), (4, 3)]
    assert len(spots) == 20
    index = {spot: i for i, spot in enumerate(spots)}
    edges = sorted(
        tuple(sorted((index[(r, c)], index[(r + dr, c + dc)])))
        for r, c in spots for dr, dc in ((0, 1), (1, 0)) if (r + dr, c + dc) in index
    )
    assert len(edges) == 30
    points = {str(index[spot]): [spot[1], spot[0]] for spot in spots}
    return {"count": 20, "edges": edges, "points": points, "numbering": "The numbers are only for this page."}


def device_target(num_qubits, couplers, sx, readout):
    """A Qiskit Target from errors: qd.pack takes one. couplers maps (a, b) to an error."""
    from qiskit.circuit import Measure, Parameter
    from qiskit.circuit.library import CZGate, RZGate, SXGate, XGate
    from qiskit.transpiler import InstructionProperties, Target

    target = Target(num_qubits=num_qubits)
    target.add_instruction(CZGate(), {pair: InstructionProperties(error=e) for pair, e in couplers.items()})
    one = lambda errors: {(q,): InstructionProperties(error=errors[q]) for q in errors}
    target.add_instruction(SXGate(), one(sx))
    target.add_instruction(XGate(), one(sx))
    target.add_instruction(RZGate(Parameter("theta")), {(q,): InstructionProperties(error=0.0) for q in sx})
    target.add_instruction(Measure(), one(readout))
    return target


def check_drawing(key, points, edges):
    """Every coupler one step long and no two qubits in one place, or the drawing is wrong: a
    line on the page is then always a coupler of the device, and a coupler always a line."""
    spots = {tuple(p) for p in points.values()}
    assert len(spots) == len(points), (key, "two qubits in one place")
    for a, b in edges:
        (x0, y0), (x1, y1) = points[str(a)], points[str(b)]
        assert abs(x0 - x1) + abs(y0 - y1) == 1, (key, a, b, points[str(a)], points[str(b)])


def gate_map_coordinates(num_qubits):
    """The [row, column] positions Qiskit's own plot_gate_map uses, for the sizes it knows."""
    import ast
    import inspect

    from qiskit.visualization import gate_map

    for node in ast.walk(ast.parse(inspect.getsource(gate_map.plot_gate_map))):
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Subscript):
            target = node.targets[0]
            if getattr(target.value, "id", "") == "qubit_coordinates_map" and ast.literal_eval(target.slice) == num_qubits:
                return ast.literal_eval(node.value)
    return None


def heavy_hex_rows(edges, num_qubits):
    """Heron numbers its qubits along long rows, with the qubits between two rows numbered
    after the first: rows are runs of consecutive numbers, and a connector sits straight
    between its two neighbours, or hangs off the last row. check_drawing tests the result."""
    pairs = {tuple(sorted(edge)) for edge in edges}
    runs, run = [], [0]
    for q in range(1, num_qubits):
        if (q - 1, q) in pairs:
            run.append(q)
        else:
            runs.append(run)
            run = [q]
    runs.append(run)
    rows = [r for r in runs if len(r) >= 6]
    connectors = [q for r in runs if len(r) < 6 for q in r]
    near = {q: set() for q in range(num_qubits)}
    for a, b in pairs:
        near[a].add(b)
        near[b].add(a)
    place = {q: (column, 0) for column, q in enumerate(rows[0])}
    for number, row in enumerate(rows[1:], 1):
        shift = None
        for c in connectors:
            above = [n for n in near[c] if n in place and place[n][1] == 2 * (number - 1)]
            below = [n for n in near[c] if n in row]
            if above and below:
                shift = place[above[0]][0] - row.index(below[0])
                break
        for column, q in enumerate(row):
            place[q] = (column + shift, 2 * number)
    for c in connectors:
        ends = [n for n in near[c] if n in place]
        height = place[ends[0]][1] + 1 if len(ends) == 1 else (place[ends[0]][1] + place[ends[1]][1]) // 2
        place[c] = (place[ends[0]][0], height)
    return {str(q): list(place[q]) for q in range(num_qubits)}


def ibm_lattice(edges, num_qubits):
    """The chip as people know it: a line, Qiskit's own positions, or Heron's rows."""
    if all(b == a + 1 for a, b in edges) and len(edges) == num_qubits - 1:
        return {str(q): [q, 0] for q in range(num_qubits)}
    known = gate_map_coordinates(num_qubits)
    if known:
        return {str(q): [known[q][1], known[q][0]] for q in range(num_qubits)}
    return heavy_hex_rows(edges, num_qubits)


def placeholder_target(count, edges):
    couplers = {}
    for a, b in edges:
        couplers[(a, b)] = couplers[(b, a)] = PLACEHOLDER["cz"]
    return device_target(count, couplers, {q: PLACEHOLDER["sx"] for q in range(count)}, {q: PLACEHOLDER["readout"] for q in range(count)})


def load_device(entry):
    """(backend or target, facts for the page, qubits, undirected couplers, drawing)."""
    from qiskit_ibm_runtime import fake_provider

    facts = {"key": entry["key"], "name": entry["name"], "family": entry["family"]}
    here = pathlib.Path(__file__).resolve().parent
    if "fake" in entry:
        backend = getattr(fake_provider, entry["fake"])()
        edges = sorted({tuple(sorted(edge)) for edge in backend.coupling_map.get_edges()})
        facts.update(kind="snapshot", snapshot=backend.properties().last_update_date.date().isoformat())
        points = ibm_lattice(edges, backend.num_qubits)
        check_drawing(entry["key"], points, edges)
        return backend, facts, backend.num_qubits, edges, points
    if "read" in entry:
        device = next(d for d in json.loads((here / "explore_ibm_layouts.json").read_text())["devices"] if d["name"] == entry["read"])
        count = device["num_qubits"]
        couplers = {tuple(int(q) for q in pair.split(",")): value["error"] for pair, value in device["edges"].items()}
        target = device_target(
            count, couplers,
            {q: device["sx_error"][str(q)] for q in range(count)},
            {q: device["readout"][str(q)] for q in range(count)},
        )
        edges = sorted({tuple(sorted(pair)) for pair in couplers})
        # Drawn from the published shape, so the qubits with no coupler in this calibration
        # still have a place. The couplers it reads must all be on that shape.
        shape = getattr(fake_provider, entry["shape_of"])()
        whole = sorted({tuple(sorted(edge)) for edge in shape.coupling_map.get_edges()})
        assert set(edges) <= set(whole), "the couplers read are not all on the published shape"
        facts.update(kind="read", snapshot=device["calibrated"][:10], missing=len(whole) - len(edges))
        points = ibm_lattice(whole, count)
        check_drawing(entry["key"], points, whole)
        return target, facts, count, edges, points
    if "drawn" in entry:
        drawn = globals()[entry["drawn"]]()
        check_drawing(entry["key"], drawn["points"], drawn["edges"])
        facts.update(kind="layout", source="the MareNostrum-Ona webpage")
        if "numbering" in drawn:
            facts["numbering"] = drawn["numbering"]
        return placeholder_target(drawn["count"], drawn["edges"]), facts, drawn["count"], drawn["edges"], drawn["points"]
    device = next(d for d in json.loads((here / "explore_devices.json").read_text())["devices"] if d["key"] == entry["layout"])
    edges = [tuple(e) for e in device["edges"]]
    check_drawing(entry["key"], device["points"], edges)
    facts.update(kind="layout", source="the " + device["source"].split(",")[0] + " package")
    return placeholder_target(device["qubits"], edges), facts, device["qubits"], edges, device["points"]


def pack_section():
    qc = QuantumCircuit(6)
    for q, a in enumerate([0.8, 0.5, 0.6, 0.4, 0.3, 0.7]):
        qc.ry(a, q)
    for q in range(5):
        qc.cz(q, q + 1)
    qc.rx(0.3, range(6))
    subs = qd.circuit_cutting_subcircuits(qc, ["CZ_9"], "Z" * 6)
    circuits = []
    for c in subs.circuits("qiskit", defer_measurements=True):
        c = c.copy()
        c.measure_all()
        circuits.append(c)
    out = {"circuits": len(circuits), "widths": [c.num_qubits for c in circuits], "placeholder": PLACEHOLDER, "devices": []}
    for entry in PACK_DEVICES:
        key = entry["key"]
        start = time.time()
        backend, facts, count, edges, points = load_device(entry)
        device = dict(facts, qubits=count, edges=[list(edge) for edge in edges], points=points, runs={})
        for buffer in (0, 1, 2):
            try:
                packed = qd.pack(circuits, backend, buffer=buffer, max_error=0.1, max_readout_error=0.2, seeds=2, optimization_level=1)
            except ValueError as error:  # the device cannot hold them with this buffer
                device["runs"][str(buffer)] = {"refused": str(error)}
                print("pack", key, buffer, "refused:", error, flush=True)
                continue
            jobs = {}
            for index, (job, _, qubits) in enumerate(packed.slots):
                jobs.setdefault(job, []).append({"circuit": index, "qubits": [int(q) for q in qubits]})
            device["runs"][str(buffer)] = {
                "text": repr(packed),
                "jobs": [jobs[k] for k in sorted(jobs)],
                "regions": [[int(q) for q in region] for region in packed.regions],
            }
            print("pack", key, count, len(edges), buffer, repr(packed), round(time.time() - start), "s", flush=True)
        out["devices"].append(device)
    return out


SECTIONS = {
    "find_cut": find_cut_section,
    "planner": planner_section,
    "lower": lower_section,
    "distributions": distribution_section,
    "errors": error_section,
    "noise": noise_section,
    "sdk": sdk_section,
    "pack": pack_section,
}


def main():
    wanted = sys.argv[1:] or list(SECTIONS)
    data = {}
    if TARGET.exists():
        data = json.loads(TARGET.read_text().split("window.QD_EXPLORE = ", 1)[1].rsplit(";", 1)[0])
    data["version"] = qd.__version__
    for name in wanted:
        start = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = SECTIONS[name]()
        print(f"{name}: {time.time() - start:.0f} s", flush=True)
        if name == "sdk":  # on Try, which loads nothing else of this
            SDK_TARGET.write_text(
                "// Generated by docs/build_explore_data.py from real Qdislib runs. Do not edit.\n"
                "window.QD_EXPLORE = window.QD_EXPLORE || {};\n"
                "window.QD_EXPLORE.sdk = " + json.dumps(result, separators=(",", ":")) + ";\n"
            )
            continue
        data[name] = result
    data.pop("sdk", None)
    TARGET.write_text(
        "// Generated by docs/build_explore_data.py from real Qdislib runs. Do not edit.\n"
        "window.QD_EXPLORE = " + json.dumps(data, separators=(",", ":")) + ";\n"
    )
    print(f"wrote {TARGET.name}: {TARGET.stat().st_size // 1024} KB")


if __name__ == "__main__":
    main()
