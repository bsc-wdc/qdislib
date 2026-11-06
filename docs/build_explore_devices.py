#!/usr/bin/env python3
"""Write docs/explore_devices.json: published device layouts that do not come with Qiskit.

    pip install cirq-google "iqm-client[qiskit]"
    python docs/build_explore_devices.py

Layouts only, from the packages that publish them: nothing connects to a device, and no
calibration is read. docs/build_explore_data.py gives these devices equal placeholder errors,
and the page says so. A device with real numbers can add "couplers": [[a, b, error], ...]
and "readout": [error, ...] in the same file.
"""

import json
import pathlib

TARGET = pathlib.Path(__file__).resolve().parent / "explore_devices.json"


def on_lattice(edges, count):
    """Coordinates on a square lattice, with every coupler one step long."""
    near = {q: set() for q in range(count)}
    for a, b in edges:
        near[a].add(b)
        near[b].add(a)
    order, seen = [], set()
    start = max(near, key=lambda q: len(near[q]))
    queue = [start]
    while queue:
        q = queue.pop(0)
        if q not in seen:
            seen.add(q)
            order.append(q)
            queue += sorted(near[q])
    place, taken = {}, {}

    def solve(i):
        if i == len(order):
            return True
        q = order[i]
        placed = [place[n] for n in near[q] if n in place]
        if placed:
            x0, y0 = placed[0]
            options = [(x0 + 1, y0), (x0, y0 + 1), (x0 - 1, y0), (x0, y0 - 1)]
        else:
            options = [(0, 0)]
        for spot in options:
            if spot in taken:
                continue
            if any(abs(spot[0] - x) + abs(spot[1] - y) != 1 for x, y in placed):
                continue
            place[q], taken[spot] = spot, q
            if solve(i + 1):
                return True
            del place[q], taken[spot]
        return False

    if not solve(0):
        raise ValueError("no square-lattice drawing")
    low_x, low_y = min(x for x, _ in place.values()), min(y for _, y in place.values())
    return {str(q): [x - low_x, y - low_y] for q, (x, y) in sorted(place.items())}


def iqm(key, name, family, cls):
    from iqm.qiskit_iqm import fake_backends

    backend = getattr(fake_backends, cls)()
    edges = sorted({tuple(sorted(edge)) for edge in backend.coupling_map.get_edges()})
    return {
        "key": key, "name": name, "family": family, "qubits": backend.num_qubits,
        "source": f"iqm-client, {cls}", "edges": [list(e) for e in edges],
        "points": on_lattice(edges, backend.num_qubits),
    }


def sycamore():
    import cirq_google

    device = cirq_google.Sycamore
    qubits = sorted(device.metadata.qubit_set)
    index = {q: i for i, q in enumerate(qubits)}
    edges = sorted({tuple(sorted((index[a], index[b]))) for a, b in device.metadata.nx_graph.edges})
    return {
        "key": "sycamore", "name": "Google Sycamore", "family": "", "qubits": len(qubits),
        "source": "cirq-google, Sycamore", "edges": [list(e) for e in edges],
        "points": {str(index[q]): [q.col, q.row] for q in qubits},
    }


def main():
    devices = [
        iqm("adonis", "IQM Adonis", "", "IQMFakeAdonis"),
        iqm("apollo", "IQM Apollo", "", "IQMFakeApollo"),
        sycamore(),
        iqm("aphrodite", "IQM Aphrodite", "", "IQMFakeAphrodite"),
    ]
    TARGET.write_text(json.dumps({"devices": devices}) + "\n")
    print(TARGET.name, [(d["name"], d["qubits"], len(d["edges"])) for d in devices])


if __name__ == "__main__":
    main()
