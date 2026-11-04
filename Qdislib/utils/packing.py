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

"""Packing small circuits side by side on a large QPU.

A cut circuit becomes many small subcircuits; a device has many more qubits than
each of them uses. :func:`pack` places several circuits on disjoint, well
calibrated groups of the device's qubits and puts them in one device-wide
circuit, so one job (and one set of shots) runs all of them at once.
:meth:`PackedCircuits.unpack` splits the results back into each circuit's
counts, keyed as if the circuit had run alone::

    subs = qd.circuit_cutting_sampling_subcircuits(qc, {"hadamard": cut})
    circuits = subs.circuits("qiskit", defer_measurements=True)
    for c in circuits:
        c.measure_all()
    packed = qd.pack(circuits, backend)               # a Qiskit BackendV2
    results = [SamplerV2(mode=backend).run([c], shots=4096).result()[0]
               for c in packed.circuits]              # one plain job each
    dist = subs.reconstruct_counts(packed.unpack(results), bit_order="q0_last")

Groups of qubits are kept apart by at least ``buffer`` idle qubits, which limits
crosstalk between neighbouring circuits; the circuits themselves are compiled to
the device's gates on their own group only.
"""

from __future__ import annotations

__all__ = ["pack", "PackedCircuits"]

#: Target operations that are not gates to compile to.
_NOT_GATES = {
    "measure",
    "reset",
    "delay",
    "barrier",
    "if_else",
    "while_loop",
    "for_loop",
    "switch_case",
    "box",
}


def _calibration(target):
    """``(edges, readout)``: the lowest two-qubit gate error of each coupled pair
    (``{frozenset(pair): error}``) and each qubit's readout error."""
    edges = {}
    for name in target.operation_names:
        if name in _NOT_GATES:
            continue
        op = target.operation_from_name(name)
        if getattr(op, "num_qubits", 0) != 2:
            continue
        for qargs, props in target[name].items():
            if qargs is None:
                continue
            error = (
                props.error if props is not None and props.error is not None else 0.0
            )
            pair = frozenset(qargs)
            edges[pair] = min(error, edges.get(pair, 1.0))
    readout = {}
    if "measure" in target.operation_names:
        for qargs, props in target["measure"].items():
            if qargs is not None:
                readout[qargs[0]] = (
                    props.error
                    if props is not None and props.error is not None
                    else 0.0
                )
    return edges, readout


def _regions(edges, readout, size, buffer, max_error, max_readout_error):
    """Disjoint connected groups of ``size`` qubits, best first: each grown from a
    seed qubit by adding the neighbour with the lowest coupler + readout error,
    then chosen greedily so that no two are closer than ``buffer`` idle qubits."""
    good = {
        pair: error
        for pair, error in edges.items()
        if error <= max_error
        and all(readout.get(q, 0.0) <= max_readout_error for q in pair)
    }
    neighbours = {}
    for pair in good:
        a, b = tuple(pair)
        neighbours.setdefault(a, set()).add(b)
        neighbours.setdefault(b, set()).add(a)

    def grow(seed):
        group, cost = {seed}, readout.get(seed, 0.0)
        while len(group) < size:
            best = None
            for q in group:
                for w in neighbours.get(q, ()):
                    if w in group:
                        continue
                    added = good[frozenset((q, w))] + readout.get(w, 0.0)
                    if best is None or added < best[0]:
                        best = (added, w)
            if best is None:
                return None
            cost += best[0]
            group.add(best[1])
        return cost, tuple(sorted(group))

    candidates = sorted({g for g in map(grow, neighbours) if g is not None})
    if size == 1:  # single qubits: every qubit with a good readout
        candidates = sorted(
            (readout.get(q, 0.0), (q,))
            for q in readout
            if readout.get(q, 0.0) <= max_readout_error
        )

    def near(group, distance):
        out, frontier = set(group), set(group)
        for _ in range(distance):
            frontier = {w for q in frontier for w in neighbours.get(q, ())} - out
            out |= frontier
        return out

    chosen, taken = [], set()
    for _, group in candidates:
        if taken & set(group):
            continue
        chosen.append(group)
        taken |= near(group, buffer)
    return chosen


def _place(circuit, group, target, good, basis, optimization_level, seeds):
    """``circuit`` compiled onto the qubits of ``group`` (region-local indices):
    the result with the fewest two-qubit gates, then the lowest total error."""
    from qiskit import transpile
    from qiskit.transpiler import CouplingMap

    local = {q: i for i, q in enumerate(group)}
    directed = []
    for name in basis:
        op = target.operation_from_name(name)
        if getattr(op, "num_qubits", 0) == 2:
            for qargs in target[name]:
                if qargs is not None and qargs[0] in local and qargs[1] in local:
                    directed.append((local[qargs[0]], local[qargs[1]]))
    cmap = CouplingMap(sorted(set(directed))) if directed else None
    best = None
    for seed in range(seeds):
        placed = transpile(
            circuit,
            coupling_map=cmap,
            basis_gates=basis,
            optimization_level=optimization_level,
            seed_transpiler=seed,
        )
        two, error = 0, 0.0
        for inst in placed.data:
            if inst.operation.num_qubits == 2 and inst.operation.name in basis:
                a, b = (group[placed.find_bit(q).index] for q in inst.qubits)
                two += 1
                error += good.get(frozenset((a, b)), 0.0)
        if best is None or (two, error) < best[0]:
            best = ((two, error), placed)
    return best[1], best[0][0]


class PackedCircuits:
    """The result of :func:`pack`: the device-wide ``circuits`` to run (one job
    each) and, for every input circuit, where it went.

    :ivar circuits: The packed circuits, on the device's qubits (ISA circuits:
        run them as they are, without transpiling again).
    :ivar slots: For each input circuit, in order: ``(packed circuit index,
        register name, device qubits)``.
    :ivar regions: The groups of device qubits used, best first.
    :ivar two_qubit_gates: For each input circuit, its two-qubit gates once
        compiled onto its group.
    """

    def __init__(self, circuits, slots, regions, two_qubit_gates, registers):
        self.circuits = circuits
        self.slots = slots
        self.regions = regions
        self.two_qubit_gates = two_qubit_gates
        self._registers = registers  # per input circuit: its own registers' sizes

    def __len__(self):
        return len(self.slots)

    def __repr__(self):
        return (
            f"PackedCircuits({len(self.slots)} circuits in {len(self.circuits)} "
            f"device circuits, {len(self.regions)} groups of qubits)"
        )

    def unpack(self, results):
        """Each input circuit's counts, from the results of :attr:`circuits` (in
        order, one per packed circuit): SamplerV2 pub results (``job.result()[0]``)
        or counts dicts (``result.get_counts(i)`` of ``backend.run`` / Aer). The
        counts are keyed as Qiskit keys the circuit run alone (classical bit 0
        rightmost; several registers separated by spaces, the last one first)."""
        results = list(results)
        if len(results) != len(self.circuits):
            raise ValueError(
                f"unpack takes one result per packed circuit ({len(self.circuits)}), "
                f"got {len(results)}"
            )
        out = []
        for (index, name, _), sizes in zip(self.slots, self._registers):
            flat = self._register_counts(results[index], index, name)
            out.append(_split_registers(flat, sizes))
        return out

    def _register_counts(self, result, index, name):
        """``{bits: count}`` of one register of one packed circuit's result."""
        data = getattr(result, "data", None)
        if data is not None:  # a SamplerV2 pub result
            return dict(getattr(data, name).get_counts())
        names = [r.name for r in self.circuits[index].cregs]
        field = len(names) - 1 - names.index(name)  # Qiskit: last register first
        counts = {}
        for key, n in result.items():
            bits = key.split()[field]
            counts[bits] = counts.get(bits, 0) + n
        return counts


def _split_registers(flat, sizes):
    """Counts over a circuit's classical bits as one register, keyed as Qiskit
    keys that circuit (its own registers, last first, space separated)."""
    if len(sizes) <= 1:
        return flat
    out = {}
    for bits, n in flat.items():
        little = bits[::-1]  # classical bit 0 first
        parts, start = [], 0
        for size in sizes:
            parts.append(little[start : start + size][::-1])
            start += size
        key = " ".join(reversed(parts))
        out[key] = out.get(key, 0) + n
    return out


def pack(
    circuits,
    backend,
    *,
    size=None,
    buffer=1,
    max_error=0.02,
    max_readout_error=0.05,
    optimization_level=3,
    seeds=8,
    per_circuit=None,
    regions=None,
):
    """Pack small Qiskit circuits side by side into device-wide circuits.

    Each circuit is compiled onto its own group of ``size`` connected device qubits
    (the best calibrated groups first, at least ``buffer`` idle qubits apart), and
    the circuits of each packed circuit run at once on disjoint groups. A circuit
    without classical bits is measured on all its qubits (bit ``i`` <- qubit
    ``i``); otherwise its own measurements are kept.

    :param circuits: Qiskit ``QuantumCircuit`` objects (e.g. a cut circuit's
        ``subs.circuits("qiskit", defer_measurements=True)``, measured).
    :param backend: A Qiskit ``BackendV2`` (or its ``Target``): its coupling map,
        gates and calibrated errors choose the groups.
    :param size: Qubits per group; default the widest circuit.
    :param buffer: Idle qubits between groups (0: only disjoint).
    :param max_error: Couplers with a larger two-qubit gate error are not used.
    :param max_readout_error: Qubits with a larger readout error are not used.
    :param optimization_level: Qiskit's, for compiling each circuit on its group.
    :param seeds: Transpiler seeds tried per circuit (the fewest two-qubit gates win).
    :param per_circuit: At most this many circuits per packed circuit (default: one
        per group).
    :param regions: The groups of device qubits to use, in order, instead of
        choosing them (e.g. another run's ``packed.regions``, to run a circuit
        alone on the qubits it had when packed). Each must be connected and hold
        the widest circuit.
    :return: A :class:`PackedCircuits`; run its ``circuits`` (one job each) and
        pass the results to its ``unpack``.
    """
    from qiskit import ClassicalRegister, QuantumCircuit
    from qiskit.circuit.library import get_standard_gate_name_mapping

    circuits = list(circuits)
    if not circuits:
        raise ValueError("pack takes at least one circuit")
    target = getattr(backend, "target", backend)
    width = max(c.num_qubits for c in circuits)
    size = width if size is None else int(size)
    if size < width:
        raise ValueError(
            f"size={size} is smaller than the widest circuit ({width} qubits)"
        )
    edges, readout = _calibration(target)
    if regions is None:
        regions = _regions(
            edges, readout, size, int(buffer), max_error, max_readout_error
        )
        if not regions:
            raise ValueError(
                f"no group of {size} connected qubits passes max_error={max_error} and "
                f"max_readout_error={max_readout_error} on this device"
            )
    else:
        regions = [tuple(int(q) for q in group) for group in regions]
        if not regions or any(len(group) < width for group in regions):
            raise ValueError(f"regions= must hold groups of at least {width} qubits")
    slots_per = (
        len(regions) if per_circuit is None else min(int(per_circuit), len(regions))
    )
    good = {p: e for p, e in edges.items() if e <= max_error}
    standard = set(get_standard_gate_name_mapping())
    basis = [n for n in target.operation_names if n in standard and n not in _NOT_GATES]

    placed, gates, registers = [], [], []
    for k, circuit in enumerate(circuits):
        if circuit.num_clbits == 0:
            circuit = circuit.copy()
            circuit.measure_all(add_bits=True)
        registers.append([r.size for r in circuit.cregs])
        group = regions[k % slots_per]
        local, two = _place(
            circuit, group, target, good, basis, optimization_level, seeds
        )
        placed.append((local, group))
        gates.append(two)

    device, slots = [], []
    for start in range(0, len(placed), slots_per):
        big = QuantumCircuit(target.num_qubits)
        for slot, (local, group) in enumerate(placed[start : start + slots_per]):
            name = f"p{slot}"
            register = ClassicalRegister(local.num_clbits, name)
            big.add_register(register)
            big.compose(
                local,
                qubits=[group[i] for i in range(local.num_qubits)],
                clbits=list(register),
                inplace=True,
            )
            slots.append((len(device), name, group))
        device.append(big)
    return PackedCircuits(device, slots, regions[:slots_per], gates, registers)
