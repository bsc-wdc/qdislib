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

"""A cut search that prices a cut the way Qdislib will execute it.

A cut is a partition of the circuit's *qubit timelines*: every two-qubit gate sits on
two (qubit, position) nodes, a gate whose two nodes land in different parts is cut,
and two consecutive nodes of one qubit in different parts are a wire cut. Every cut
then separates its two sides by construction, and a part's width is the number of
qubit lines it touches plus the Hadamard ancillas it carries.

The price is the sampling overhead gamma (shots grow as gamma**2), as a sum of logs:

- a gate cut costs the cheaper of measure-and-prepare (3 per CZ of the gate's
  decomposition) and the Hadamard test (the 1-norm of its KAK decomposition, so a
  weak rotation is cheap) -- after Fröhler et al., arXiv:2608.05287, and CiFold,
  arXiv:2412.18705;
- the Hadamard-cut gates crossing the same pair of parts can be cut jointly with one
  shared ancilla per side, ``2 * prod ||c_i||_1**2 - prod ||c_i||_2**2`` (``2**(k+1) - 1``
  for k CZ/CX; Schmitt, Piveteau & Sutter, Quantum 9, 1634 (2025)), priced per
  boundary as in Brandhofer et al., arXiv:2308.09567;
- a wire cut costs 4 (its joint measurement saves circuits, not gamma).

The search is anytime and bounded by ``time_limit``: it starts from the classic
partitioners and, for small circuits, an exact gate-cut MILP (SciPy / HiGHS), then
improves by local moves -- whole qubits first, then pieces of a qubit's timeline,
which is what a wire cut buys -- in the spirit of the two-stage search of Fröhler et
al.
"""

import itertools
import logging
import math
import random
import time

import networkx
import numpy as np

logger = logging.getLogger(__name__)

_LOG4 = math.log(4.0)
_INF = float("inf")


class _Gate:
    __slots__ = ("id", "a", "b", "ia", "ib", "mp", "had", "l1", "l2", "forbidden")


def _gate_costs(dag, node, cache):
    """``(log gamma measure-and-prepare, log gamma Hadamard, ||c||_1, ||c||_2**2)``."""
    from Qdislib.core.cutting_algorithms.hadamard_cut import (
        local_unitary_decomposition,
    )

    data = dag.nodes[node]
    name = str(data.get("gate") or "").upper()
    if name == "UNITARY":  # its one parameter is the gate's matrix
        params = tuple(data["parameters"])
        key = (name, np.asarray(params[0]).tobytes())
    else:
        params = tuple(float(p) for p in data.get("parameters") or ())
        key = (name, params)
    if key in cache:
        return cache[key]
    # measure-and-prepare, as circuit_cutting runs it: CZ/CX at gamma 3, any other
    # gate as its rzz rotations, each at its own angle (1 + 2|sin phi|)
    try:
        from Qdislib.core.cutting_algorithms.gate_cutting import _rotation_gamma

        mp = math.log(_rotation_gamma(name, params))
    except Exception:  # noqa: BLE001 -- unknown gate: assume the worst case
        mp = 3 * math.log(3.0)
    try:
        terms = local_unitary_decomposition(name, params)
        l1 = sum(abs(c) for c, _, _ in terms)
        l2 = sum(abs(c) ** 2 for c, _, _ in terms)
        had = math.log(max(1.0, 2 * l1**2 - l2))
    except Exception:  # noqa: BLE001 -- not cuttable with the Hadamard test
        l1, l2, had = _INF, _INF, _INF
    cache[key] = (mp, had, l1, l2)
    return cache[key]


class _Problem:
    """The circuit as qubit timelines, and what each cut on it costs."""

    def __init__(
        self,
        circuit,
        max_qubits,
        max_cuts,
        hadamard,
        joint,
        wires,
        gates,
        forbidden,
        gate_scale=None,
    ):
        from Qdislib.utils.backends import convert_to_dag, iter_gates

        self.dag = convert_to_dag(circuit)
        self.n = 1 + max(
            (q for _, d in self.dag.nodes(data=True) for q in d.get("qubits") or ()),
            default=0,
        )
        self.max_qubits = max_qubits
        self.max_cuts = max_cuts
        self.hadamard = hadamard
        self.joint = joint
        self.wires = wires
        self.allow_gates = gates
        forbidden_gates, forbidden_qubits = forbidden
        self.forbidden_qubits = set(forbidden_qubits or ())

        # Every multi-qubit gate in circuit order, and each qubit's sequence of them.
        self.seq = {q: [] for q in range(self.n)}
        self.gates, self.locks = [], []  # cuttable 2-qubit gates; >2-qubit gates
        cache = {}
        for index, (name, qubits) in enumerate(iter_gates(circuit), start=1):
            if len(qubits) < 2:
                continue
            gid = f"{name}_{index}"
            positions = []
            for q in qubits:
                positions.append((q, len(self.seq[q])))
                self.seq[q].append(gid)
            if len(qubits) > 2 or gid not in self.dag:
                self.locks.append(positions)
                continue
            g = _Gate()
            g.id = gid
            (g.a, g.ia), (g.b, g.ib) = positions
            g.mp, g.had, g.l1, g.l2 = _gate_costs(self.dag, gid, cache)
            if not hadamard:
                g.had = _INF
            if gate_scale and gid in gate_scale:
                g.mp *= gate_scale[gid]
                g.had *= gate_scale[gid]
            g.forbidden = (
                not gates
                or gid in (forbidden_gates or ())
                or bool(self.forbidden_qubits & {g.a, g.b})
            )
            self.gates.append(g)
        self.length = {q: max(1, len(s)) for q, s in self.seq.items()}

    # --------------------------------------------------------------- pricing
    def _joint(self, members):
        l1 = math.prod(g.l1**2 for g in members)
        l2 = math.prod(g.l2 for g in members)
        return math.log(max(1.0, 2 * l1 - l2))

    def _group_options(self, members):
        """``[(log gamma, ancillas per side, methods), ...]`` for the gates cut
        between one pair of fragments, cheapest first -- only what
        ``circuit_cutting`` can run: measure-and-prepare, and Hadamard cuts that
        (with ``joint``) share one ancilla per side across that pair."""
        options = [(sum(g.mp for g in members), 0, ["gate"] * len(members))]
        if not self.hadamard:
            return options
        cheaper = [g for g in members if g.had < g.mp]
        candidates = [cheaper]
        if self.joint and all(g.l1 < _INF for g in members):
            candidates.append(members)
        for chosen in candidates:
            if not chosen:
                continue
            ids = {id(g) for g in chosen}
            rest = sum(g.mp for g in members if id(g) not in ids)
            if self.joint:
                cost, anc = self._joint(chosen) + rest, 1
            else:
                cost, anc = sum(g.had for g in chosen) + rest, len(chosen)
            methods = ["hadamard" if id(g) in ids else "gate" for g in members]
            options.append((cost, anc, methods))
        return sorted(options, key=lambda option: option[0])

    def evaluate(self, part):
        """``(log gamma, solution)`` of an assignment ``part[q][i]``, or ``(inf,
        None)`` when it breaks a constraint.

        The fragments are the connected pieces of each part (union-find over the
        uncut timeline and gate edges), so cuts are grouped, priced and their
        ancillas counted per pair of fragments, exactly as execution does.
        """
        for lock in self.locks:  # a gate on 3+ qubits is never cut
            if len({part[q][i] for q, i in lock}) > 1:
                return _INF, None
        parent = {}

        def find(x):
            root = x
            while parent.get(root, root) != root:
                root = parent[root]
            while parent.get(x, x) != root:
                parent[x], x = root, parent[x]
            return root

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry

        wires = []
        for q in range(self.n):
            row = part[q]
            for i in range(len(row) - 1):
                if row[i] != row[i + 1]:
                    if not self.wires or q in self.forbidden_qubits:
                        return _INF, None
                    wires.append((q, i))
                else:
                    union((q, i), (q, i + 1))
        crossing = []
        for g in self.gates:
            if part[g.a][g.ia] != part[g.b][g.ib]:
                if g.forbidden:
                    return _INF, None
                crossing.append(g)
            else:
                union((g.a, g.ia), (g.b, g.ib))
        for lock in self.locks:
            for node in lock[1:]:
                union(lock[0], node)
        n_cuts = len(crossing) + len(wires)
        if self.max_cuts is not None and n_cuts > self.max_cuts:
            return _INF, None
        qubits = {}
        for q in range(self.n):
            for i in range(self.length[q]):
                qubits.setdefault(find((q, i)), set()).add(q)
        width = {f: len(qs) for f, qs in qubits.items()}
        groups = {}
        for g in crossing:
            fa, fb = find((g.a, g.ia)), find((g.b, g.ib))
            groups.setdefault((fa, fb) if fa <= fb else (fb, fa), []).append(g)
        choice = {key: self._group_options(m) for key, m in groups.items()}
        picked = dict.fromkeys(groups, 0)
        ancillas = dict.fromkeys(width, 0)
        for (fa, fb), options in choice.items():
            ancillas[fa] += options[0][1]
            ancillas[fb] += options[0][1]
        cap = self.max_qubits
        for f in sorted(width, key=lambda f: -(width[f] + ancillas[f])):
            if width[f] > cap:
                return _INF, None
            while width[f] + ancillas[f] > cap:
                # too wide with its ancillas: fall back, boundary by boundary, to
                # options with fewer ancillas (measure-and-prepare has none)
                best = None
                for key in choice:
                    if f not in key:
                        continue
                    current = choice[key][picked[key]]
                    for index, option in enumerate(choice[key]):
                        if option[1] < current[1]:
                            delta = option[0] - current[0]
                            if best is None or delta < best[0]:
                                best = (delta, key, index)
                if best is None:
                    return _INF, None
                _, key, index = best
                saved = choice[key][picked[key]][1] - choice[key][index][1]
                for side in key:
                    ancillas[side] -= saved
                picked[key] = index
        cost = sum(choice[key][picked[key]][0] for key in groups) + _LOG4 * len(wires)
        methods = {}
        for key, members in groups.items():
            for g, method in zip(members, choice[key][picked[key]][2]):
                methods[g.id] = method
        return cost, (methods, wires)

    # --------------------------------------------------------------- results
    def to_cuts(self, solution):
        """The ``{"gate", "hadamard", "wire"}`` dict ``circuit_cutting`` takes."""
        from Qdislib.core.find_cut.find_cut import _find_nodes_with_qubit

        methods, wires = solution
        cuts = {"gate": [], "hadamard": [], "wire": []}
        for gid, method in methods.items():
            cuts[method].append(gid)
        for q, i in wires:
            source = self.seq[q][i]
            nxt = _find_nodes_with_qubit(
                self.dag, source, qubit=q, direction="successor"
            )
            if not nxt:
                return None
            cuts["wire"].append((source, nxt[0]))
        return cuts

    def check(self, cuts):
        """The gamma of ``cuts`` as ``circuit_cutting`` will run them, or ``None``
        when a cut does not separate or a fragment is wider than ``max_qubits``.

        Fragments come from the cut plan (which checks that every cut separates);
        Hadamard cuts are grouped per pair of fragments, as execution groups them
        with ``joint=True``, and priced and counted in closed form -- a joint group
        of k gates has ``2**k`` product terms, which are never written out.
        """
        from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan

        try:
            plan = _CutPlan(self.circuit, cuts, None, False, False, False)
        except Exception as exc:  # noqa: BLE001 -- includes InvalidCutError
            logger.debug("find_cut: candidate %s rejected: %s", cuts, exc)
            return None
        by_id = {g.id: g for g in self.gates}
        fragment_of = {
            node: index
            for index, fragment in enumerate(plan.fragments)
            for node in fragment["nodes"]
        }
        groups, ancillas = {}, [0] * len(plan.fragments)
        for k, gid in enumerate(plan.had_cuts, start=1):
            a = fragment_of[f"SUBS1_{k}"]
            b = fragment_of[f"SUBS2_{k}"]
            key = (min(a, b), max(a, b)) if self.joint else ("gate", k)
            groups.setdefault(key, [a, b, []])[2].append(by_id[gid])
        # the measure-and-prepare cuts as the plan runs them (each rotation at its
        # angle)
        log_gamma = sum(
            math.log(sum(abs(w) for w in table[2])) for table in plan.gate_tables
        )
        # 4 per wire, 3 at a golden cut point (the plan finds them)
        log_gamma += sum(math.log(len(d) / 2) for d in plan.wire_digits)
        for a, b, members in groups.values():
            ancillas[a] += 1
            ancillas[b] += 1
            l1 = math.prod(g.l1**2 for g in members)
            l2 = math.prod(g.l2 for g in members)
            log_gamma += math.log(max(1.0, 2 * l1 - l2))
        if self.max_qubits is not None:
            for index, fragment in enumerate(plan.fragments):
                qubits = {
                    q
                    for node in fragment["nodes"]
                    for q in plan.graph.nodes[node].get("qubits") or ()
                }
                if len(qubits) + ancillas[index] > self.max_qubits:
                    return None
        return math.exp(log_gamma)


# ------------------------------------------------------------------- search
def _qubit_assignment(problem, groups):
    """``part[q][i]`` from groups of qubits (each qubit whole in one part)."""
    part = {}
    for label, qubits in enumerate(groups):
        for q in qubits:
            part[q] = [label] * problem.length[q]
    for q in range(problem.n):
        part.setdefault(q, [0] * problem.length[q])
    return part


def _starts(problem, parts, seed, exact, deadline):
    """Initial partitions: the classic partitioners, an exact MILP, region growing."""
    from Qdislib.core.find_cut import find_cut as fc

    graph = networkx.MultiGraph()
    graph.add_nodes_from(range(problem.n))
    for g in problem.gates:
        graph.add_edge(g.a, g.b, gate=g.id, weight=1)
    starts = []
    if exact:
        groups = _exact_gate_partition(problem, parts, deadline)
        if groups:
            starts.append(("exact", groups))
    for name, fn in (
        ("metis", lambda: fc._metis_partition(graph, parts, None, None)),
        ("metis+1", lambda: fc._metis_partition(graph, parts + 1, None, None)),
        (
            "spectral",
            lambda: fc._find_best_spectral_clustering(graph, parts, None, None, seed),
        ),
        (
            "kernighan-lin",
            lambda: fc._find_best_kernighan_lin_partition(graph, None, None, seed=seed),
        ),
    ):
        try:
            groups, _ = fn()
        except Exception:  # noqa: BLE001
            groups = []
        if groups:
            starts.append((name, [list(c) for c in groups]))
    rng = random.Random(seed)
    # Breadth-first region growing from random seeds, also one qubit short of the
    # limit, which leaves a fragment room for a Hadamard ancilla.
    rooms = (0, 1) if problem.hadamard else (0,)
    for k, room in itertools.product(range(3), rooms):
        order = list(range(problem.n))
        rng.shuffle(order)
        size = max(1, (problem.max_qubits or -(-problem.n // 2)) - room)
        assigned, groups = set(), []
        for root in order:
            if root in assigned:
                continue
            group, frontier = [], [root]
            while frontier and len(group) < size:
                q = frontier.pop(0)
                if q in assigned:
                    continue
                assigned.add(q)
                group.append(q)
                frontier.extend(sorted(set(graph.neighbors(q)) - assigned))
            groups.append(group)
        starts.append((f"grow-{k}-{room}", groups))
    return starts


def _exact_gate_partition(problem, parts, deadline):
    """Gate cuts only, each gate at its cheapest separate gamma: the optimal qubit
    partition into parts of at most ``max_qubits`` (MILP, SciPy/HiGHS)."""
    try:
        import numpy as np
        from scipy.optimize import Bounds, LinearConstraint, milp
    except ImportError:
        return None
    n, cap = problem.n, problem.max_qubits or -(-problem.n // 2)
    gates = problem.gates
    best = None
    for K in range(max(2, parts), max(2, parts) + 2):
        left = deadline - time.time()
        if left <= 0.2:
            break
        nx_, ng = n * K, len(gates)
        X = lambda q, k: q * K + k  # noqa: E731
        Y = lambda i: nx_ + i  # noqa: E731
        cost = np.zeros(nx_ + ng)
        upper = np.ones(nx_ + ng)
        for i, g in enumerate(gates):
            cost[Y(i)] = min(g.mp, g.had)
            if g.forbidden:
                upper[Y(i)] = 0
        rows, lo, hi = [], [], []

        def row(entries, low, high):
            r = np.zeros(nx_ + ng)
            for col, val in entries:
                r[col] = val
            rows.append(r)
            lo.append(low)
            hi.append(high)

        for q in range(n):
            row([(X(q, k), 1) for k in range(K)], 1, 1)
        for k in range(K):
            row([(X(q, k), 1) for q in range(n)], 0, cap)
        for i, g in enumerate(gates):
            for k in range(K):
                row([(X(g.a, k), 1), (X(g.b, k), -1), (Y(i), -1)], -np.inf, 0)
                row([(X(g.b, k), 1), (X(g.a, k), -1), (Y(i), -1)], -np.inf, 0)
        for lock in problem.locks:
            qubits = sorted({q for q, _ in lock})
            for q in qubits[1:]:
                for k in range(K):
                    row([(X(qubits[0], k), 1), (X(q, k), -1)], 0, 0)
        if problem.max_cuts is not None:
            row([(Y(i), 1) for i in range(ng)], 0, problem.max_cuts)
        row([(X(0, 0), 1)], 1, 1)  # symmetry
        res = milp(
            cost,
            constraints=LinearConstraint(np.array(rows), lo, hi),
            integrality=np.ones(nx_ + ng),
            bounds=Bounds(0, upper),
            options={"time_limit": max(0.1, min(left - 0.1, 0.6 * left))},
        )
        if res.x is None:
            continue
        if best is None or res.fun < best[0] - 1e-9:
            groups = [[q for q in range(n) if res.x[X(q, k)] > 0.5] for k in range(K)]
            best = (res.fun, [g for g in groups if g], res.status == 0)
    if best:
        logger.debug(
            "find_cut: exact gate partition, log gamma %.3f (%s)",
            best[0],
            "optimal" if best[2] else "time limit",
        )
    return best[1] if best else None


def _repair(problem, part, parts):
    """Move qubits out of parts wider than ``max_qubits`` (ancillas aside)."""
    if problem.max_qubits is None:
        return part
    labels = sorted(set(range(parts + 2)) | {part[q][0] for q in range(problem.n)})
    for _ in range(problem.n):
        width = {p: 0 for p in labels}
        for q in range(problem.n):
            width[part[q][0]] += 1
        over = [p for p in labels if width[p] > problem.max_qubits]
        if not over:
            return part
        p = over[0]
        target = min(labels, key=lambda r: width[r])
        q = next(q for q in range(problem.n) if part[q][0] == p)
        part[q] = [target] * problem.length[q]
    return part


def _local_search(problem, part, labels, deadline, rng, segments):
    """First-improvement moves of whole qubits (and, with ``segments``, of a qubit's
    timeline before or after a position -- a wire cut) until none helps."""
    cost, solution = problem.evaluate(part)
    improved = True
    while improved and time.time() < deadline:
        improved = False
        moves = []
        for q in range(problem.n):
            for r in labels:
                moves.append((q, None, r))
            if segments and problem.wires and q not in problem.forbidden_qubits:
                for i in range(1, problem.length[q]):
                    for r in labels:
                        moves.append((q, ("before", i), r))
                        moves.append((q, ("after", i), r))
        rng.shuffle(moves)
        for q, segment, r in moves:
            if time.time() >= deadline:
                break
            row = part[q]
            if segment is None:
                if all(p == r for p in row):
                    continue
                new_row = [r] * len(row)
            else:
                side, i = segment
                new_row = list(row)
                rng_ = range(0, i) if side == "before" else range(i, len(row))
                if all(row[j] == r for j in rng_):
                    continue
                for j in rng_:
                    new_row[j] = r
            part[q] = new_row
            new_cost, new_solution = problem.evaluate(part)
            if new_cost < cost - 1e-12:
                cost, solution, improved = new_cost, new_solution, True
            else:
                part[q] = row
    return cost, solution


def cost_aware_cut(
    circuit,
    max_qubits=None,
    max_cuts=None,
    *,
    gates=True,
    hadamard=True,
    joint=True,
    wires=True,
    exact="auto",
    time_limit=10.0,
    seed=None,
    forbidden_gates=None,
    forbidden_qubits=None,
    balanced=True,
    gate_scale=None,
):
    """The cheapest cut found within ``time_limit``, as ``{"gate", "hadamard",
    "wire"}``, every cut separating, every fragment within ``max_qubits`` (its
    Hadamard ancillas included). Empty when no cut is needed or none fits.

    ``max_qubits=None`` asks for halves (``balanced``) or for any split into two or
    more parts (``balanced=False``). ``gate_scale`` multiplies the log-cost of the
    named gates (SparseCut makes long-routing gates cheaper to cut).

    :return: ``(cuts, gamma)``.
    """
    start = time.time()
    deadline = start + float(time_limit)
    problem = _Problem(
        circuit,
        max_qubits,
        max_cuts,
        hadamard,
        joint,
        wires,
        gates,
        (forbidden_gates, forbidden_qubits),
        gate_scale,
    )
    problem.circuit = circuit
    empty = {"gate": [], "hadamard": [], "wire": []}
    if max_cuts == 0 or (max_qubits is not None and problem.n <= max_qubits):
        return empty, 1.0
    if problem.max_qubits is None:  # halves, or any split
        problem.max_qubits = -(-problem.n // 2) if balanced else max(1, problem.n - 1)
    parts = max(2, -(-problem.n // problem.max_qubits))
    labels = list(range(min(problem.n, parts + 2)))
    if exact == "auto":
        exact = problem.n <= 40 and len(problem.gates) <= 400
    rng = random.Random(seed)

    results = []  # (cost, name, solution, part)
    starts = _starts(problem, parts, seed, exact, start + 0.5 * time_limit)
    for name, groups in starts:
        if time.time() >= deadline:
            break
        part = _repair(problem, _qubit_assignment(problem, groups), parts)
        budget = min(
            deadline, time.time() + 0.6 * (deadline - time.time()) / max(1, len(starts))
        )
        cost, solution = _local_search(
            problem, part, labels, budget, rng, segments=False
        )
        if solution is not None:
            results.append(
                (cost, name, solution, {q: list(r) for q, r in part.items()})
            )
    results.sort(key=lambda item: item[0])
    if problem.wires and results and time.time() < deadline:
        cost, name, _, part = results[0]
        cost2, solution2 = _local_search(
            problem, part, labels, deadline, rng, segments=True
        )
        if solution2 is not None and cost2 < cost - 1e-12:
            results.insert(0, (cost2, name + "+wires", solution2, part))

    checked = []
    for cost, name, solution, _ in results:
        cuts = problem.to_cuts(solution)
        gamma = problem.check(cuts) if cuts is not None and any(cuts.values()) else None
        if gamma is not None:
            checked.append((gamma, name, cuts))
    if checked:
        gamma, name, cuts = min(checked, key=lambda item: item[0])
        logger.info(
            "find_cut: %s cut, gamma %.4g, %d gate / %d Hadamard / %d wire cuts, %.2f s",
            name,
            gamma,
            len(cuts["gate"]),
            len(cuts["hadamard"]),
            len(cuts["wire"]),
            time.time() - start,
        )
        return cuts, gamma
    logger.warning(
        "find_cut found no cut within max_qubits=%s, max_cuts=%s in %.1f s.",
        max_qubits,
        max_cuts,
        time.time() - start,
    )
    return empty, _INF
