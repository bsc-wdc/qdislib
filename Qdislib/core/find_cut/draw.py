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

"""draw_cut: the circuit as a diagram, coloured by the fragments a cut leaves."""

import networkx

__all__ = ["draw_cut"]

#: One colour per fragment (cycled beyond ten).
_PALETTE = (
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#B279A2",
    "#72B7B2",
    "#EECA3B",
    "#9D755D",
    "#FF9DA6",
    "#79706E",
    "#BAB0AC",
)
_GATE_CUT = "#D62728"  # measure-and-prepare gate cut
_HADAMARD_CUT = "#6A3D9A"  # Hadamard gate cut
_WIRE_CUT = "#D62728"
_INVALID = "#222222"
_IDLE = "#C8C8C8"


def _columns(dag):
    """ASAP time step of every gate; a multi-qubit gate also blocks the qubits
    its connector crosses, so nothing is drawn underneath it."""
    nxt, column = {}, {}
    for node in networkx.topological_sort(dag):
        qubits = dag.nodes[node].get("qubits") or ()
        if not qubits:
            continue
        span = range(min(qubits), max(qubits) + 1)
        col = max(nxt.get(q, 0) for q in span)
        column[node] = col
        for q in span:
            nxt[q] = col + 1
    return column


def _label(node, data):
    name = str(data.get("gate") or node).upper()
    return name if len(name) <= 4 else name[:4]


def _index(node):
    tail = str(node).rsplit("_", 1)[-1]
    return tail if tail.isdigit() else ""


def draw_cut(circuit, cuts, filename=None, *, hadamard=False, ax=None, show=None):
    """Draw a circuit with a cut: every fragment in its own colour, every cut marked.

    The circuit is drawn as a diagram -- qubit lines left to right, gates in time
    order, each with its index above it (gate ``CZ`` numbered 7 is ``"CZ_7"``, the
    name to cut). The qubit lines and gates take the colour of the fragment they
    end up in, so you see the independent pieces the cut produces. A
    measure-and-prepare gate cut is a red dashed gate with ``✂``, a Hadamard gate cut
    a purple one, a wire cut a red ``✂`` on the wire. A cut that does not separate
    its two sides is drawn in black with ``✗`` and named in the title (it would be
    refused by the cutting functions). The title gives the number of fragments,
    the cuts and the sampling overhead ``gamma``.

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: As in :func:`~Qdislib.circuit_cutting`: a list of gate names and
        ``(source, target)`` wire edges, or ``{"gate": [...], "hadamard": [...],
        "wire": [...]}`` (what :func:`~Qdislib.find_cut` returns with
        ``mixed=True``).
    :param filename: Save the figure there (any format matplotlib writes).
    :param hadamard: Draw the gate cuts of a list as Hadamard cuts.
    :param ax: A matplotlib ``Axes`` to draw into (a new figure otherwise).
    :param show: Show the figure; by default only when it is neither saved nor
        drawn into a given ``ax``.
    :return: The matplotlib ``Axes``.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch, Patch

    from Qdislib.utils.backends import convert_to_dag
    from Qdislib.core.cutting_algorithms.mixed_cutting import _CutPlan
    from Qdislib.core.cutting_algorithms.hadamard_cut import _num_circuit_qubits

    dag = convert_to_dag(circuit)
    if not isinstance(dag, networkx.DiGraph):
        raise TypeError("draw_cut requires a circuit convertible to a Qdislib DAG.")
    plan = _CutPlan(circuit, cuts, None, hadamard, False, True, validate=False)
    invalid = plan.invalid_cuts()

    # Fragment of every node of the cut circuit, and of each side of a cut gate.
    fragment_of = {}
    for index, fragment in enumerate(plan.fragments):
        for node in fragment["nodes"]:
            fragment_of[node] = index
    gate_kind, gate_sides = {}, {}
    for kind, names, marker in (
        ("gate", plan.gate_cuts, "GATE"),
        ("hadamard", plan.had_cuts, "SUBS"),
    ):
        for k, name in enumerate(names, start=1):
            gate_kind[name] = kind
            if name in dag:
                qubits = dag.nodes[name]["qubits"]
                gate_sides[name] = {
                    qubits[0]: fragment_of.get(f"{marker}1_{k}"),
                    qubits[1]: fragment_of.get(f"{marker}2_{k}"),
                }
    # A general gate cut by measure-and-prepare is cut as the CZs it decomposes
    # into, under other names: mark the gate the user named.
    for name in plan.requested_gate_cuts:
        gate_kind.setdefault(name, "gate")

    column = _columns(dag)
    width = max(
        _num_circuit_qubits(circuit),
        1
        + max(
            (q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()),
            default=0,
        ),
    )
    num_columns = 1 + max(column.values(), default=0)

    def side_fragment(node, qubit):
        """The fragment ``node`` belongs to on ``qubit``."""
        if node in gate_sides and gate_sides[node].get(qubit) is not None:
            return gate_sides[node][qubit]
        if node in fragment_of:
            return fragment_of[node]
        # a cut gate under another name: take its neighbour's on that qubit
        for neighbours in (dag.predecessors(node), dag.successors(node)):
            for other in neighbours:
                if (
                    qubit in (dag.nodes[other].get("qubits") or ())
                    and other in fragment_of
                ):
                    return fragment_of[other]
        return None

    def colour(fragment):
        return _IDLE if fragment is None else _PALETTE[fragment % len(_PALETTE)]

    created = ax is None
    if created:
        fig_w = min(40.0, max(6.0, 0.8 * num_columns + 2.2))
        fig_h = max(2.4, 0.66 * width + 1.7)
        _, ax = plt.subplots(figsize=(fig_w, fig_h))
        if show is None:
            show = filename is None

    # Wire cuts: where each one sits on its qubit line.
    wire_marks = {}  # qubit -> [(x, cut, invalid)]
    for cut in plan.wire_cuts:
        source, target, qubit = cut
        x = (column[source] + column[target]) / 2
        edge_invalid = any(
            cut == bad or (isinstance(bad, tuple) and bad[:2] == cut[:2])
            for bad in invalid
        )
        wire_marks.setdefault(qubit, []).append((x, cut, edge_invalid))

    # Qubit lines, segment by segment, in their fragments' colours.
    by_qubit = {q: [] for q in range(width)}
    for node, col in column.items():
        for q in dag.nodes[node]["qubits"]:
            by_qubit[q].append((col, node))
    left, right = -0.8, num_columns - 0.2
    for q in range(width):
        y = -q
        events = sorted(by_qubit[q])
        if not events:
            ax.plot([left, right], [y, y], color=_IDLE, lw=1.6, zorder=1)
            continue
        cuts_here = sorted(wire_marks.get(q, []))
        points = [(left, side_fragment(events[0][1], q))]
        for col, node in events:
            points.append((col, side_fragment(node, q)))
        points.append((right, points[-1][1]))
        for (x0, f0), (x1, f1) in zip(points, points[1:]):
            splits = [x for x, _, _ in cuts_here if x0 < x < x1]
            if splits:
                ax.plot([x0, splits[0]], [y, y], color=colour(f0), lw=2.4, zorder=1)
                ax.plot([splits[0], x1], [y, y], color=colour(f1), lw=2.4, zorder=1)
            else:
                ax.plot([x0, x1], [y, y], color=colour(f0), lw=2.4, zorder=1)
        ax.text(left - 0.15, y, f"q{q}", ha="right", va="center", fontsize=9)

    # Gates.
    box_w, box_h = 0.62, 0.5
    for node, col in column.items():
        data = dag.nodes[node]
        qubits = list(data["qubits"])
        kind = gate_kind.get(node)
        is_invalid = node in invalid
        edge = (
            _INVALID
            if is_invalid
            else {"gate": _GATE_CUT, "hadamard": _HADAMARD_CUT}.get(kind, "#333333")
        )
        if len(qubits) > 1:
            ax.plot(
                [col, col],
                [-min(qubits), -max(qubits)],
                color=edge,
                lw=2.0 if kind else 1.3,
                ls=(0, (3, 2)) if kind else "-",
                zorder=2,
            )
        for q in qubits:
            ax.add_patch(
                FancyBboxPatch(
                    (col - box_w / 2, -q - box_h / 2),
                    box_w,
                    box_h,
                    boxstyle="round,pad=0.02,rounding_size=0.08",
                    facecolor=colour(side_fragment(node, q)),
                    edgecolor=edge,
                    lw=2.2 if kind else 0.8,
                    ls=(0, (3, 2)) if kind else "-",
                    zorder=3,
                )
            )
            ax.text(
                col,
                -q,
                _label(node, data),
                ha="center",
                va="center",
                fontsize=7.5,
                color="white",
                fontweight="bold",
                zorder=4,
            )
        ax.text(
            col,
            -min(qubits) + box_h / 2 + 0.08,
            _index(node),
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#666666",
            zorder=4,
        )
        if kind:
            # Between the gate's top qubit and the next line, so the mark never
            # lands on a qubit line it does not cut.
            y_mark = -min(qubits) - 0.5 if len(qubits) > 1 else -qubits[0]
            ax.text(
                col + (0.0 if len(qubits) > 1 else box_w / 2 + 0.12),
                y_mark,
                "✗" if is_invalid else "✂",
                ha="center",
                va="center",
                fontsize=13,
                color=edge,
                zorder=5,
                bbox=dict(boxstyle="circle,pad=0.08", fc="white", ec="none"),
            )

    # Wire cuts.
    for q, marks in wire_marks.items():
        for x, _, bad in marks:
            mark = _INVALID if bad else _WIRE_CUT
            ax.plot(
                [x, x],
                [-q - 0.32, -q + 0.32],
                color=mark,
                lw=1.6,
                ls=(0, (2, 1.5)),
                zorder=5,
            )
            ax.text(
                x,
                -q,
                "✗" if bad else "✂",
                ha="center",
                va="center",
                fontsize=13,
                color=mark,
                zorder=6,
                bbox=dict(boxstyle="circle,pad=0.08", fc="white", ec="none"),
            )

    # Title and legend.
    counts = [
        (len(plan.requested_gate_cuts), "gate cut"),
        (len(plan.had_cuts), "Hadamard cut"),
        (len({c[:2] for c in plan.wire_cuts}), "wire cut"),
    ]
    parts = [f"{n} {what}{'s' if n != 1 else ''}" for n, what in counts if n]
    title = f"{len(plan.fragments)} fragment{'s' if len(plan.fragments) != 1 else ''}"
    title += " · " + (", ".join(parts) if parts else "no cuts")
    if parts:
        title += f" · γ = {plan.gamma():g}"
    if invalid:
        title += f"\nNot a valid cut (its two sides stay connected): {invalid}"
    ax.set_title(
        title, fontsize=11, loc="left", color="#8B0000" if invalid else "black"
    )

    handles = []
    for index, fragment in enumerate(plan.fragments):
        qubits = {
            q
            for node in fragment["nodes"]
            for q in plan.graph.nodes[node].get("qubits") or ()
        }
        handles.append(
            Patch(
                color=colour(index),
                label=f"fragment {index + 1}: {len(qubits)} qubit{'s' if len(qubits) != 1 else ''}",
            )
        )
    for present, colour_, label in (
        (plan.requested_gate_cuts, _GATE_CUT, "gate cut (measure-and-prepare)"),
        (plan.had_cuts, _HADAMARD_CUT, "gate cut (Hadamard)"),
        (plan.wire_cuts, _WIRE_CUT, "wire cut"),
    ):
        if present:
            handles.append(
                Line2D(
                    [],
                    [],
                    color=colour_,
                    ls=(0, (3, 2)),
                    lw=2,
                    marker="$✂$",
                    markersize=11,
                    label=label,
                )
            )
    if invalid:
        handles.append(
            Line2D(
                [],
                [],
                color=_INVALID,
                ls=(0, (3, 2)),
                lw=2,
                marker="$✗$",
                markersize=9,
                label="not a valid cut",
            )
        )
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        fontsize=8.5,
    )

    ax.set_xlim(left - 0.9, right + 0.3)
    ax.set_ylim(-(width - 1) - 0.8, 0.9)
    ax.set_aspect("equal")
    ax.axis("off")
    if filename:
        ax.figure.savefig(filename, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    elif created and filename:
        plt.close(ax.figure)
    return ax
