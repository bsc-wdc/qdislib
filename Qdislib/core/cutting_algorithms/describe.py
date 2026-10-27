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

"""Readable views of the three-step objects: a subcircuit list
(``subs.describe(values)``), one subcircuit (its ``repr``) and one result
(:class:`FragmentResult`).

A result's values are Z-correlators named by the qubits they read --
``Z[q0,q1,q2]`` is ``<Z0 Z1 Z2>`` of that subcircuit, with the gate cut's
measurement sign applied -- and qubits are named as in the circuit: ``q2'`` is a
wire cut's prepared line carrying qubit 2 on, ``a0`` a Hadamard cut's ancilla.
"""

import html

import numpy as np

__all__ = ["FragmentResult", "Description"]

#: The measure-and-prepare gate cut's operations on one side of a term.
_GATE_OP = {
    "remove": "identity",
    "rz": "RZ(-π/2)",
    "z": "Z, RZ(-π/2)",
    "measure": "measure, RZ(-π/2)",
    "rzpi": "RZ(-π)",
    "pz": "Z",
    "pmeasure": "measure",
    "rz+": "RZ(+π/2)",
    "rz-": "RZ(-π/2)",
}
_STATE = {"0": "|0⟩", "1": "|1⟩", "+": "|+⟩", "-": "|-⟩", "+i": "|+i⟩", "-i": "|-i⟩"}
#: A distribution's result has 2**width correlators: show only the single-qubit ones.
_MAX_SHOWN = 6


class FragmentResult(tuple):
    """What running one subcircuit returns: ``(corr, qubits)``.

    ``corr`` maps a Z mask to its correlator (bit ``j`` = Z on ``qubits[j]``), or
    is the full vector over every mask (distribution lists, exact runs); ``qubits``
    are the circuit qubits the bits refer to. It unpacks as the pair it is; its
    ``repr`` names the values by their qubits (:meth:`labelled`).
    """

    __slots__ = ()

    def __new__(cls, corr, qubits):
        return tuple.__new__(cls, (corr, qubits))

    def __reduce__(self):
        return (FragmentResult, (self[0], self[1]))

    @property
    def corr(self):
        return self[0]

    @property
    def qubits(self):
        return self[1]

    def labelled(self, names=None, limit=None):
        """``{"Z[q0,q1]": value, ...}``; for a full vector (a distribution) only
        ``<Z>`` of each qubit. ``names`` renames the qubits."""
        return _labelled(self[0], list(self[1]), names, limit)

    def __repr__(self):
        values = self.labelled(limit=_MAX_SHOWN)
        shown = ", ".join(f"{k}: {v:+.6g}" for k, v in values.items())
        return f"FragmentResult({shown}; {_size(self[0], self[1])})"


def _size(corr, qubits):
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

    if isinstance(corr, _JointCounts):
        shots = sum(corr.joint.values())
        return f"{len(corr.joint)} outcomes of {shots} shots, any correlator"
    if isinstance(corr, dict):
        return f"{len(corr)} correlator{'s' if len(corr) != 1 else ''}"
    return f"all 2**{len(qubits)} correlators"


def _mask_label(mask, names, letters=None):
    """``Z[q0,q1]`` for a mask; with ``letters`` (the Pauli each qubit was measured
    in) the qubits are grouped by it: ``X[q0] Z[q2]``."""
    picked = [j for j in range(len(names)) if (mask >> j) & 1]
    if not picked:
        return "I"
    letters = letters or ["Z"] * len(names)
    parts = {}
    for j in picked:
        parts.setdefault(letters[j], []).append(names[j])
    return " ".join(f"{p}[{','.join(qs)}]" for p, qs in parts.items())


def _labelled(corr, qubits, names=None, limit=None, letters=None):
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

    names = names or [f"q{q}" for q in qubits]
    if isinstance(corr, _JointCounts):  # a cached run: show the single-qubit ones
        corr = corr.correlators()
        masks = [1 << j for j in range(len(qubits))]
    elif isinstance(corr, dict):
        masks = sorted(corr)
    else:  # every correlator: a distribution's run
        corr = np.asarray(corr).reshape(-1)
        masks = [1 << j for j in range(len(qubits)) if (1 << j) < len(corr)]
    if limit is not None:
        masks = masks[:limit]
    # + 0.0 turns an exact -0.0 into 0.0
    return {_mask_label(m, names, letters): float(corr[m]) + 0.0 for m in masks}


def _qubit_names(template, qubits):
    """Circuit names of a subcircuit's qubits: ``q3``, ``q2'`` (a wire cut's
    prepared line carrying qubit 2), ``a0`` (a Hadamard ancilla)."""
    names = []
    for q in qubits:
        if q >= template.ancilla_base:
            names.append(f"a{q - template.ancilla_base}")
        elif q in template.prepared_to_orig:
            names.append(f"q{template.prepared_to_orig[q]}'")
        else:
            names.append(f"q{q}")
    return names


def _basis_letters(template, qubits, basis):
    """The Pauli each result qubit was measured in (its observable's ``X`` / ``Y``
    was rotated into ``Z``); a prepared wire line follows its qubit."""
    if not basis:
        return None
    letters = []
    for q in qubits:
        original = template.prepared_to_orig.get(q, q)
        letter = basis[original] if original < len(basis) else "Z"
        letters.append(letter if letter in "XY" else "Z")
    return letters


def _template_qubits(template):
    used = {
        q for _, d in template.graph.nodes(data=True) for q in d.get("qubits") or ()
    }
    ancillas = set(template.fragment["ancilla_of"].values())
    return sorted(used) + sorted(ancillas - used)


def cut_terms(circuit):
    """What one subcircuit runs of each cut it touches, in words."""
    plan = circuit.template  # it carries the cuts' names and term tables
    fragment = plan.fragment
    digits, config = circuit.spec[0]
    states, setting = circuit.spec[1], circuit.spec[2]
    parts = []
    for k, side in fragment["gate_sides"]:
        op = plan.gate_tables[k - 1][side - 1][digits[k - 1]]
        parts.append(f"{plan.gate_labels[k - 1]} side {side}: {_GATE_OP.get(op, op)}")
    for g in fragment["had_groups"]:
        names = ",".join(plan.had_cuts[k - 1] for k in plan.groups[g][0])
        s, t, parity = config[g]
        parts.append(f"Hadamard {names}: term ({s},{t}), parity {parity}")
    for i, state in zip(fragment["preps"], states):
        source, target, _ = plan.wire_cuts[i - 1]
        parts.append(f"wire {source}→{target}: prepare {_STATE.get(state, state)}")
    if fragment["measures"]:
        wires = ",".join(
            f"{plan.wire_cuts[i - 1][0]}→{plan.wire_cuts[i - 1][1]}"
            for i in fragment["measures"]
        )
        if setting and setting[0] == "joint":
            parts.append(f"wire {wires}: joint measurement {setting[1]}")
        else:
            parts.append(f"wire {wires}: measure {''.join(setting)}")
    return "; ".join(parts) or "no cut"


def rows(subs, values=None, start=0, group=None):
    """One dict per subcircuit of a ``CutCircuits`` list."""
    basis = getattr(subs, "basis", None)
    out = []
    for n, circuit in enumerate(subs):
        qubits = _template_qubits(circuit.template)
        row = {
            "#": start + n,
            "group": group,
            "fragment": circuit.template.f,
            "qubits": " ".join(_qubit_names(circuit.template, qubits)),
            "terms": cut_terms(circuit),
            "draws": circuit.draws,
        }
        if values is not None:
            corr, result_qubits = values[n]
            names = _qubit_names(circuit.template, list(result_qubits))
            letters = _basis_letters(circuit.template, result_qubits, basis)
            row["values"] = _labelled(corr, list(result_qubits), names, None, letters)
        out.append(row)
    for column, unused in (
        ("group", group is None),
        ("draws", all(row["draws"] == 1 for row in out)),  # only with num_samples
    ):
        if unused:
            for row in out:
                del row[column]
    return out


class Description:
    """A table of a subcircuit list: ``.rows`` (dicts, e.g. for pandas), printed
    as text and shown as HTML in notebooks."""

    def __init__(self, title, table, footer=""):
        self.title = title
        self.rows = table
        self.footer = footer

    def _columns(self):
        return list(dict.fromkeys(k for row in self.rows for k in row))

    @staticmethod
    def _cell(value):
        if isinstance(value, dict):
            shown = list(value.items())[:_MAX_SHOWN]
            text = ", ".join(
                f"{k}: {(0.0 if abs(v) < 1e-12 else v):+.4g}" for k, v in shown
            )
            return text + (", …" if len(value) > _MAX_SHOWN else "")
        return "" if value is None else str(value)

    def __repr__(self):
        columns = self._columns()
        cells = [[self._cell(row.get(c)) for c in columns] for row in self.rows]
        widths = [
            max([len(c)] + [len(r[i]) for r in cells]) for i, c in enumerate(columns)
        ]
        line = lambda items: "  ".join(  # noqa: E731
            str(x).ljust(w) for x, w in zip(items, widths)
        ).rstrip()
        lines = [self.title, line(columns), line("-" * w for w in widths)]
        lines += [line(r) for r in cells]
        if self.footer:
            lines.append(self.footer)
        return "\n".join(lines)

    __str__ = __repr__

    def _repr_html_(self):
        columns = self._columns()
        head = "".join(f"<th>{html.escape(c)}</th>" for c in columns)
        body = "".join(
            "<tr>"
            + "".join(
                f"<td>{html.escape(self._cell(row.get(c)))}</td>" for c in columns
            )
            + "</tr>"
            for row in self.rows
        )
        foot = f"<p>{html.escape(self.footer)}</p>" if self.footer else ""
        return (
            f"<p><b>{html.escape(self.title)}</b></p>"
            f"<table><tr>{head}</tr>{body}</table>{foot}"
        )


def summary(subs):
    """The one-line ``repr`` of a ``CutCircuits`` list."""
    plan = subs.plan
    cuts = []
    if plan.gate_labels:
        cuts.append(f"gate {sorted(set(plan.gate_labels))}")
    if plan.had_cuts:
        cuts.append(f"Hadamard {list(plan.had_cuts)}")
    if plan.wire_cuts:
        cuts.append(f"wire {[(s, t) for s, t, _ in plan.wire_cuts]}")
    if subs.mode == "distribution":
        what = "distribution"
    elif isinstance(subs.observable, list):
        n = len(subs.observable)
        what = f"{n} commuting terms" if n > 1 else "1 term"
        if getattr(subs, "basis", None):
            what += f" (basis {subs.basis})"
    else:
        what = f"observable {subs.observable}"
    return (
        f"<{type(subs).__name__}: {len(subs)} subcircuits, "
        f"{len(plan.fragments)} fragments, {'; '.join(cuts) or 'no cut'}, "
        f"gamma {plan.gamma():.4g}, {what}>"
    )


def describe(subs, values=None):
    """The :class:`Description` of a list (``CutCircuits`` or
    ``CombinedCutSubcircuits``), with each subcircuit's values when given."""
    if values is not None:
        values = list(values)
        if len(values) != len(subs):
            raise ValueError(
                f"expected {len(subs)} results, one per subcircuit, got {len(values)}"
            )
    parts = getattr(subs, "parts", None)
    if parts is None:
        return Description(summary(subs), rows(subs, values), _footer(subs, values))
    table, start = [], 0
    for g, part in enumerate(parts):
        chunk = None if values is None else values[start : start + len(part)]
        table += rows(part, chunk, start, group=g)
        start += len(part)
    titles = "; ".join(f"group {g}: {summary(p)}" for g, p in enumerate(parts))
    return Description(
        f"<{type(subs).__name__}: {len(subs)} subcircuits in {len(parts)} parts>  "
        + titles,
        table,
        _footer(subs, values),
    )


def item_repr(circuit):
    """The ``repr`` of one subcircuit."""
    template = circuit.template
    qubits = " ".join(_qubit_names(template, _template_qubits(template)))
    return (
        f"<Subcircuit: fragment {template.f}, {template.num_qubits} qubits "
        f"[{qubits}], {cut_terms(circuit)}>"
    )


def _footer(subs, values):
    if values is None:
        return ""
    from Qdislib.utils.runtime import compss_wait_on

    try:
        value = compss_wait_on(subs.reconstruct(values))
    except Exception as exc:  # noqa: BLE001 -- a view must not fail
        return f"reconstruct: {type(exc).__name__}: {exc}"
    if isinstance(value, dict):
        top = sorted(value.items(), key=lambda kv: -kv[1])[:3]
        return "reconstructed distribution, most likely: " + ", ".join(
            f"{k}: {p:.4g}" for k, p in top
        )
    return f"reconstructed value: {value:+.10g}"
