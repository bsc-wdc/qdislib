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

"""Generate the Qdislib tutorial handout / cheat sheet as a PDF.

A two-column, digestible reference for a hands-on tutorial: what circuit cutting
is, the two cut types and their cost, the core API, backends, caching, and the
common gotchas. Requires reportlab (``pip install reportlab``)::

    python docs/handout/build_cheatsheet.py   # -> Qdislib_cheatsheet.pdf (beside this file)
"""

import os

from reportlab.lib.colors import HexColor, white
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    ListFlowable,
    ListItem,
    PageTemplate,
    Paragraph,
    Preformatted,
    Spacer,
    Table,
    TableStyle,
)

# --------------------------------------------------------------------------- #
# Palette + geometry
# --------------------------------------------------------------------------- #
INK = HexColor("#1a2233")
ACCENT = HexColor("#7d6ab0")  # Qdislib purple (gate cutting) -- soft, not electric
ACCENT2 = HexColor("#3f8f9b")  # muted teal (wire cutting)
CODE_BG = HexColor("#f4f2fa")
RULE = HexColor("#ddd7ee")
MUTED = HexColor("#55607a")

_LOGO = ImageReader(os.path.join(os.path.dirname(__file__), "..", "logo_favicon.png"))

PAGE_W, PAGE_H = A4
MARGIN_X = 11 * mm
MARGIN_TOP = 27 * mm  # room for the banner
MARGIN_BOT = 12 * mm
GUTTER = 7 * mm
COL_W = (PAGE_W - 2 * MARGIN_X - GUTTER) / 2
FRAME_H = PAGE_H - MARGIN_TOP - MARGIN_BOT

# --------------------------------------------------------------------------- #
# Styles
# --------------------------------------------------------------------------- #
BODY = ParagraphStyle(
    "body", fontName="Helvetica", fontSize=8, leading=10.4, textColor=INK
)
BODY_MUTED = ParagraphStyle("bodym", parent=BODY, textColor=MUTED)
HDR = ParagraphStyle(
    "hdr", fontName="Helvetica-Bold", fontSize=9.5, leading=11, textColor=white
)
LEAD = ParagraphStyle(
    "lead", parent=BODY, fontSize=8.4, leading=11, spaceAfter=2, textColor=INK
)
CODE = ParagraphStyle(
    "code", fontName="Courier", fontSize=6.7, leading=8.2, textColor=INK
)
KEY = ParagraphStyle(
    "key", fontName="Courier-Bold", fontSize=7, leading=9, textColor=ACCENT
)
VAL = ParagraphStyle("val", parent=BODY, fontSize=7.5, leading=9)


def section(title, color=ACCENT):
    """A full-width coloured section-header bar."""
    bar = Table([[Paragraph(title.upper(), HDR)]], colWidths=[COL_W])
    bar.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), color),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    return bar


def body(text, style=BODY):
    return Paragraph(text, style)


def code(src):
    """A shaded, monospaced code box (lines must fit the column width)."""
    box = Table([[Preformatted(src, CODE)]], colWidths=[COL_W])
    box.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), CODE_BG),
                ("BOX", (0, 0), (-1, -1), 0.5, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return box


def api_table(rows, key_frac=0.42):
    """Two-column reference table: monospaced key -> description.

    ``key_frac`` is the fraction of the column width given to the (monospaced)
    key; widen it for long function signatures (e.g. the Core API table).
    """
    data = [[Paragraph(k, KEY), Paragraph(v, VAL)] for k, v in rows]
    t = Table(data, colWidths=[COL_W * key_frac, COL_W * (1 - key_frac)])
    t.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LINEBELOW", (0, 0), (-1, -2), 0.3, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
            ]
        )
    )
    return t


def bullets(items):
    return ListFlowable(
        [ListItem(Paragraph(t, BODY), leftIndent=8, value="•") for t in items],
        bulletType="bullet",
        start="•",
        leftIndent=8,
        bulletFontSize=7,
    )


def gap(h=4):
    return Spacer(1, h)


# --------------------------------------------------------------------------- #
# Page furniture (banner + footer), drawn under the flowing columns
# --------------------------------------------------------------------------- #
def _decorate(canvas, doc):
    canvas.saveState()
    # Top banner
    canvas.setFillColor(ACCENT)
    canvas.rect(0, PAGE_H - 20 * mm, PAGE_W, 20 * mm, stroke=0, fill=1)
    canvas.setFillColor(ACCENT2)
    canvas.rect(0, PAGE_H - 20 * mm, PAGE_W, 1.6 * mm, stroke=0, fill=1)
    canvas.setFillColor(white)
    canvas.setFont("Helvetica-Bold", 15)
    canvas.drawString(
        MARGIN_X, PAGE_H - 12 * mm, "Qdislib — Circuit Cutting Cheat Sheet"
    )
    canvas.setFont("Helvetica", 8.2)
    canvas.drawString(
        MARGIN_X,
        PAGE_H - 17 * mm,
        "Run circuits bigger than your device: cut, execute the pieces, recombine.",
    )
    # Logo emblem on a white chip (the emblem is purple-on-transparent, so it
    # needs a light backing to read against the banner), top-right of the banner.
    lw = 15 * mm
    lx = PAGE_W - MARGIN_X - lw
    ly = PAGE_H - 17.5 * mm
    canvas.setFillColor(white)
    canvas.roundRect(
        lx - 1.6 * mm,
        ly - 1.6 * mm,
        lw + 3.2 * mm,
        lw + 3.2 * mm,
        2 * mm,
        stroke=0,
        fill=1,
    )
    canvas.drawImage(_LOGO, lx, ly, lw, lw, mask="auto", preserveAspectRatio=True)
    # Column divider
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.4)
    x = MARGIN_X + COL_W + GUTTER / 2
    canvas.line(x, MARGIN_BOT, x, PAGE_H - MARGIN_TOP + 3 * mm)
    # Footer
    canvas.setFillColor(MUTED)
    canvas.setFont("Helvetica", 6.8)
    canvas.drawString(
        MARGIN_X,
        7 * mm,
        "gitlab.bsc.es/wdc/quantum/qdislib  ·  Barcelona Supercomputing Center",
    )
    canvas.drawRightString(PAGE_W - MARGIN_X, 7 * mm, f"p. {doc.page}")
    canvas.restoreState()


def build(path):
    doc = BaseDocTemplate(
        path,
        pagesize=A4,
        leftMargin=MARGIN_X,
        rightMargin=MARGIN_X,
        topMargin=MARGIN_TOP,
        bottomMargin=MARGIN_BOT,
        title="Qdislib Circuit Cutting Cheat Sheet",
        author="Barcelona Supercomputing Center",
    )
    left = Frame(
        MARGIN_X,
        MARGIN_BOT,
        COL_W,
        FRAME_H,
        id="left",
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
    )
    right = Frame(
        MARGIN_X + COL_W + GUTTER,
        MARGIN_BOT,
        COL_W,
        FRAME_H,
        id="right",
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
    )
    doc.addPageTemplates(
        [PageTemplate(id="two-col", frames=[left, right], onPage=_decorate)]
    )
    doc.build(_story())


# --------------------------------------------------------------------------- #
# Content
# --------------------------------------------------------------------------- #
def _story():
    s = []

    s.append(
        body(
            "<b>Qdislib</b> is a distributed circuit-cutting library: it splits a "
            "quantum circuit into smaller subcircuits that fit your hardware, runs "
            "them independently (in parallel, even across backends or QPUs) and "
            "recombines the results classically. It is <b>backend-agnostic</b> — the "
            "same code drives Qibo, Qiskit and CUDA-Q — and scales with PyCOMPSs.",
            LEAD,
        )
    )
    s.append(gap(4.5))

    s.append(section("Two ways to cut"))
    s.append(gap(3))
    s.append(
        body(
            "<b><font color='#7d6ab0'>Gate cutting</font></b> replaces a two-qubit "
            "gate (CZ / CX) by a sum of single-qubit <i>measure-and-prepare</i> "
            "operations, splitting the circuit where that gate was. "
            "<b>Cost: 6<super>k</super></b> subcircuit terms for <i>k</i> cuts."
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "<b><font color='#3f8f9b'>Wire cutting</font></b> cuts a qubit's wire "
            "(its timeline) into “measure here / re-prepare there”. "
            "<b>Cost: 8<super>k</super></b> subcircuit terms for <i>k</i> cuts."
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "Every cut multiplies the work <b>exponentially</b>, so use as few, and "
            "as well-placed, as possible. Preview the blow-up with "
            "<font face='Courier'>cut_cost(cut)</font> before you run.",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("The workflow"))
    s.append(gap(3))
    s.append(
        body(
            "Get one <b>expectation value</b> &lt;O&gt; in a single call, or use the "
            "<b>three-step</b> form to control how each piece runs (cluster, QPU):"
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "cut = qd.find_cut(circ, gate_cut=True)\n"
            "val = qd.gate_cutting(circ, cut,\n"
            "          observables='ZZZZ', shots=8192)\n"
            "\n"
            "subs = qd.gate_cutting_subcircuits(\n"
            "         circ, cut, observables='ZZZZ')\n"
            "vals = [qd.expectation_value(s)\n"
            "        for s in subs]     # run anywhere\n"
            "val = subs.reconstruct(vals)"
        )
    )
    s.append(gap(3))
    s.append(
        body(
            "<font face='Courier'>estimate(circ, 'ZZZZ', max_qubits=n)</font> does "
            "find-cut + cut + reconstruct in one call.",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("Quickstart"))
    s.append(gap(3))
    s.append(
        code(
            "import Qdislib as qd\n"
            "from qibo import Circuit, gates\n"
            "\n"
            "circ = Circuit(4)\n"
            "circ.add(gates.H(0))\n"
            "circ.add(gates.CZ(0, 1))\n"
            "circ.add(gates.CZ(1, 2))\n"
            "circ.add(gates.CZ(2, 3))\n"
            "\n"
            "cut = qd.find_cut(circ, gate_cut=True,\n"
            "        wire_cut=False, seed=0)\n"
            "val = qd.gate_cutting(circ, cut,\n"
            "        observables='ZZZZ', shots=8192)"
        )
    )
    s.append(gap(4.5))

    s.append(section("Core API"))
    s.append(gap(3))
    s.append(
        api_table(
            [
                ("find_cut(circ, ...)", "Choose where to cut (returns it)."),
                ("gate_cutting(circ, cut)", "&lt;O&gt; via gate cutting."),
                ("wire_cutting(circ, cut)", "&lt;O&gt; via wire cutting."),
                ("..._sampling(circ, cut)", "Full output distribution."),
                ("..._subcircuits(...)", "3-step: pieces + .reconstruct()."),
                ("expectation_value(sub)", "Run one subcircuit."),
                ("estimate(circ, ...)", "One-call cut + run + reconstruct."),
                ("cut_cost(cut)", "Preview terms = base**#cuts."),
                ("analytical_solution(...)", "Exact &lt;O&gt; (to check)."),
                ("draw_cut(circ, cut, path)", "Visualise the cut."),
                ("from_qasm(qasm, ...)", "Import any OpenQASM 2."),
            ],
            key_frac=0.55,
        )
    )
    s.append(gap(4.5))

    s.append(section("Sampling the distribution", ACCENT2))
    s.append(gap(3))
    s.append(
        body(
            "Reconstruct the full <b>output distribution</b> "
            "<font face='Courier'>{bitstring: prob}</font> — to read out an image "
            "or inspect a state — instead of a single number:"
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "dist = qd.gate_cutting_sampling(\n"
            "         circ, cut, shots=200_000)\n"
            "# three-step (parallel under PyCOMPSs):\n"
            "subs = qd.gate_cutting_sampling_subcircuits(\n"
            "         circ, cut)\n"
            "dist = subs.reconstruct(\n"
            "  [qd.sample_subcircuit(s) for s in subs])"
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "<font face='Courier'>wire_cutting_sampling</font> too. Sampling carries "
            "a quasiprobability overhead, so it needs <b>more shots</b> than an "
            "expectation value.",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("find_cut options"))
    s.append(gap(3))
    s.append(
        api_table(
            [
                ("max_qubits=n", "Max qubits per subcircuit (device size)."),
                ("max_cuts=k", "Cap the number of cuts."),
                ("gate_cut / wire_cut", "Which cut types to consider."),
                ("implementation=", "'qdislib' (default) or 'ibm-ckt'."),
                ("seed=n", "Reproducible cut selection."),
            ]
        )
    )
    s.append(gap(4.5))

    s.append(section("Observables"))
    s.append(gap(3))
    s.append(
        body(
            "Per-qubit Pauli string over <b>I X Y Z</b> (qubit-0-first): "
            "<font face='Courier'>'ZZZZ'</font>, "
            "<font face='Courier'>'ZIZIZ'</font>, "
            "<font face='Courier'>'XZZY'</font>. A Hamiltonian is a weighted sum:"
        )
    )
    s.append(gap(2))
    s.append(
        code("[(0.5, 'ZZII'), (-1.0, 'IIZZ')]\n" "# or  {'ZZII': 0.5, 'IIZZ': -1.0}")
    )
    s.append(gap(2))
    s.append(
        body(
            "<font face='Courier'>I</font> ignores a qubit — and works even on a "
            "<i>cut</i> qubit (CPU &amp; GPU).",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("Backends", ACCENT2))
    s.append(gap(3))
    s.append(
        api_table(
            [
                ("software=", "'qibo' | 'qiskit' | 'cudaq'."),
                ("gpu=True", "Aer GPU (needs CUDA)."),
                ("qpu=True, qpu_dict=", "{'MN_Ona': n} / {'IBM_Quantum': 127}."),
                ("available_backends()", "List registered backends."),
                ("build_execution_backend", "How software/gpu/qpu pick one."),
                ("QdBackEnd(pool)", "Load-balance across devices."),
                ("register_backend(...)", "Plug in your own (BaseBackend)."),
            ]
        )
    )
    s.append(gap(4.5))

    s.append(section("Caching"))
    s.append(gap(3))
    s.append(
        code(
            "cache = qd.default_circuit_cache(\n"
            "          'memory', max_entries=5000)\n"
            "qd.gate_cutting(circ, cut, cache=cache)"
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "The <b>semantic cache</b> reuses functionally-equivalent subcircuits "
            "(e.g. H·H·H = H). Stores: "
            "<font face='Courier'>'memory'</font> / "
            "<font face='Courier'>'lmdb'</font> / "
            "<font face='Courier'>'redis'</font>. "
            "<font face='Courier'>set_log_level('DEBUG')</font> for diagnostics.",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("Distributed with PyCOMPSs", ACCENT2))
    s.append(gap(3))
    s.append(
        body(
            "Qdislib is built for HPC: every subcircuit is a PyCOMPSs "
            "<font face='Courier'>@task</font>, so the pieces run <b>in parallel "
            "across cores and nodes</b>. The three-step workflow fans them out "
            "automatically — the same code scales from a laptop to a supercomputer."
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "import pycompss.interactive as ipycompss\n"
            "ipycompss.start()   # or:  runcompss app.py\n"
            "subs = qd.gate_cutting_subcircuits(circ, cut)\n"
            "vals = [qd.expectation_value(s) for s in subs]\n"
            "val = subs.reconstruct(vals)   # auto-parallel\n"
            "ipycompss.stop()"
        )
    )
    s.append(gap(4.5))

    s.append(section("Gotchas", HexColor("#b5342b")))
    s.append(gap(3))
    s.append(
        bullets(
            [
                "Cost is exponential in #cuts (6<super>k</super> / 8<super>k</super>) — "
                "fewer, well-placed cuts win; preview with <font face='Courier'>cut_cost</font>.",
                "Gate cutting cuts only <b>CZ / CX / CNOT</b> — decompose other "
                "two-qubit gates first.",
                "Sampling needs <b>more shots</b> than an expectation value "
                "(quasiprobability overhead).",
                "Identity on a cut qubit: <b>CPU &amp; GPU yes</b>, <b>QPU no</b> "
                "(hardware returns only aggregated counts).",
                "Use <font face='Courier'>seed=</font> for a reproducible "
                "<font face='Courier'>find_cut</font>.",
            ]
        )
    )
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "Qdislib_cheatsheet.pdf")
    build(out)
    print(f"wrote {out}")
