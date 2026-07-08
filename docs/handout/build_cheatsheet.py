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

A two-column, two-page reference for a hands-on tutorial: what circuit cutting is,
the cut methods (gate, Hadamard-test, wire, mixed) and their cost, the core API,
backends, parallelism knobs, caching, and the common gotchas. Requires reportlab
(``pip install reportlab``)::

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

    s.append(section("Ways to cut"))
    s.append(gap(3))
    s.append(
        body(
            "<b><font color='#7d6ab0'>Gate cutting</font></b> replaces a two-qubit "
            "gate by a sum of single-qubit <i>measure-and-prepare</i> operations, "
            "splitting the circuit where that gate was. Cuts <b>any two-qubit gate</b> "
            "(CZ/CX natively; others auto-decomposed to CZ). "
            "<b>Cost: 6<super>k</super></b> terms for <i>k</i> cuts."
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "<b><font color='#3f8f9b'>Wire cutting</font></b> cuts a qubit's wire "
            "(its timeline) into “measure here / re-prepare there”. "
            "<b>Cost: 8<super>k</super></b> terms for <i>k</i> cuts."
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "<b><font color='#7d6ab0'>Hadamard gate cutting</font></b> "
            "(<font face='Courier'>gate_cutting_hadamard</font>) cuts a gate with one "
            "<b>local ancilla</b> per fragment instead of mid-circuit measurement — "
            "cuts <b>any</b> gate with the <b>optimal</b> overhead "
            "(&#947;(CZ)=3, &#8594;1 for weak gates). Great for superconducting QPUs."
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
                ("gate_cutting_hadamard(...)", "&lt;O&gt; via the ancilla cut."),
                ("gate_wire_cutting(...)", "Mixed gate + wire cut."),
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

    s.append(section("Hadamard (ancilla) gate cut"))
    s.append(gap(3))
    s.append(
        body(
            "An alternative gate cut (Harrow–Lowe, <i>arXiv:2403.01018</i>): each "
            "fragment carries <b>one extra ancilla</b> running a Hadamard test — <b>no "
            "mid-circuit measurement</b>, and it cuts <b>any</b> two-qubit gate "
            "(iSWAP, fSim, …) at the <b>optimal</b> overhead."
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "val = qd.gate_cutting_hadamard(\n"
            "        circ, ['ISWAP_3'], 'ZZ', software='qibo')\n"
            "# Monte-Carlo (bounded cost for many cuts):\n"
            "val = qd.gate_cutting_hadamard_sampling(\n"
            "        circ, cut, 'ZZ', num_samples=5000)\n"
            "# three-step, runs on any backend / QPU:\n"
            "subs = qd.gate_cutting_hadamard_subcircuits(\n"
            "         circ, cut, 'ZZ')\n"
            "val = subs.reconstruct(\n"
            "  [qd.expectation_value(s) for s in subs])"
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "Fragments are independent (own ancilla, no feed-forward) → they fan out "
            "as PyCOMPSs tasks like every other cut.",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("Mixed wire + gate cutting", ACCENT2))
    s.append(gap(3))
    s.append(
        body(
            "Cut some <b>gates</b> and some <b>wires</b> on the <b>same circuit</b> in "
            "one workflow. The two quasiprobability decompositions tensor, and every "
            "fragment across the combined space runs in parallel."
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "val = qd.gate_wire_cutting(circ,\n"
            "        gates_cut=['CZ_2'],\n"
            "        wire_cut=[('RY_3', 'CZ_4')],\n"
            "        observable='ZZZ', software='qibo')"
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

    s.append(section("Choosing a method"))
    s.append(gap(3))
    s.append(
        api_table(
            [
                ("gate_cutting", "Any 2-qubit gate; distribution too; no ancilla."),
                (
                    "gate_cutting_hadamard",
                    "Any gate at optimal &#947;; no mid-circuit "
                    "measure (best on superconducting QPUs); expectation only.",
                ),
                ("wire_cutting", "Cut a qubit's wire between two subcircuits."),
                ("gate_wire_cutting", "Both cut types on one circuit at once."),
                ("*_sampling", "Reconstruct the full output distribution."),
            ],
            key_frac=0.44,
        )
    )
    s.append(gap(4.5))

    s.append(section("Overhead at a glance", ACCENT2))
    s.append(gap(3))
    s.append(
        api_table(
            [
                ("gate cut / cut", "6 terms  (&#947; = 3 for CZ/CX)."),
                (
                    "wire cut / cut",
                    "8 terms  (or 3<super>n</super>+4<super>n</super> "
                    "with optimized=True).",
                ),
                (
                    "hadamard / cut",
                    "&#947;(CZ)=3, &#947;(SWAP)=7, &#8594;1 for weak "
                    "gates; sampling caps cost at num_samples.",
                ),
                (
                    "general gate",
                    "&#8804; 3 CZ &#8594; up to 6<super>3</super> terms "
                    "(one cut with the Hadamard method).",
                ),
            ],
            key_frac=0.40,
        )
    )
    s.append(gap(4.5))

    s.append(section("Scale the parallelism", ACCENT2))
    s.append(gap(3))
    s.append(
        body("Every subcircuit is an independent task. Extra knobs for big fan-outs:")
    )
    s.append(gap(2))
    s.append(
        api_table(
            [
                ("sync=False", "Return futures; overlap several cuts (wait once)."),
                ("reduce_chunk=n", "Tree reduction instead of one fan-in task."),
                (
                    "optimized=True",
                    "Wire cut: 3<super>n</super>+4<super>n</super> "
                    "fragments, not 8<super>n</super>.",
                ),
                (
                    "num_samples=n",
                    "Hadamard sampling: bound the cost, not 6<super>k</super>.",
                ),
            ],
            key_frac=0.36,
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
                "Any two-qubit gate is cuttable now; a non-CZ/CX gate becomes "
                "<b>3 CZ cuts</b>, so prefer <font face='Courier'>gate_cutting_hadamard</font> "
                "(1 cut, optimal) for exotic gates.",
                "Sampling needs <b>more shots</b> than an expectation value "
                "(quasiprobability overhead).",
                "Identity on a cut qubit: <b>CPU &amp; GPU yes</b>, <b>QPU no</b>; for a "
                "general gate use <font face='Courier'>software='qibo'</font> or the "
                "Hadamard cut.",
                "Use <font face='Courier'>seed=</font> for a reproducible "
                "<font face='Courier'>find_cut</font>.",
            ]
        )
    )
    s.append(gap(4.5))

    s.append(section("Real hardware", ACCENT2))
    s.append(gap(3))
    s.append(
        body("Route the subcircuits to GPUs or QPUs — the cutting code is unchanged:")
    )
    s.append(gap(2))
    s.append(
        code(
            "qd.gate_cutting(circ, cut, gpu=True)      # Aer GPU\n"
            "qd.gate_cutting(circ, cut, qpu=True,\n"
            "        qpu_dict={'MN_Ona': 4})           # BSC QPU\n"
            "qd.gate_cutting(circ, cut, qpu=True,\n"
            "        qpu_dict={'IBM_Quantum': 127})    # IBM"
        )
    )
    s.append(gap(2))
    s.append(
        body(
            "Under PyCOMPSs different fragments can even run on <b>different devices "
            "at once</b> (<font face='Courier'>QdBackEnd(pool)</font>).",
            BODY_MUTED,
        )
    )
    s.append(gap(4.5))

    s.append(section("Bring your own circuit"))
    s.append(gap(3))
    s.append(
        body(
            "Pass a <b>Qibo, Qiskit, CUDA-Q or PennyLane</b> circuit — Qdislib converts "
            "it through a neutral DAG that speaks every backend (U3/CU3/… included). "
            "Or import OpenQASM 2:"
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "circ = qd.from_qasm(open('c.qasm').read())\n"
            "qd.gate_cutting(circ, cut, software='qiskit')"
        )
    )
    s.append(gap(4.5))

    s.append(section("Learn more"))
    s.append(gap(3))
    s.append(
        body(
            "Full docs, the nine tutorial notebooks (Bell states → 25-qubit lattices) "
            "and the API reference:"
        )
    )
    s.append(gap(2))
    s.append(
        api_table(
            [
                ("Docs", "qdislib.readthedocs.io"),
                ("Source", "gitlab.bsc.es/wdc/quantum/qdislib"),
                ("Notebooks", "examples/ — cutting, sampling, scaling, QPUs."),
                ("estimate(...)", "One call: find-cut + cut + reconstruct."),
            ],
            key_frac=0.30,
        )
    )
    s.append(gap(4.5))

    s.append(section("Sanity check", HexColor("#2f8f6b")))
    s.append(gap(3))
    s.append(
        body(
            "Validate a cut against the exact value on a small instance before you "
            "scale it up (cutting reconstructs the same &lt;O&gt;, within shot noise):"
        )
    )
    s.append(gap(2))
    s.append(
        code(
            "exact = qd.analytical_solution(circ, 'ZZ')\n"
            "val   = qd.gate_cutting(circ, cut, 'ZZ',\n"
            "                        shots=8192)\n"
            "assert abs(val - exact) < 0.05\n"
            "qd.draw_cut(circ, cut, 'cut.png')   # see it"
        )
    )
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "Qdislib_cheatsheet.pdf")
    build(out)
    print(f"wrote {out}")
