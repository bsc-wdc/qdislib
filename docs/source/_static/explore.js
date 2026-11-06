/* The Explore page. The circuit playground computes everything here, with explore_engine.js;
   the other sections replay window.QD_EXPLORE, stored runs of real Qdislib made by
   docs/build_explore_data.py. No other network request. */
(function () {
  "use strict";
  var E = window.QD_ENGINE, D = window.QD_EXPLORE;

  function el(tag, attrs, text) {
    var node = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (key) { node.setAttribute(key, attrs[key]); });
    if (text !== undefined) node.textContent = text;
    return node;
  }
  function esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); }
  function num(x, digits) { return (x >= 0 ? "+" : "−") + Math.abs(x).toFixed(digits); }
  function count(n) {
    if (n < 1e6) return Math.round(n).toLocaleString("en-US");
    var e = Math.floor(Math.log10(n));
    return (n / Math.pow(10, e)).toFixed(1) + "e" + e;
  }
  function small(x) { return x < 1e-12 ? "below 1e-12, rounding only" : x.toExponential(0); }
  function gammaText(g) { return g >= 1000 ? count(g) : String(Math.round(g * 1000) / 1000); }
  function stat(box, value, label) {
    var item = box.appendChild(el("div", { class: "qd-stat" }));
    item.appendChild(el("b", {}, value));
    item.appendChild(el("span", {}, label));
  }
  function toggles(box, items, active, onPick) {
    box.textContent = "";
    items.forEach(function (item) {
      var b = box.appendChild(el("button", { type: "button", "aria-pressed": String(item.value === active) }, item.label));
      b.addEventListener("click", function () { onPick(item.value); });
    });
  }
  function fragColor(f) { return "var(--qdx-f" + (f % 8) + ")"; }
  function shortName(op) { return op.name.toUpperCase(); }

  /* ---- A circuit as SVG, coloured by fragment ---- */
  var CELL = 50, LEFT = 40, TOP = 18, ROW = 46;

  function columns(n, ops) {
    var free = [], out = [];
    for (var q = 0; q < n; q++) free[q] = 0;
    ops.forEach(function (op) {
      var low = Math.min.apply(null, op.qubits), high = Math.max.apply(null, op.qubits), c = 0;
      for (var q = low; q <= high; q++) c = Math.max(c, free[q]);
      for (q = low; q <= high; q++) free[q] = c + 1;
      out.push(c);
    });
    return out;
  }

  /* opts: n, ops, cols (per op), cuts, plan, width (columns), labels (per qubit), hits (bool) */
  function drawCircuit(o) {
    var width = o.width, W = LEFT + width * CELL + 12, H = TOP + o.n * ROW + 6;
    var y = function (q) { return TOP + q * ROW + ROW / 2; };
    var x = function (c) { return LEFT + c * CELL + CELL / 2; };
    var s = '<svg viewBox="0 0 ' + W + " " + H + '" width="' + W + '" height="' + H + '" style="max-width:' + W + 'px" role="img" aria-label="' + esc(o.aria || "A circuit") + '">';
    var p = o.plan;
    var hasPlan = p && !p.invalid.length;

    /* the lines, one colour per segment */
    var wireX = {};
    (o.cuts || []).forEach(function (cut) {
      if (cut.kind === "wire") wireX[cut.qubit + ":" + cut.from] = (o.cols[cut.from] + o.cols[cut.to] + 1) / 2 * CELL + LEFT;
    });
    for (var q = 0; q < o.n; q++) {
      var label = o.labels ? o.labels[q] : "q" + q;
      s += '<text class="q" x="6" y="' + (y(q) + 4) + '">' + esc(label) + "</text>";
      var segs = hasPlan ? p.segments.filter(function (g) { return g.qubit === q; }) : [];
      var start = LEFT - 8;
      if (!segs.length) {
        s += '<line x1="' + start + '" x2="' + (W - 6) + '" y1="' + y(q) + '" y2="' + y(q) + '" stroke="var(--qd-line)" stroke-width="2"/>';
      }
      segs.forEach(function (g) {
        var end = W - 6;
        if (g.outgoing) end = wireX[q + ":" + g.outgoing.from];
        var colour = g.fragment === undefined ? "var(--qd-line)" : fragColor(g.fragment);
        s += '<line x1="' + start + '" x2="' + end + '" y1="' + y(q) + '" y2="' + y(q) + '" stroke="' + colour + '" stroke-width="2.5"/>';
        start = end;
      });
    }
    /* the gates */
    var cutOf = {};
    (o.cuts || []).forEach(function (cut) { if (cut.kind !== "wire") cutOf[cut.op] = cut; });
    o.ops.forEach(function (op, i) {
      var c = x(o.cols[i]), cut = cutOf[i];
      var fragOf = function (side) {
        if (!hasPlan) return null;
        var seg = p.segments[p.opSegments[i][side]];
        return seg && seg.fragment !== undefined ? seg.fragment : null;
      };
      if (op.qubits.length === 2) {
        var a = y(op.qubits[0]), b = y(op.qubits[1]);
        var stroke = cut ? (cut.kind === "hadamard" ? "var(--qdx-had)" : "var(--qdx-cut)") : "var(--color-foreground-secondary, #555)";
        s += '<line x1="' + c + '" x2="' + c + '" y1="' + a + '" y2="' + b + '" stroke="' + stroke + '" stroke-width="2"' + (cut ? ' stroke-dasharray="4 3"' : "") + "/>";
        if (cut) s += '<text class="scissors" x="' + (c + 3) + '" y="' + ((a + b) / 2 + 5) + '" style="fill:' + stroke + '">✂</text>';
      }
      op.qubits.forEach(function (qubit, side) {
        var f = fragOf(side), fill = f === null ? "#667" : fragColor(f);
        var cy = y(qubit);
        s += '<rect x="' + (c - 17) + '" y="' + (cy - 14) + '" width="34" height="28" rx="6" fill="' + fill + '"' +
          (cut ? ' stroke="' + (cut.kind === "hadamard" ? "var(--qdx-had)" : "var(--qdx-cut)") + '" stroke-width="3" stroke-dasharray="4 2"' : "") + "/>";
        s += '<text class="g" x="' + c + '" y="' + (cy + 4) + '" text-anchor="middle">' + esc(shortName(op)) + "</text>";
      });
      s += '<text class="i" x="' + c + '" y="' + (y(Math.min.apply(null, op.qubits)) - 17) + '" text-anchor="middle">' + (i + 1) + "</text>";
    });
    /* wire cuts */
    (o.cuts || []).forEach(function (cut) {
      if (cut.kind !== "wire") return;
      var cx = wireX[cut.qubit + ":" + cut.from];
      s += '<text class="scissors" x="' + cx + '" y="' + (y(cut.qubit) + 5) + '" text-anchor="middle">✂</text>';
    });
    if (o.overlay) s += o.overlay({ x: x, y: y });
    return s + "</svg>";
  }

  function board(host, svg, fit) {
    host.className = "qdx-board" + (fit ? " qdx-fit" : "");
    host.innerHTML = svg;
  }

  /* ---- 1. The playground ---- */
  function playground(root) {
    var COLS = 12, MAX_CUTS = 4;
    var PALETTE = [
      ["h", "H"], ["x", "X"], ["s", "S"], ["t", "T"], ["rx", "RX"], ["ry", "RY"], ["rz", "RZ"],
      ["cz", "CZ"], ["cx", "CX"], ["rzz", "RZZ"]
    ];
    var ANGLED = { rx: 1, ry: 1, rz: 1, rzz: 1 }, TWO = { cz: 1, cx: 1, rzz: 1 };
    var state = { n: 5, ops: [], cuts: [], observable: "ZZZZZ", tool: "cut", angle: 0.8, pending: null, shots: 10, seed: 1 };
    var uid = 1;

    function add(name, qubits, col, theta) {
      var op = { id: uid++, name: name, qubits: qubits, col: col };
      if (ANGLED[name]) op.theta = theta === undefined ? state.angle : theta;
      state.ops.push(op);
      return op;
    }
    function presets() {
      return {
        chain: function () {
          state.n = 5; state.observable = "ZZZZZ"; state.ops = []; state.cuts = [];
          [0.8, 0.5, 0.6, 0.4, 0.3].forEach(function (a, q) { add("ry", [q], 0, a); });
          var mid = null;
          for (var q = 0; q < 4; q++) { var op = add("cz", [q, q + 1], q + 1); if (q === 1) mid = op; }
          [0.2, 0.1, 0.7, 0.3, 0.9].forEach(function (a, q) { add("rx", [q], 5 + q % 2, a); });
          state.cuts = [{ kind: "gate", opId: mid.id }];
        },
        blocks: function () {
          state.n = 6; state.observable = "ZZZZZZ"; state.ops = []; state.cuts = [];
          [0, 1, 2, 3, 4, 5].forEach(function (q) { add("h", [q], 0); });
          add("rzz", [0, 1], 1, 1.2); add("rzz", [3, 4], 1, 1.2);
          add("rzz", [1, 2], 2, 1.2); add("rzz", [4, 5], 2, 1.2);
          var a = add("rzz", [2, 3], 3, 0.2), b = add("rzz", [0, 5], 4, 0.2);
          [0, 1, 2, 3, 4, 5].forEach(function (q) { add("rx", [q], 5, 0.5); });
          state.cuts = [{ kind: "gate", opId: a.id }, { kind: "gate", opId: b.id }];
        },
        wire: function () {
          state.n = 4; state.observable = "ZZZZ"; state.ops = []; state.cuts = [];
          [0, 1, 2, 3].forEach(function (q) { add("ry", [q], 0, 0.5 + 0.2 * q); });
          add("cz", [0, 1], 1);
          var a = add("cz", [1, 2], 2), b = add("cz", [2, 3], 3);
          state.cuts = [{ kind: "wire", qubit: 2, fromId: a.id, toId: b.id }];
        },
        ghz: function () {
          state.n = 6; state.observable = "ZZZZZZ"; state.ops = []; state.cuts = [];
          add("h", [0], 0);
          for (var q = 0; q < 5; q++) add("cx", [q, q + 1], q + 1);
          state.cuts = [{ kind: "hadamard", opId: state.ops[3].id }];
        },
        random: function () {
          var rng = E.mulberry((state.seed = state.seed + 1)), n = 4 + Math.floor(rng() * 3);
          state.n = n; state.ops = []; state.cuts = [];
          state.observable = new Array(n + 1).join("Z");
          for (var c = 0; c < 8; c++) {
            var used = {};
            for (var t = 0; t < n; t++) {
              var q = Math.floor(rng() * n);
              if (used[q]) continue;
              if (rng() < 0.45) {
                var r = (q + 1 + Math.floor(rng() * (n - 1))) % n, low = Math.min(q, r), high = Math.max(q, r), clear = true;
                for (var k = low; k <= high; k++) if (used[k]) clear = false;
                if (!clear) continue;
                for (k = low; k <= high; k++) used[k] = 1;
                var name = ["cz", "cx", "rzz"][Math.floor(rng() * 3)];
                add(name, [q, r], c, 0.3 + 2.4 * rng());
              } else {
                used[q] = 1;
                var one = PALETTE[Math.floor(rng() * 7)][0];
                add(one, [q], c, 0.3 + 2.4 * rng());
              }
            }
          }
        }
      };
    }

    /* The circuit as Qdislib numbers it: by column, then from the top. */
    function ordered() {
      return state.ops.slice().sort(function (a, b) {
        return a.col - b.col || Math.min.apply(null, a.qubits) - Math.min.apply(null, b.qubits);
      });
    }
    function indexOf(ops, id) {
      for (var i = 0; i < ops.length; i++) if (ops[i].id === id) return i;
      return -1;
    }
    /* Consecutive gates on a qubit that share only it: what a wire cut can name. */
    function wireEdges(ops) {
      var last = {}, out = [];
      ops.forEach(function (op, i) {
        op.qubits.forEach(function (q) {
          if (q in last) {
            var other = ops[last[q]];
            var shared = other.qubits.filter(function (x) { return op.qubits.indexOf(x) >= 0; }).length;
            if (shared === 1) out.push({ qubit: q, from: last[q], to: i });
          }
          last[q] = i;
        });
      });
      return out;
    }
    function dropStaleCuts(ops) {
      var edges = wireEdges(ops);
      state.cuts = state.cuts.filter(function (cut) {
        if (cut.opId !== undefined) return indexOf(ops, cut.opId) >= 0;
        var a = indexOf(ops, cut.fromId), b = indexOf(ops, cut.toId);
        return edges.some(function (e) { return e.qubit === cut.qubit && e.from === a && e.to === b; });
      });
    }
    function occupied(q, col) {
      return state.ops.filter(function (op) {
        return op.col === col && q >= Math.min.apply(null, op.qubits) && q <= Math.max.apply(null, op.qubits);
      })[0];
    }

    root.textContent = "";
    root.appendChild(el("span", { class: "qd-label" }, "Start from"));
    var presetBox = root.appendChild(el("div", { class: "qdx-tools" }));
    root.appendChild(el("span", { class: "qd-label" }, "Then tap a tool, and tap the circuit"));
    var toolBox = root.appendChild(el("div", { class: "qdx-tools" }));
    var help = root.appendChild(el("p", { class: "qd-note" }));
    var angleRow = root.appendChild(el("div", { class: "qdx-row" }));
    var angleLabel = angleRow.appendChild(el("label", { for: "qdx-angle" }));
    var angle = angleRow.appendChild(el("input", { id: "qdx-angle", type: "range", min: "0.1", max: "3", step: "0.1", value: "0.8", "aria-label": "Angle of the next rotation" }));
    angle.addEventListener("input", function () { state.angle = Number(angle.value); angleLabel.textContent = "Angle of the next rotation: " + state.angle.toFixed(1); });
    angleLabel.textContent = "Angle of the next rotation: 0.8";
    var boardHost = root.appendChild(el("div"));
    var chips = root.appendChild(el("div", { class: "qdx-chips" }));
    var warn = root.appendChild(el("div", { class: "qdx-warn", hidden: "" }));
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var check = root.appendChild(el("div", { class: "qd-check" }));
    var note = root.appendChild(el("p", { class: "qd-note" }));

    root.appendChild(el("span", { class: "qd-label" }, "With shots"));
    var shotRow = root.appendChild(el("div", { class: "qd-shots" }));
    var shotSlider = shotRow.appendChild(el("input", { type: "range", min: "6", max: "16", step: "1", value: "10", "aria-label": "Shots per subcircuit" }));
    var again = shotRow.appendChild(el("button", { type: "button" }, "Draw again"));
    var shotStats = root.appendChild(el("div", { class: "qd-stats" }));
    var strip = root.appendChild(el("div"));
    var shotNote = root.appendChild(el("p", { class: "qd-note" },
      "Shots are simulated here: each distinct fragment circuit gives a +1 or −1 per shot, drawn with its exact probability. Dots are 30 independent runs."));

    var termsBox = root.appendChild(el("details"));
    termsBox.appendChild(el("summary", {}, "The terms of each cut"));
    var termsBody = termsBox.appendChild(el("div"));
    var codeBox = root.appendChild(el("details"));
    codeBox.appendChild(el("summary", {}, "The same thing with Qdislib"));
    var codeBody = codeBox.appendChild(el("pre"));

    toggles(presetBox, [
      { value: "chain", label: "Chain" }, { value: "blocks", label: "Weak gates" }, { value: "wire", label: "Wire cut" },
      { value: "ghz", label: "GHZ" }, { value: "random", label: "Random" }, { value: "clear", label: "Clear" }
    ], null, function (v) {
      if (v === "clear") { state.ops = []; state.cuts = []; state.pending = null; }
      else presets()[v]();
      state.pending = null;
      render();
    });
    var toolItems = PALETTE.map(function (p) { return { value: p[0], label: p[1] }; }).concat([
      { value: "cut", label: "✂ Cut" }, { value: "erase", label: "Erase" }
    ]);
    function drawTools() {
      toggles(toolBox, toolItems, state.tool, function (v) { state.tool = v; state.pending = null; render(); });
      var more = toolBox.appendChild(el("button", { type: "button" }, "+ qubit"));
      more.addEventListener("click", function () {
        if (state.n < 6) { state.n += 1; state.observable += "Z"; render(); }
      });
      var less = toolBox.appendChild(el("button", { type: "button" }, "− qubit"));
      less.addEventListener("click", function () {
        if (state.n > 2) {
          state.n -= 1;
          state.observable = state.observable.slice(0, state.n);
          state.ops = state.ops.filter(function (op) { return op.qubits.every(function (q) { return q < state.n; }); });
          render();
        }
      });
    }

    function cycleCut(op) {
      var i = -1;
      state.cuts.forEach(function (c, k) { if (c.opId === op.id) i = k; });
      if (i < 0) {
        if (state.cuts.length >= MAX_CUTS) return;
        state.cuts.push({ kind: "gate", opId: op.id });
      } else if (state.cuts[i].kind === "gate") state.cuts[i].kind = "hadamard";
      else state.cuts.splice(i, 1);
    }
    function toggleWire(edge, ops) {
      var fromId = ops[edge.from].id, toId = ops[edge.to].id, i = -1;
      state.cuts.forEach(function (c, k) { if (c.kind === "wire" && c.fromId === fromId && c.toId === toId && c.qubit === edge.qubit) i = k; });
      if (i >= 0) state.cuts.splice(i, 1);
      else if (state.cuts.length < MAX_CUTS) state.cuts.push({ kind: "wire", qubit: edge.qubit, fromId: fromId, toId: toId });
    }

    function tapCell(q, col) {
      var tool = state.tool, here = occupied(q, col);
      if (tool === "erase") {
        if (here) state.ops.splice(state.ops.indexOf(here), 1);
      } else if (tool === "cut") {
        if (here && here.qubits.length === 2) cycleCut(here);
      } else if (TWO[tool]) {
        if (!state.pending) {
          if (!here) state.pending = { q: q, col: col };
        } else if (state.pending.col !== col || state.pending.q === q) {
          state.pending = here ? null : { q: q, col: col };
        } else {
          var low = Math.min(q, state.pending.q), high = Math.max(q, state.pending.q), free = true;
          for (var k = low; k <= high; k++) if (occupied(k, col)) free = false;
          if (free) add(tool, [state.pending.q, q], col);
          state.pending = null;
        }
      } else if (!here) {
        add(tool, [q], col);
      }
      render();
    }

    function render() {
      drawTools();
      var ops = ordered();
      dropStaleCuts(ops);
      var plain = ops.map(function (op) {
        var o = { name: op.name, qubits: op.qubits };
        if (op.theta !== undefined) o.theta = op.theta;
        return o;
      });
      var cols = ops.map(function (op) { return op.col; });
      var cuts = state.cuts.map(function (c) {
        return c.opId !== undefined
          ? { kind: c.kind, op: indexOf(ops, c.opId) }
          : { kind: "wire", qubit: c.qubit, from: indexOf(ops, c.fromId), to: indexOf(ops, c.toId) };
      });
      var circuit = { n: state.n, ops: plain };
      var plan = E.plan(circuit, cuts, state.observable);
      var edges = wireEdges(ops);

      /* The drawing, with its tap targets */
      var labels = state.observable.split("").map(function (p, q) { return p + " q" + q; });
      var svg = drawCircuit({
        n: state.n, ops: plain, cols: cols, cuts: cuts, plan: plan, width: COLS,
        labels: state.observable.split(""),
        aria: "The circuit. Each colour is a fragment.",
        overlay: function (g) {
          var s = "";
          if (state.pending) s += '<circle class="pending" cx="' + g.x(state.pending.col) + '" cy="' + g.y(state.pending.q) + '" r="16"/>';
          for (var q = 0; q < state.n; q++) {
            s += '<rect class="qdx-hit" data-obs="' + q + '" x="0" y="' + (g.y(q) - ROW / 2) + '" width="' + LEFT + '" height="' + ROW + '"/>';
            for (var c = 0; c < COLS; c++) {
              s += '<rect class="qdx-hit" data-q="' + q + '" data-c="' + c + '" x="' + (g.x(c) - CELL / 2) + '" y="' + (g.y(q) - ROW / 2) + '" width="' + CELL + '" height="' + ROW + '"/>';
            }
          }
          if (state.tool === "cut") {
            edges.forEach(function (edge, i) {
              var cx = LEFT + ((cols[edge.from] + cols[edge.to] + 1) / 2) * CELL;
              var on = cuts.some(function (k) { return k.kind === "wire" && k.qubit === edge.qubit && k.from === edge.from && k.to === edge.to; });
              if (!on) s += '<circle cx="' + cx + '" cy="' + g.y(edge.qubit) + '" r="5" fill="var(--color-background-primary, #fff)" stroke="var(--qdx-cut)" stroke-width="2"/>';
              s += '<rect class="qdx-hit" data-edge="' + i + '" x="' + (cx - 13) + '" y="' + (g.y(edge.qubit) - 16) + '" width="26" height="32"/>';
            });
          }
          return s;
        }
      });
      board(boardHost, svg);

      /* Chips: the same cuts with a large target */
      chips.textContent = "";
      if (state.tool === "cut") {
        ops.forEach(function (op, i) {
          if (op.qubits.length !== 2) return;
          var cut = cuts.filter(function (k) { return k.op === i; })[0];
          var text = shortName(op) + "_" + (i + 1) + (cut ? (cut.kind === "gate" ? ": gate cut" : ": Hadamard cut") : "");
          var b = chips.appendChild(el("button", { type: "button", "aria-pressed": String(!!cut) }, text));
          b.addEventListener("click", function () { cycleCut(op); render(); });
        });
        edges.forEach(function (edge) {
          var on = cuts.some(function (k) { return k.kind === "wire" && k.qubit === edge.qubit && k.from === edge.from && k.to === edge.to; });
          var text = "wire q" + edge.qubit + ": " + (edge.from + 1) + "→" + (edge.to + 1);
          var b = chips.appendChild(el("button", { type: "button", "aria-pressed": String(on) }, text));
          b.addEventListener("click", function () { toggleWire(edge, ops); render(); });
        });
      }
      help.textContent = state.tool === "cut"
        ? "Tap a two-qubit gate once for a gate cut, again for a Hadamard cut, again to remove it. Tap a small circle on a line for a wire cut. At most " + MAX_CUTS + " cuts. Tap a letter on the left to change what is measured on that qubit."
        : state.tool === "erase" ? "Tap a gate to remove it."
        : TWO[state.tool] ? "Tap one qubit, then another in the same column. The first is the control of a CX."
        : "Tap an empty place on the circuit.";
      angleRow.style.display = ANGLED[state.tool] ? "" : "none";

      /* What it costs, and the result */
      var observable = state.observable;
      stats.textContent = "";
      check.textContent = "";
      note.textContent = "";
      var exact = E.exact(circuit, observable);
      var bad = plan.invalid.length > 0;
      warn.hidden = !bad;
      if (bad) {
        warn.textContent = "This cut does not separate its two sides: another path still joins them. Qdislib refuses it with InvalidCutError. Cut the path as well, or remove it.";
        check.appendChild(el("p", {}, "Uncut circuit, statevector: <" + observable + "> = " + num(exact, 6)));
        shotStats.textContent = "";
        strip.textContent = "";
        termsBody.textContent = "";
      } else {
        var prepared = E.prepare(circuit, cuts, observable);
        var widths = plan.widths.length ? plan.widths.join(" + ") : "0";
        stat(stats, String(plan.fragments.length), plan.fragments.length === 1 ? "fragment" : "fragments, of " + widths + " qubits");
        stat(stats, count(plan.subcircuits), plan.configurations + " configurations × " + plan.fragments.length + " fragment" + (plan.fragments.length === 1 ? "" : "s"));
        if (plan.optimized !== null && cuts.length) stat(stats, count(plan.optimized), "subcircuits with optimized=True");
        stat(stats, gammaText(plan.gamma), "gamma, the sampling overhead");
        check.appendChild(el("p", {}, "Uncut circuit, statevector: <" + observable + "> = " + num(exact, 6)));
        check.appendChild(el("p", {}, "Reconstructed from the fragments: " + num(prepared.value, 6)));
        check.appendChild(el("p", {}, "Difference " + small(Math.abs(prepared.value - exact)) + ", from " + prepared.runs + " distinct fragment circuits computed in your browser."));
        var goldenCount = Object.keys(plan.golden).length, notes = [];
        if (goldenCount) notes.push(goldenCount + " wire cut" + (goldenCount > 1 ? "s are" : " is") + " golden: the state that sends it is real, so its Y terms are zero and are left out, with 6 terms and gamma 3.");
        if (plan.joint && plan.joint.shared) notes.push("Qdislib's default, joint=True, would share ancillas between the " + plan.joint.shared + " Hadamard cuts across one boundary: " + count(plan.joint.subcircuits) + " subcircuits and gamma " + gammaText(plan.joint.gamma) + ", against " + count(plan.subcircuits) + " and " + gammaText(plan.gamma) + " here. The value is the same.");
        note.textContent = notes.join(" ");

        /* Shots */
        var shots = Math.pow(2, Number(shotSlider.value));
        var rng = E.mulberry(state.seed * 7919 + shots);
        var draws = [];
        for (var i = 0; i < 30; i++) draws.push(prepared.draw(shots, rng));
        var run = draws[0], mean = draws.reduce(function (a, b) { return a + b; }, 0) / draws.length;
        var spread = Math.sqrt(draws.reduce(function (a, b) { return a + (b - mean) * (b - mean); }, 0) / (draws.length - 1));
        shotStats.textContent = "";
        stat(shotStats, count(shots), "shots per subcircuit");
        stat(shotStats, num(run, 4), "one run, exact " + num(exact, 4));
        stat(shotStats, spread.toFixed(4), "spread of 30 runs");
        stat(shotStats, count(shots * plan.subcircuits), "shots in total, per run");
        strip.innerHTML = stripPlot(draws, exact);

        /* The terms */
        termsBody.textContent = "";
        if (!cuts.length) termsBody.appendChild(el("p", { class: "qd-note" }, "No cuts yet."));
        cuts.forEach(function (cut, index) {
          var title = cut.kind === "wire"
            ? "Wire cut on q" + cut.qubit + (plan.golden[index] ? " (golden)" : "")
            : (cut.kind === "gate" ? "Gate cut of " : "Hadamard cut of ") + shortName(plain[cut.op]) + "_" + (cut.op + 1);
          termsBody.appendChild(el("p", {}, title));
          var rows = E.terms(plan, cut, index).map(function (t) {
            var a, b;
            if (cut.kind === "wire") { a = t.read === "I" ? "read nothing" : "measure " + t.read; b = "prepare |" + t.prep + "⟩"; }
            else if (cut.kind === "gate") { a = termName(t.a); b = termName(t.b); }
            else { a = hadName(t.a.had); b = hadName(t.b.had); }
            return "<tr><td>" + (t.w >= 0 ? "+" : "−") + Math.abs(t.w).toFixed(4) + "</td><td>" + esc(a) + "</td><td>" + esc(b) + "</td></tr>";
          }).join("");
          var table = termsBody.appendChild(el("table", { class: "qdx-t" }));
          table.innerHTML = "<thead><tr><th>Weight</th><th>Side 1</th><th>Side 2</th></tr></thead><tbody>" + rows + "</tbody>";
        });
      }

      /* The code */
      codeBody.textContent = pythonFor(ops, cuts, observable);
    }
    function termName(a) {
      return { i: "nothing", z: "Z", m: "measure, sign", "r+": "RZ(π/2)", "r-": "RZ(−π/2)" }[a];
    }
    function hadName(h) {
      return h.i === h.j ? (h.i ? "Z" : "nothing") : "ancilla, " + (h.g ? "imaginary" : "real") + " part";
    }
    function pythonFor(ops, cuts, observable) {
      var lines = ["from qiskit import QuantumCircuit", "import Qdislib as qd", "", "qc = QuantumCircuit(" + state.n + ")"];
      ops.forEach(function (op) {
        var args = (op.theta !== undefined ? [String(op.theta)] : []).concat(op.qubits);
        lines.push("qc." + op.name + "(" + args.join(", ") + ")");
      });
      var name = function (i) { return shortName(ops[i]) + "_" + (i + 1); };
      var parts = { gate: [], hadamard: [], wire: [] };
      cuts.forEach(function (c) {
        if (c.kind === "wire") parts.wire.push('("' + name(c.from) + '", "' + name(c.to) + '")');
        else parts[c.kind].push('"' + name(c.op) + '"');
      });
      lines.push("");
      lines.push('cuts = {"gate": [' + parts.gate.join(", ") + '], "hadamard": [' + parts.hadamard.join(", ") + '], "wire": [' + parts.wire.join(", ") + "]}");
      lines.push("print(qd.cutting_overhead(qc, cuts, joint=False))");
      lines.push('print(qd.circuit_cutting(qc, cuts, "' + observable + '", shots=None, joint=False))');
      lines.push('print(qd.analytical_solution(qc, "' + observable + '"))');
      return lines.join("\n");
    }
    function stripPlot(draws, exact) {
      var W = 340, H = 64, L = 8, R = 8;
      var extent = Math.max.apply(null, draws.map(function (d) { return Math.abs(d - exact); }).concat([1e-3]));
      var x = function (v) { return L + ((v - exact + extent * 1.15) / (2.3 * extent)) * (W - L - R); };
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" style="max-width:560px" role="img" aria-label="30 runs against the exact value">';
      s += '<line x1="' + x(exact) + '" x2="' + x(exact) + '" y1="6" y2="42" stroke="var(--qd-b)" stroke-width="2"/>';
      draws.forEach(function (d, i) { s += '<circle cx="' + x(d) + '" cy="' + (12 + (i % 5) * 6) + '" r="3.2" fill="var(--qd-a)" opacity="0.8"/>'; });
      s += '<text x="' + x(exact) + '" y="56" text-anchor="middle">exact ' + num(exact, 3) + "</text>";
      s += '<text x="' + L + '" y="56">' + num(exact - extent * 1.15, 3) + "</text>";
      s += '<text x="' + (W - R) + '" y="56" text-anchor="end">' + num(exact + extent * 1.15, 3) + "</text>";
      return s + "</svg>";
    }

    boardHost.addEventListener("click", function (event) {
      var t = event.target;
      if (t.dataset.obs !== undefined) {
        var q = Number(t.dataset.obs), order = "ZXYI", now = state.observable.charAt(q);
        var next = order.charAt((order.indexOf(now) + 1) % 4);
        state.observable = state.observable.slice(0, q) + next + state.observable.slice(q + 1);
        render();
      } else if (t.dataset.edge !== undefined) {
        var ops = ordered();
        toggleWire(wireEdges(ops)[Number(t.dataset.edge)], ops);
        render();
      } else if (t.dataset.q !== undefined) {
        tapCell(Number(t.dataset.q), Number(t.dataset.c));
      }
    });
    shotSlider.addEventListener("input", render);
    again.addEventListener("click", function () { state.seed += 1; render(); });
    presets().chain();
    render();
  }

  /* ---- What lowers the cost ---- */
  function lowerSection(root) {
    var data = D.lower, k = 3, angle = 1;
    root.textContent = "";
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-low" }));
    var slider = row.appendChild(el("input", { id: "qdx-low", type: "range", min: "1", max: String(data.cuts.length), step: "1", value: String(k), "aria-label": "Number of cuts" }));
    var body = root.appendChild(el("div"));
    root.appendChild(el("span", { class: "qd-label" }, "A gate cut follows the gate's angle"));
    var row2 = root.appendChild(el("div", { class: "qdx-row" }));
    var label2 = row2.appendChild(el("label", { for: "qdx-ang" }));
    var slider2 = row2.appendChild(el("input", { id: "qdx-ang", type: "range", min: "0", max: String(data.angles.length - 1), step: "1", value: String(angle), "aria-label": "Angle of the rzz gates that are cut" }));
    var angles = root.appendChild(el("div", { class: "qd-stats" }));
    root.appendChild(el("p", { class: "qd-note" }, "Every count and gamma is stored from qd.cutting_overhead. The circuits run with shots=None up to 3 cuts, and differ from the uncut value by " + small(data.largest_difference) + " at most. Golden wires are on a chain of fragments, as each has to send one wire alone; the rest on two blocks joined by k gates."));
    function update() {
      label.textContent = "Cuts: " + k;
      body.textContent = "";
      data.groups.forEach(function (group) {
        body.appendChild(el("span", { class: "qd-label" }, group.label));
        body.appendChild(el("p", { class: "qd-note" }, group.note));
        var most = Math.max.apply(null, group.rows.map(function (r) { return r.values[k - 1].subcircuits; }));
        group.rows.forEach(function (r) {
          var v = r.values[k - 1], line = body.appendChild(el("div", { class: "qdx-bar" }));
          line.appendChild(el("span", {}, r.label));
          var track = line.appendChild(el("div"));
          track.appendChild(el("div", { class: "qd-bar", style: "background:var(--qd-a);height:12px;width:" + Math.max(1.5, 100 * Math.log10(v.subcircuits + 1) / Math.log10(most + 1)) + "%" }));
          line.appendChild(el("span", { class: "qdx-num" }, count(v.subcircuits) + " circuits, gamma " + gammaText(v.gamma)));
        });
      });
      var a = data.angles[angle];
      label2.textContent = k + " rzz(" + a.phi + ") gates cut";
      angles.textContent = "";
      stat(angles, gammaText(a.gamma[k - 1]), "gamma, gate cuts");
      stat(angles, count(Math.pow(a.gamma[k - 1], 2)), "shots factor, gamma squared");
      stat(angles, gammaText(a.joint[k - 1]), "gamma, joint Hadamard cuts");
    }
    slider.addEventListener("input", function () { k = Number(slider.value); update(); });
    slider2.addEventListener("input", function () { angle = Number(slider2.value); update(); });
    update();
  }

  /* ---- 2. find_cut ---- */
  function findCutSection(root) {
    var data = D.find_cut, family = data[0].key, kind = "mixed", max = 4;
    root.textContent = "";
    root.appendChild(el("span", { class: "qd-label" }, "A circuit"));
    var familyBox = root.appendChild(el("div", { class: "qd-choices qd-methods" }));
    root.appendChild(el("span", { class: "qd-label" }, "What find_cut may use"));
    var kindBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-max" }));
    var slider = row.appendChild(el("input", { id: "qdx-max", type: "range", step: "1", "aria-label": "max_qubits" }));
    var host = root.appendChild(el("div"));
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var code = root.appendChild(el("pre"));
    var note = root.appendChild(el("p", { class: "qd-note" }));

    function update() {
      var f = data.filter(function (x) { return x.key === family; })[0];
      toggles(familyBox, data.map(function (x) { return { value: x.key, label: x.label }; }), family, function (v) { family = v; update(); });
      toggles(kindBox, [{ value: "gate", label: "Gate" }, { value: "wire", label: "Wire" }, { value: "mixed", label: "Any, by cost" }], kind, function (v) { kind = v; update(); });
      slider.min = "2"; slider.max = String(f.n - 1);
      max = Math.min(max, f.n - 1);
      slider.value = String(max);
      label.textContent = "max_qubits = " + max;
      var run = f.runs[kind + ":" + max];
      var cuts = run.cuts || [];
      var plan = E.plan({ n: f.n, ops: f.ops }, cuts.map(function (c) { return Object.assign({}, c); }), new Array(f.n + 1).join("Z"));
      var cols = columns(f.n, f.ops);
      board(host, drawCircuit({ n: f.n, ops: f.ops, cols: cols, cuts: cuts, plan: plan, width: Math.max.apply(null, cols) + 1, aria: f.label }), true);
      stats.textContent = "";
      if (cuts.length) {
        stat(stats, String(cuts.length), cuts.length === 1 ? "cut" : "cuts");
        stat(stats, String(run.fragments), "fragments");
        stat(stats, count(run.subcircuits), "subcircuits");
        stat(stats, gammaText(run.gamma), "gamma");
      } else {
        stat(stats, "none", "find_cut found no cut within this limit");
      }
      var options = kind === "mixed" ? "mixed=True" : kind === "gate" ? "wire_cut=False" : "gate_cut=False";
      code.textContent = "qd.find_cut(circuit, max_qubits=" + max + ", " + options + ", seed=0)\n# " + (cuts.length ? JSON.stringify(run.spec) : "[]");
      note.textContent = run.difference !== undefined
        ? "Reconstructed with shots=None, this cut differs from the uncut value by " + small(run.difference) + ". Stored from a real run."
        : cuts.length ? "Stored from a real run. Too many configurations to reconstruct here." : "Stored from a real run.";
    }
    slider.addEventListener("input", function () { max = Number(slider.value); update(); });
    update();
  }

  /* ---- 3. Planner ---- */
  function plannerSection(root) {
    var data = D.planner, k = 3, exp = 2;
    root.textContent = "";
    var row1 = root.appendChild(el("div", { class: "qdx-row" }));
    var label1 = row1.appendChild(el("label", { for: "qdx-k" }));
    var s1 = row1.appendChild(el("input", { id: "qdx-k", type: "range", min: "1", max: String(data.cuts.length), step: "1", value: String(k), "aria-label": "Number of cuts" }));
    var row2 = root.appendChild(el("div", { class: "qdx-row" }));
    var label2 = row2.appendChild(el("label", { for: "qdx-eps" }));
    var s2 = row2.appendChild(el("input", { id: "qdx-eps", type: "range", min: "1", max: "3", step: "0.25", value: String(exp), "aria-label": "Target precision" }));
    var table = root.appendChild(el("table", { class: "qdx-t" }));
    root.appendChild(el("p", { class: "qd-note" },
      "Stored from qd.cutting_overhead and qd.shots_for, with delta = " + data.delta + ". The shots are Hoeffding's worst case for any circuit with this cut, usually far more than a run needs: see the error bars below."));
    function update() {
      var epsilon = Math.pow(10, -exp);
      label1.textContent = "Cuts: " + k;
      label2.textContent = "Target error: " + (epsilon < 0.01 ? epsilon.toFixed(4) : epsilon.toFixed(3));
      var scale = Math.pow(data.epsilon / epsilon, 2);
      var html = "<thead><tr><th>Cut</th><th>Circuits</th><th>gamma</th><th>Shots, worst case</th></tr></thead><tbody>";
      Object.keys(data.kinds).forEach(function (key) {
        var kind = data.kinds[key], r = kind.rows[k - 1];
        html += "<tr><td>" + esc(kind.label) + "</td><td>" + count(r.subcircuits) + "</td><td>" + gammaText(r.gamma) + "</td><td>" + count(r.shots * scale) + "</td></tr>";
      });
      table.innerHTML = html + "</tbody>";
    }
    s1.addEventListener("input", function () { k = Number(s1.value); update(); });
    s2.addEventListener("input", function () { exp = Number(s2.value); update(); });
    update();
  }

  /* ---- 4. Distributions ---- */
  function bars(host, keys, exact, values, labels) {
    host.textContent = "";
    var legend = host.appendChild(el("div", { class: "qd-legend" }));
    [["var(--qd-a)", labels[0]], ["var(--qd-b)", labels[1]]].forEach(function (x) {
      var item = legend.appendChild(el("span"));
      item.appendChild(el("span", { class: "qd-swatch", style: "background:" + x[0] }));
      item.appendChild(document.createTextNode(x[1]));
    });
    var scale = Math.max.apply(null, exact.concat(values).concat([1e-9]));
    keys.forEach(function (key, i) {
      var row = host.appendChild(el("div", { class: "qd-row" }));
      row.appendChild(el("code", {}, key));
      var cell = row.appendChild(el("div"));
      [["var(--qd-a)", exact[i]], ["var(--qd-b)", values[i]]].forEach(function (x) {
        cell.appendChild(el("div", { class: "qd-bar", style: "background:" + x[0] + ";width:" + Math.max(0, 100 * x[1] / scale) + "%" }));
      });
      row.appendChild(el("span", {}, values[i].toFixed(3)));
    });
  }
  function distributionSection(root) {
    var data = D.distributions, tab = "full", level = 2, seed = 0, a = 0, b = 7, pair = 10;
    root.textContent = "";
    var tabs = root.appendChild(el("div", { class: "qdx-tabs" }));
    var body = root.appendChild(el("div"));
    function update() {
      toggles(tabs, [{ value: "full", label: "Full distribution" }, { value: "marginal", label: "qubits= marginal" }, { value: "big", label: "40 qubits" }], tab, function (v) { tab = v; update(); });
      body.textContent = "";
      if (tab === "full") {
        var f = data.full;
        var row = body.appendChild(el("div", { class: "qd-shots" }));
        var slider = row.appendChild(el("input", { type: "range", min: "0", max: String(f.shots.length - 1), step: "1", value: String(level), "aria-label": "Shots per subcircuit" }));
        var again = row.appendChild(el("button", { type: "button" }, "Another run (" + (seed + 1) + " of 3)"));
        var run = f.runs[level][seed];
        var stats = body.appendChild(el("div", { class: "qd-stats" }));
        stat(stats, count(f.shots[level]), "shots per subcircuit");
        stat(stats, run.tvd.toFixed(4), "total variation distance from exact");
        bars(body.appendChild(el("div")), f.keys, f.exact, run.top, ["exact", "circuit_cutting_sampling"]);
        body.appendChild(el("p", { class: "qd-note" }, "The 12 likeliest of 256 outcomes, 8 qubits with one gate cut, " + data.cut + ". Stored from real runs."));
        slider.addEventListener("input", function () { level = Number(slider.value); seed = 0; update(); });
        again.addEventListener("click", function () { seed = (seed + 1) % 3; update(); });
      } else if (tab === "marginal") {
        body.appendChild(el("span", { class: "qd-label" }, "Keep two qubits"));
        var pick = body.appendChild(el("div", { class: "qdx-tools" }));
        for (var q = 0; q < 8; q++) (function (q) {
          var on = q === a || q === b;
          var btn = pick.appendChild(el("button", { type: "button", "aria-pressed": String(on) }, "q" + q));
          btn.addEventListener("click", function () {
            if (q === a || q === b) return;
            a = b; b = q;
            update();
          });
        })(q);
        var lo = Math.min(a, b), hi = Math.max(a, b), m = data.marginals[lo + "," + hi];
        var code = body.appendChild(el("pre"));
        code.textContent = 'qd.circuit_cutting_sampling(circuit, ["' + data.cut + '"], shots=4096, seed=1, qubits=[' + lo + ", " + hi + "])";
        bars(body.appendChild(el("div")), ["00", "01", "10", "11"], m.exact, m.run, ["exact", "qubits=[" + lo + ", " + hi + "]"]);
        body.appendChild(el("p", { class: "qd-note" }, "Each fragment circuit returns only what the marginal reads, so the cost follows the kept qubits. Bit j of a key is qubits[j]. Stored from real runs."));
      } else {
        var big = data.factored;
        var s = body.appendChild(el("div", { class: "qd-stats" }));
        stat(s, "2^40", "outcomes, which no dict can hold");
        stat(s, big.text.match(/([\d,]+) stored/)[1], "values stored by the FactoredDistribution");
        stat(s, gammaText(big.gamma), "gamma, two gate cuts");
        var r = body.appendChild(el("div", { class: "qdx-row" }));
        var lab = r.appendChild(el("label", { for: "qdx-pair" }));
        var sl = r.appendChild(el("input", { id: "qdx-pair", type: "range", min: "0", max: "38", step: "1", value: String(pair), "aria-label": "Which neighbouring qubits to read" }));
        lab.textContent = "d.marginal([" + pair + ", " + (pair + 1) + "])";
        var values = big.pairs[pair + "," + (pair + 1)];
        var host = body.appendChild(el("div"));
        bars(host, ["00", "01", "10", "11"], values, values, ["contracted when read", "the same"]);
        host.querySelector(".qd-legend").remove();
        sl.addEventListener("input", function () { pair = Number(sl.value); update(); });
        body.appendChild(el("span", { class: "qd-label" }, "d.most_likely(5)"));
        var pre = body.appendChild(el("pre"));
        pre.textContent = Object.keys(big.most_likely).map(function (k) { return k + "  " + big.most_likely[k].toExponential(2); }).join("\n");
        body.appendChild(el("p", { class: "qd-note" }, big.text + ". Exact, in " + big.seconds + " s on one machine. Stored from a real run; no circuit runs when it is read. On qubits 2–7, its marginals differ from a 10-qubit statevector of the same gates by " + small(big.check) + "."));
      }
    }
    update();
  }

  /* ---- 5. Error bars ---- */
  function errorSection(root) {
    var data = D.errors, key = "gate", level = 1;
    root.textContent = "";
    var tabs = root.appendChild(el("div", { class: "qdx-tabs" }));
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-eshots" }));
    var slider = row.appendChild(el("input", { id: "qdx-eshots", type: "range", min: "0", max: String(data.shots.length - 1), step: "1", value: String(level), "aria-label": "Shots per subcircuit" }));
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var plot = root.appendChild(el("div"));
    root.appendChild(el("p", { class: "qd-note" }, "Each dot is a separate run, with its own seed. The band is the ±1 standard error that one run reported with error=True, from a bootstrap over its shots, at no device cost. Stored from real runs, 24 per point."));
    function update() {
      toggles(tabs, Object.keys(data.cases).map(function (k) { return { value: k, label: data.cases[k].label }; }), key, function (v) { key = v; update(); });
      var r = data.cases[key].rows[level];
      label.textContent = count(data.shots[level]) + " shots per subcircuit";
      stats.textContent = "";
      stat(stats, r.reported.toFixed(4), "error=True, from one run");
      stat(stats, r.spread.toFixed(4), "spread of 24 runs");
      stat(stats, count(r.total), "shots in one run");
      var W = 340, H = 92, L = 8, R = 8, extent = Math.max.apply(null, r.values.map(function (v) { return Math.abs(v - data.exact); }).concat([r.reported * 1.3]));
      var x = function (v) { return L + ((v - data.exact + extent * 1.15) / (2.3 * extent)) * (W - L - R); };
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" style="max-width:560px" role="img" aria-label="Runs against the reported error">';
      s += '<rect x="' + x(r.value - r.reported) + '" y="6" width="' + (x(r.value + r.reported) - x(r.value - r.reported)) + '" height="60" fill="var(--qd-b)" opacity="0.2"/>';
      s += '<line x1="' + x(data.exact) + '" x2="' + x(data.exact) + '" y1="4" y2="68" stroke="var(--qd-b)" stroke-width="2"/>';
      r.values.forEach(function (v, i) { s += '<circle cx="' + x(v) + '" cy="' + (12 + (i % 8) * 7) + '" r="3.2" fill="var(--qd-a)"/>'; });
      s += '<text x="' + x(data.exact) + '" y="82" text-anchor="middle">exact ' + num(data.exact, 3) + "</text>";
      s += '<text x="' + L + '" y="82">' + num(data.exact - extent * 1.15, 3) + "</text>";
      s += '<text x="' + (W - R) + '" y="82" text-anchor="end">' + num(data.exact + extent * 1.15, 3) + "</text>";
      plot.innerHTML = s + "</svg>";
    }
    slider.addEventListener("input", function () { level = Number(slider.value); update(); });
    update();
  }

  /* ---- 6. Noise ---- */
  function noiseSection(root) {
    var data = D.noise, level = 2;
    root.textContent = "";
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-noise-p" }));
    var slider = row.appendChild(el("input", { id: "qdx-noise-p", type: "range", min: "0", max: String(data.rows.length - 1), step: "1", value: String(level), "aria-label": "Two-qubit depolarizing error" }));
    var host = root.appendChild(el("div"));
    root.appendChild(el("pre", {}, 'noise = NoiseModel()\nnoise.add_all_qubit_quantum_error(depolarizing_error(p, 2), ["cz", "cx"])\nbackend = qd.build_backend("AER_BACKEND", qd.BackendConfig(extra={"noise_model": noise}))\nqd.gate_cutting(circuit, cut, "ZZZZZZ", backend=backend, shots=' + data.shots + ", seed=1, error=True)"));
    root.appendChild(el("p", { class: "qd-note" }, "A simulation on 6 qubits, cut at 2 CZ gates into two fragments of 3, with a depolarizing error on two-qubit gates only and " + count(data.shots) + " shots per circuit, stored from real runs. It has no readout or single-qubit noise and no idle time, so it says nothing about any device. Hadamard ancillas add two-qubit gates of their own."));
    function update() {
      var r = data.rows[level];
      label.textContent = "Two-qubit error p = " + r.p;
      var rows = [["Uncut", r.uncut, null], ["Gate cuts", r.gate[0], r.gate[1]], ["Hadamard cuts", r.hadamard[0], r.hadamard[1]]];
      var html = '<table class="qdx-t"><thead><tr><th></th><th>&lt;ZZZZZZ&gt;</th><th>Error bar</th><th>Off by</th></tr></thead><tbody>';
      rows.forEach(function (x) {
        html += "<tr><td>" + x[0] + "</td><td>" + num(x[1], 4) + "</td><td>" + (x[2] === null ? "" : "±" + x[2].toFixed(4)) + "</td><td>" + Math.abs(x[1] - data.exact).toFixed(4) + "</td></tr>";
      });
      host.innerHTML = html + "<tr><td>Exact</td><td>" + num(data.exact, 4) + "</td><td></td><td></td></tr></tbody></table>";
    }
    slider.addEventListener("input", function () { level = Number(slider.value); update(); });
    update();
  }

  /* ---- 7. SDKs ---- */
  function sdkSection(root) {
    var data = D.sdk, active = data[0].name;
    root.textContent = "";
    var tabs = root.appendChild(el("div", { class: "qdx-tabs" }));
    var body = root.appendChild(el("div"));
    function update() {
      toggles(tabs, data.map(function (s) { return { value: s.name, label: s.name }; }), active, function (v) { active = v; update(); });
      var s = data.filter(function (x) { return x.name === active; })[0];
      body.textContent = "";
      body.appendChild(el("pre", {}, s.install));
      body.appendChild(el("pre", {}, s.code));
      var stats = body.appendChild(el("div", { class: "qd-stats" }));
      stat(stats, s.cut, "the cut find_cut chose");
      stat(stats, num(s.value, 4), "cut, 8192 shots on " + s.name + "'s simulator");
      stat(stats, num(s.exact, 4), "exact");
    }
    root.appendChild(el("p", { class: "qd-note" }, "Each snippet ran as shown, as its own process, with Qdislib " + D.version + "."));
    update();
  }

  /* ---- 8. pack ---- */
  function packSection(root) {
    var data = D.pack, buffer = "1", job = 0;
    root.textContent = "";
    root.appendChild(el("span", { class: "qd-label" }, "Idle qubits between groups: buffer="));
    var bufferBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-job" }));
    var slider = row.appendChild(el("input", { id: "qdx-job", type: "range", min: "0", step: "1", "aria-label": "Which job" }));
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var map = root.appendChild(el("div", { class: "qdx-map" }));
    var chipBox = root.appendChild(el("div", { class: "qd-legend" }));
    root.appendChild(el("p", { class: "qd-note" }, "qd.pack on a fake device, " + data.device + ", a published calibration snapshot of " + data.qubits + " qubits: the " + data.circuits + " subcircuits of a 6-qubit circuit with one cut. Each colour is one subcircuit's group. Placement only, stored from a real call; nothing ran on any device."));
    function update() {
      var run = data.runs[buffer];
      toggles(bufferBox, Object.keys(data.runs).map(function (b) { return { value: b, label: b }; }), buffer, function (v) { buffer = v; job = 0; update(); });
      slider.max = String(run.jobs.length - 1);
      job = Math.min(job, run.jobs.length - 1);
      slider.value = String(job);
      label.textContent = "Job " + (job + 1) + " of " + run.jobs.length;
      stats.textContent = "";
      stat(stats, String(run.jobs.length), "jobs for " + data.circuits + " subcircuits");
      stat(stats, String(run.regions.length), "groups of qubits on the device");
      stat(stats, String(run.jobs[job].length), "subcircuits in this job");
      var xs = [], ys = [];
      Object.keys(data.points).forEach(function (q) { xs.push(data.points[q][0]); ys.push(data.points[q][1]); });
      var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs), minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
      var W = 340, H = 300, P = 16;
      var px = function (q) { var p = data.points[q]; return [P + (p[0] - minX) / (maxX - minX) * (W - 2 * P), P + (p[1] - minY) / (maxY - minY) * (H - 2 * P)]; };
      var owner = {};
      run.jobs[job].forEach(function (slot, i) { slot.qubits.forEach(function (q) { owner[q] = i; }); });
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" style="max-width:520px;margin:auto" role="img" aria-label="The device, with this job\'s groups of qubits">';
      data.edges.forEach(function (e) {
        var a = px(e[0]), b = px(e[1]);
        s += '<line x1="' + a[0] + '" y1="' + a[1] + '" x2="' + b[0] + '" y2="' + b[1] + '" stroke="var(--qd-line)" stroke-width="2"/>';
      });
      Object.keys(data.points).forEach(function (q) {
        var p = px(q), on = owner[q] !== undefined;
        s += '<circle cx="' + p[0] + '" cy="' + p[1] + '" r="9.5" fill="' + (on ? fragColor(owner[q]) : "var(--color-background-primary, #fff)") + '" stroke="' + (on ? "none" : "var(--color-foreground-secondary, #777)") + '"/>';
        s += '<text x="' + p[0] + '" y="' + (p[1] + 3) + '" text-anchor="middle" style="fill:' + (on ? "#fff" : "var(--color-foreground-secondary, #555)") + '">' + q + "</text>";
      });
      map.innerHTML = s + "</svg>";
      chipBox.textContent = "";
      run.jobs[job].forEach(function (slot, i) {
        var item = chipBox.appendChild(el("span"));
        item.appendChild(el("span", { class: "qd-swatch", style: "background:" + fragColor(i) }));
        item.appendChild(document.createTextNode("subcircuit " + slot.circuit + ": " + slot.qubits.join(", ")));
      });
    }
    slider.addEventListener("input", function () { job = Number(slider.value); update(); });
    update();
  }

  var SECTIONS = {
    "qdx-play": playground, "qdx-low": lowerSection, "qdx-find": findCutSection, "qdx-plan": plannerSection, "qdx-dist": distributionSection,
    "qdx-err": errorSection, "qdx-noise": noiseSection, "qdx-sdk": sdkSection, "qdx-pack": packSection
  };
  Object.keys(SECTIONS).forEach(function (id) {
    var root = document.getElementById(id);
    if (root) SECTIONS[id](root);
  });
})();
