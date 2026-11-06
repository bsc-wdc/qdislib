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
  function num(x, digits) { return (Number(x.toFixed(digits)) < 0 ? "−" : "+") + Math.abs(x).toFixed(digits); }
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
  /* A press shows on every button, also on a touch screen, where :active barely lasts. */
  document.addEventListener("click", function (event) {
    var b = event.target.closest && event.target.closest(".qd-try button:not(.copybtn)");
    if (!b) return;
    b.classList.add("qdx-flash");
    setTimeout(function () { b.classList.remove("qdx-flash"); }, 220);
  });
  /* A code cell with the docs' own copy button (sphinx-copybutton): its markup and icon, and
     its ClipboardJS, which listens on the whole page and so also serves cells drawn later.
     Not a div: that script gives a button to each div.highlight it finds when the page
     loads, and these would get a second one. */
  var COPY_ICON = '<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-copy" width="44" height="44" viewBox="0 0 24 24" stroke-width="1.5" stroke="#000000" fill="none" stroke-linecap="round" stroke-linejoin="round"><title>Copy to clipboard</title><path stroke="none" d="M0 0h24v24H0z" fill="none"/><rect x="8" y="8" width="12" height="12" rx="2"/><path d="M16 8v-2a2 2 0 0 0 -2 -2h-8a2 2 0 0 0 -2 2v8a2 2 0 0 0 2 2h2"/></svg>';
  var cells = 0;
  function codeCell(parent, text) {
    var box = parent.appendChild(el("section", { class: "highlight qdx-code" }));
    var id = "qdx-code-" + (++cells);
    var pre = box.appendChild(el("pre", { id: id }, text));
    var button = box.appendChild(el("button", { type: "button", class: "copybtn o-tooltip--left", "data-tooltip": "Copy", "data-clipboard-target": "#" + id }));
    button.innerHTML = COPY_ICON;
    return pre;
  }
  /* The page's dropdown: the circuit picker's look, for every list of things to pick. */
  function picker(root, id, label, items, value, onChange) {
    var row = root.appendChild(el("div", { class: "qdx-start" }));
    row.appendChild(el("label", { for: id }, label));
    var menu = row.appendChild(el("select", { id: id }));
    items.forEach(function (item) { menu.appendChild(el("option", { value: item.value }, item.label)); });
    menu.value = value;
    menu.addEventListener("change", function () { onChange(menu.value); });
    return menu;
  }
  function fragColor(f) { return "var(--qdx-f" + (f % 8) + ")"; }
  function shortName(op) { return op.name.toUpperCase(); }

  /* Scissors as a drawing: a glyph can come out as a coloured emoji on a phone. */
  function scissors(x, y, colour, size, key) {
    var u = (size || 12) / 12;
    return '<g transform="translate(' + x + " " + y + ") scale(" + u + ')"><g' + (key ? ' data-k="' + key + '"' : "") + ' fill="none" stroke="' + colour + '" stroke-width="1.9" stroke-linecap="round">' +
      '<circle cx="0" cy="0" r="7" fill="var(--color-background-primary, #fff)" stroke="none"/>' +
      '<circle cx="-5" cy="-3.5" r="2.4"/><circle cx="-5" cy="3.5" r="2.4"/>' +
      '<path d="M-3 -2.2 L6 4 M-3 2.2 L6 -4"/></g></g>';
  }
  function cutColour(kind) { return kind === "hadamard" ? "var(--qdx-had)" : kind === "gate" ? "var(--qdx-gate)" : "var(--qdx-wire)"; }
  function legendHtml() {
    return [["gate", "gate cut"], ["hadamard", "Hadamard cut"], ["wire", "wire cut"]].map(function (x) {
      return '<span><svg width="22" height="18" viewBox="-11 -9 22 18" style="display:inline-block;width:22px;vertical-align:-4px">' + scissors(0, 0, cutColour(x[0]), 15) + "</svg> " + x[1] + "</span>";
    }).join("");
  }

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
    var width = o.width, W = LEFT + width * CELL + (o.ends ? 50 : 12), H = TOP + o.n * ROW + 6, lineEnd = o.ends ? W - 38 : W - 6;
    var y = function (q) { return TOP + q * ROW + ROW / 2; };
    var x = function (c) { return LEFT + c * CELL + CELL / 2; };
    var s = '<svg viewBox="0 0 ' + W + " " + H + '" width="' + W + '" height="' + H + '" style="max-width:' + (o.keys ? Math.round(W * 1.35) : W) + "px" + (o.keys ? ";min-width:" + Math.round(W * 0.7) + "px" : "") + '" role="img" aria-label="' + esc(o.aria || "A circuit") + '">';
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
        s += '<line x1="' + start + '" x2="' + lineEnd + '" y1="' + y(q) + '" y2="' + y(q) + '" stroke="var(--qd-line)" stroke-width="2"/>';
      }
      segs.forEach(function (g) {
        var end = lineEnd;
        if (g.outgoing) end = wireX[q + ":" + g.outgoing.from];
        var colour = g.fragment === undefined ? "var(--qd-line)" : fragColor(g.fragment);
        s += '<line' + (o.keys ? ' data-k="l' + q + ":" + g.ordinal + '"' : "") + ' x1="' + start + '" x2="' + end + '" y1="' + y(q) + '" y2="' + y(q) + '" stroke="' + colour + '" stroke-width="2.5"/>';
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
        var stroke = cut ? cutColour(cut.kind) : "var(--color-foreground-secondary, #555)";
        s += '<line x1="' + c + '" x2="' + c + '" y1="' + a + '" y2="' + b + '" stroke="' + stroke + '" stroke-width="2"' + (cut ? ' stroke-dasharray="4 3"' : "") + "/>";
        if (cut) s += scissors(c, (a + b) / 2, stroke, 12, o.keys ? "xg" + op.id + cut.kind : null);
      }
      op.qubits.forEach(function (qubit, side) {
        var f = fragOf(side), fill = f === null ? "#667" : fragColor(f);
        var cy = y(qubit);
        s += '<rect' + (o.keys ? ' data-k="g' + op.id + ":" + side + '"' : "") + ' x="' + (c - 17) + '" y="' + (cy - 14) + '" width="34" height="28" rx="6" fill="' + fill + '"' +
          (cut ? ' stroke="' + cutColour(cut.kind) + '" stroke-width="3" stroke-dasharray="4 2"' : "") + "/>";
        s += '<text class="g" x="' + c + '" y="' + (cy + 4) + '" text-anchor="middle">' + esc(shortName(op)) + "</text>";
      });
      s += '<text class="i" x="' + c + '" y="' + (y(Math.min.apply(null, op.qubits)) - 17) + '" text-anchor="middle">' + (i + 1) + "</text>";
    });
    /* wire cuts */
    (o.cuts || []).forEach(function (cut) {
      if (cut.kind !== "wire") return;
      var cx = wireX[cut.qubit + ":" + cut.from];
      s += scissors(cx, y(cut.qubit), cutColour("wire"), 12, o.keys ? "xw" + cut.qubit + ":" + cut.from : null);
    });
    if (o.ends) {
      o.ends.forEach(function (letter, q) {
        s += '<g class="qdx-obs"><rect x="' + (W - 34) + '" y="' + (y(q) - 13) + '" width="26" height="26" rx="6"/><text class="q" x="' + (W - 21) + '" y="' + (y(q) + 4.5) + '" text-anchor="middle">' + esc(letter) + "</text></g>";
      });
    }
    if (o.overlay) s += o.overlay({ x: x, y: y, W: W });
    return s + "</svg>";
  }

  function board(host, svg, fit) {
    host.className = "qdx-board" + (fit ? " qdx-fit" : "");
    host.innerHTML = svg;
  }

  /* ---- 1. The playground ---- */
  function playground(root) {
    var COLS = 10, MAX_CUTS = 4;
    var PALETTE = [
      ["h", "H"], ["x", "X"], ["s", "S"], ["t", "T"], ["rx", "RX"], ["ry", "RY"], ["rz", "RZ"],
      ["cz", "CZ"], ["cx", "CX"], ["rzz", "RZZ"]
    ];
    var ANGLED = { rx: 1, ry: 1, rz: 1, rzz: 1 }, TWO = { cz: 1, cx: 1, rzz: 1 };
    var state = { n: 5, ops: [], cuts: [], observable: "ZZZZZ", tool: null, angle: 0.8, pending: null, seed: 1, said: "", example: "chain", auto: false, editing: false };
    var uid = 1, painted = {}, keyed = null;

    function add(name, qubits, col, theta) {
      var op = { id: uid++, name: name, qubits: qubits, col: col };
      if (ANGLED[name]) op.theta = theta === undefined ? state.angle : theta;
      state.ops.push(op);
      return op;
    }
    function fresh(n) {
      state.n = n;
      state.observable = new Array(n + 1).join("Z");
      state.ops = [];
      state.cuts = [];
      state.pending = null;
      state.said = "";
      state.auto = false;
    }
    var EXAMPLES = {
      chain: ["Chain", function () {
        fresh(5);
        [0.8, 0.5, 0.6, 0.4, 0.3].forEach(function (a, q) { add("ry", [q], 0, a); });
        for (var q = 0; q < 4; q++) add("cz", [q, q + 1], q + 1);
        [0.2, 0.1, 0.7, 0.3, 0.9].forEach(function (a, q) { add("rx", [q], 5 + q % 2, a); });
      }],
      blocks: ["Weak links", function () {
        fresh(6);
        [0, 1, 2, 3, 4, 5].forEach(function (q) { add("h", [q], 0); });
        add("rzz", [0, 1], 1, 1.2); add("rzz", [3, 4], 1, 1.2);
        add("rzz", [1, 2], 2, 1.2); add("rzz", [4, 5], 2, 1.2);
        add("rzz", [2, 3], 3, 0.2); add("rzz", [0, 5], 4, 0.2);
        [0, 1, 2, 3, 4, 5].forEach(function (q) { add("rx", [q], 5, 0.5); });
      }],
      real: ["Real gates", function () {
        fresh(4);
        [0, 1, 2, 3].forEach(function (q) { add("ry", [q], 0, 0.5 + 0.2 * q); });
        add("cz", [0, 1], 1); add("cz", [1, 2], 2); add("cz", [2, 3], 3);
      }],
      ghz: ["GHZ", function () {
        fresh(6);
        add("h", [0], 0);
        for (var q = 0; q < 5; q++) add("cx", [q, q + 1], q + 1);
      }],
      qaoa: ["QAOA ring", function () {
        fresh(6);
        [0, 1, 2, 3, 4, 5].forEach(function (q) { add("h", [q], 0); });
        add("rzz", [0, 1], 1, 0.7); add("rzz", [2, 3], 1, 0.7); add("rzz", [4, 5], 1, 0.7);
        add("rzz", [1, 2], 2, 0.7); add("rzz", [3, 4], 2, 0.7);
        add("rzz", [0, 5], 3, 0.7);
        [0, 1, 2, 3, 4, 5].forEach(function (q) { add("rx", [q], 4, 0.6); });
      }],
      ising: ["Ising steps", function () {
        fresh(5);
        [0, 1, 2, 3, 4].forEach(function (q) { add("ry", [q], 0, 0.4); });
        [1, 4].forEach(function (c) {
          add("rzz", [0, 1], c, 0.5); add("rzz", [2, 3], c, 0.5);
          add("rzz", [1, 2], c + 1, 0.5); add("rzz", [3, 4], c + 1, 0.5);
          [0, 1, 2, 3, 4].forEach(function (q) { add("rx", [q], c + 2, 0.3); });
        });
      }],
      brick: ["Brickwork", function () {
        fresh(6);
        [0, 3, 6].forEach(function (c) {
          [0, 1, 2, 3, 4, 5].forEach(function (q) { add("ry", [q], c, 0.3 + 0.1 * q + 0.05 * c); });
          if (c === 6) return;
          add("cz", [0, 1], c + 1); add("cz", [2, 3], c + 1); add("cz", [4, 5], c + 1);
          add("cz", [1, 2], c + 2); add("cz", [3, 4], c + 2);
        });
      }],
      linked: ["Linked GHZ", function () {
        fresh(6);
        add("h", [0], 0); add("h", [2], 0);
        add("cx", [0, 1], 1); add("cx", [2, 3], 1);
        add("cx", [3, 4], 2);
        add("cx", [4, 5], 3);
        add("cx", [1, 2], 4);
        [0, 1, 2, 3, 4, 5].forEach(function (q) { add("rx", [q], 5, 0.4); });
      }],
      star: ["Star", function () {
        fresh(5);
        [0, 1, 2, 3, 4].forEach(function (q) { add("ry", [q], 0, 0.5 + 0.1 * q); });
        [1, 2, 3, 4].forEach(function (q) { add("cz", [0, q], q); });
        [0, 1, 2, 3, 4].forEach(function (q) { add("rx", [q], 5, 0.3); });
      }],
      clifford: ["Clifford + T", function () {
        fresh(4);
        [0, 1, 2, 3].forEach(function (q) { add("h", [q], 0); add("t", [q], 1); });
        add("cx", [0, 1], 2); add("cx", [2, 3], 2);
        add("cx", [1, 2], 3);
        [0, 1, 2, 3].forEach(function (q) { add("t", [q], 4); add("h", [q], 5); });
        add("cx", [0, 1], 6); add("cx", [2, 3], 6);
        [0, 1, 2, 3].forEach(function (q) { add("t", [q], 7); add("h", [q], 8); });
      }],
      random: ["Random", function () {
        var rng = E.mulberry((state.seed = state.seed + 1)), n = 4 + Math.floor(rng() * 3);
        fresh(n);
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
              add(["cz", "cx", "rzz"][Math.floor(rng() * 3)], [q, r], c, 0.3 + 2.4 * rng());
            } else {
              used[q] = 1;
              add(PALETTE[Math.floor(rng() * 7)][0], [q], c, 0.3 + 2.4 * rng());
            }
          }
        }
      }],
      empty: ["Empty", function () { fresh(4); }]
    };

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
    /* What one tap can cut: a line between two gates next to each other on it. Qdislib names
       a wire cut by its two gates, and that name cuts every line on which they are next to
       each other, so a pair that is neighbours on both lines is two wires, not offered here. */
    function wireEdges(ops) {
      var last = {}, out = [];
      ops.forEach(function (op, i) {
        var before = {};
        Object.keys(last).forEach(function (q) { before[q] = last[q]; });
        op.qubits.forEach(function (q) {
          if (q in before) {
            var lines = ops[before[q]].qubits.filter(function (x) { return op.qubits.indexOf(x) >= 0 && before[x] === before[q]; });
            if (lines.length === 1) out.push({ qubit: q, from: before[q], to: i });
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

    /* ---- The page ---- */
    root.textContent = "";
    var startRow = root.appendChild(el("div", { class: "qdx-start" }));
    startRow.appendChild(el("label", { for: "qdx-circuit" }, "Circuit"));
    var menu = startRow.appendChild(el("select", { id: "qdx-circuit" }));
    Object.keys(EXAMPLES).forEach(function (key) { menu.appendChild(el("option", { value: key }, EXAMPLES[key][0])); });
    var yours = menu.appendChild(el("option", { value: "yours" }, "Your circuit"));
    var editButton = startRow.appendChild(el("button", { type: "button", "aria-pressed": "false" }, "Edit the circuit"));
    editButton.addEventListener("click", function () {
      state.editing = !state.editing;
      state.tool = null;
      state.pending = null;
      state.said = "";
      render();
    });
    var another = startRow.appendChild(el("button", { type: "button" }, "Another random one"));
    function pickCircuit(key) {
      state.example = key;
      EXAMPLES[key][1]();
      state.tool = null;
      state.editing = key === "empty";
      render();
    }
    menu.addEventListener("change", function () { if (menu.value !== "yours") pickCircuit(menu.value); });
    another.addEventListener("click", function () { pickCircuit("random"); });
    function drawExamples() {
      yours.hidden = yours.disabled = state.example !== null;
      menu.value = state.example === null ? "yours" : state.example;
      another.style.display = state.example === "random" ? "" : "none";
    }
    var hint = root.appendChild(el("p", { class: "qdx-hint", "aria-live": "polite" }));
    var editBody = root.appendChild(el("div", { class: "qdx-edit" }));
    var toolBox = editBody.appendChild(el("div", { class: "qdx-tools" }));
    var angleRow = editBody.appendChild(el("div", { class: "qdx-row" }));
    var angleLabel = angleRow.appendChild(el("label", { for: "qdx-angle" }, "Angle of the next rotation: 0.8"));
    var angle = angleRow.appendChild(el("input", { id: "qdx-angle", type: "range", min: "0.1", max: "3", step: "0.1", value: "0.8", "aria-label": "Angle of the next rotation" }));
    angle.addEventListener("input", function () { state.angle = Number(angle.value); angleLabel.textContent = "Angle of the next rotation: " + state.angle.toFixed(1); });

    var boardHost = root.appendChild(el("div", { id: "qdx-play-board" }));
    var legend = root.appendChild(el("div", { class: "qd-legend" }));
    legend.innerHTML = legendHtml();
    root.appendChild(el("p", { class: "qd-note" }, "The letters at the end of the lines are what is measured. Tap one to change it: Z, X, Y, or I for nothing."));
    var actions = root.appendChild(el("div", { class: "qdx-actions" }));
    var auto = actions.appendChild(el("button", { type: "button", class: "qdx-primary" }, "Auto cut"));
    var clear = actions.appendChild(el("button", { type: "button" }, "Clear cuts"));
    var warn = root.appendChild(el("div", { class: "qdx-warn", hidden: "" }));
    var payoff = root.appendChild(el("div", { class: "qdx-payoff" }));

    var shotBody = root.appendChild(el("div", { class: "qdx-withshots" }));
    shotBody.appendChild(el("span", { class: "qd-label" }, "With shots"));
    var shotRow = shotBody.appendChild(el("div", { class: "qd-shots" }));
    var shotSlider = shotRow.appendChild(el("input", { type: "range", min: "6", max: "16", step: "1", value: "10", "aria-label": "Shots per subcircuit" }));
    var again = shotRow.appendChild(el("button", { type: "button" }, "Draw again"));
    var shotStats = shotBody.appendChild(el("div", { class: "qd-stats" }));
    var strip = shotBody.appendChild(el("div"));
    shotBody.appendChild(el("p", { class: "qd-note" },
      "Slide to more shots and the dots close in on the exact value: the axis stays the same, set by the fewest shots. Shots are simulated here: each distinct fragment circuit gives a +1 or −1 per shot, drawn with its exact probability. Dots are 30 independent runs."));

    var how = root.appendChild(el("details"));
    how.appendChild(el("summary", {}, "How it works"));
    how.appendChild(el("div", {})).innerHTML =
      "<p>A cut replaces part of a circuit with local operations, so that it falls apart into <b>fragments</b>, each in its own colour, which a smaller device can run.</p>" +
      "<ul><li>A <b>gate cut</b> replaces a two-qubit gate with single-qubit operations and a measurement on each side. Tap the gate.</li>" +
      "<li>A <b>Hadamard cut</b> does it with one extra qubit per side instead of the measurement. Tap the gate again.</li>" +
      "<li>A <b>wire cut</b> measures a qubit's line at one point and prepares it again later. Tap the circle on the line.</li></ul>" +
      "<p>Each fragment runs several times, once per term of its cuts. Those are the <b>subcircuits</b>. Adding up their results with the right weights gives back the uncut result. The price is <b>gamma</b>: for the same worst-case precision, the shots grow as gamma squared.</p>";

    var termsBox = root.appendChild(el("details"));
    termsBox.appendChild(el("summary", {}, "The terms of each cut"));
    var termsBody = termsBox.appendChild(el("div"));
    var codeSection = root.appendChild(el("details"));
    codeSection.appendChild(el("summary", {}, "The same thing with Qdislib"));
    var codeBody = codeCell(codeSection);

    function nameOf(op, ops) { return shortName(op) + "_" + (ops.indexOf(op) + 1); }

    function drawTools() {
      var items = PALETTE.map(function (p) { return { value: p[0], label: p[1] }; }).concat([{ value: "erase", label: "Erase" }]);
      toggles(toolBox, items, state.tool, function (v) {
        state.tool = state.tool === v ? null : v;
        state.pending = null;
        render();
      });
      var more = toolBox.appendChild(el("button", { type: "button" }, "+ qubit"));
      more.addEventListener("click", function () {
        if (state.n < 6) { state.n += 1; state.observable += "Z"; state.example = null; render(); }
      });
      var less = toolBox.appendChild(el("button", { type: "button" }, "− qubit"));
      less.addEventListener("click", function () {
        if (state.n > 2) {
          state.n -= 1;
          state.example = null;
          state.observable = state.observable.slice(0, state.n);
          state.ops = state.ops.filter(function (op) { return op.qubits.every(function (q) { return q < state.n; }); });
          render();
        }
      });
    }

    function tapGate(op, ops) {
      state.auto = false;
      var i = -1;
      state.cuts.forEach(function (c, k) { if (c.opId === op.id) i = k; });
      var name = nameOf(op, ops);
      if (i < 0) {
        if (state.cuts.length >= MAX_CUTS) { state.said = "At most " + MAX_CUTS + " cuts here. Tap a pair of scissors to remove one."; return; }
        state.cuts.push({ kind: "gate", opId: op.id });
        state.said = name + " has a gate cut. Tap it again for a Hadamard cut.";
      } else if (state.cuts[i].kind === "gate") {
        state.cuts[i].kind = "hadamard";
        state.said = name + " has a Hadamard cut. Tap it again to join it back.";
      } else {
        state.cuts.splice(i, 1);
        state.said = name + " is whole again.";
      }
    }
    function tapWire(edge, ops) {
      state.auto = false;
      var fromId = ops[edge.from].id, toId = ops[edge.to].id, i = -1;
      state.cuts.forEach(function (c, k) { if (c.kind === "wire" && c.fromId === fromId && c.toId === toId && c.qubit === edge.qubit) i = k; });
      if (i >= 0) {
        state.cuts.splice(i, 1);
        state.said = "The line of q" + edge.qubit + " is whole again.";
      } else if (state.cuts.length < MAX_CUTS) {
        state.cuts.push({ kind: "wire", qubit: edge.qubit, fromId: fromId, toId: toId });
        state.said = "A wire cut on q" + edge.qubit + ". Tap the scissors to join it back.";
      } else {
        state.said = "At most " + MAX_CUTS + " cuts here. Tap a pair of scissors to remove one.";
      }
    }
    function tapCell(q, col) {
      var tool = state.tool, here = occupied(q, col);
      state.example = null;
      if (tool === "erase") {
        if (here) state.ops.splice(state.ops.indexOf(here), 1);
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

    function buildCircuit(ops) {
      return {
        n: state.n,
        ops: ops.map(function (op) {
          var o = { name: op.name, qubits: op.qubits, id: op.id };
          if (op.theta !== undefined) o.theta = op.theta;
          return o;
        })
      };
    }
    function cutsFor(ops) {
      return state.cuts.map(function (c) {
        return c.opId !== undefined
          ? { kind: c.kind, op: indexOf(ops, c.opId) }
          : { kind: "wire", qubit: c.qubit, from: indexOf(ops, c.fromId), to: indexOf(ops, c.toId) };
      });
    }

    /* Colours and new scissors ease in, keeping the page light on a tablet. */
    function paint() {
      var seen = {};
      Array.prototype.forEach.call(boardHost.querySelectorAll("[data-k]"), function (node) {
        var key = node.getAttribute("data-k"), value = node.getAttribute("fill") || node.getAttribute("stroke");
        seen[key] = value;
        if (key.charAt(0) === "x") {
          if (!(key in painted)) node.classList.add("qdx-pop");
        } else if (painted[key] && painted[key] !== value) {
          var property = node.getAttribute("fill") ? "fill" : "stroke";
          node.style[property] = painted[key];
          node.getBoundingClientRect();
          node.style[property] = "";
        }
      });
      painted = seen;
    }

    function render() {
      editBody.hidden = !state.editing;
      editButton.setAttribute("aria-pressed", String(state.editing));
      if (state.editing) drawTools();
      drawExamples();
      auto.setAttribute("aria-pressed", String(state.auto));
      var ops = ordered();
      dropStaleCuts(ops);
      var circuit = buildCircuit(ops), cols = ops.map(function (op) { return op.col; }), cuts = cutsFor(ops);
      var plan = E.plan(circuit, cuts, state.observable), edges = wireEdges(ops), editing = !!state.tool;
      var cutOps = {};
      cuts.forEach(function (c) { if (c.kind !== "wire") cutOps[c.op] = true; });

      var svg = drawCircuit({
        n: state.n, ops: circuit.ops, cols: cols, cuts: cuts, plan: plan, keys: true,
        width: state.editing ? COLS : Math.max(5, Math.max.apply(null, cols.concat([0])) + 1),
        ends: state.observable.split(""),
        aria: "The circuit. Each colour is a fragment. Tap a two-qubit gate or a line to cut there.",
        overlay: function (g) {
          var s = "";
          if (state.pending) s += '<circle class="pending" cx="' + g.x(state.pending.col) + '" cy="' + g.y(state.pending.q) + '" r="16"/>';
          for (var q = 0; q < state.n; q++) {
            s += '<rect class="qdx-hit" role="button" tabindex="0" aria-label="Qubit ' + q + ' measures ' + state.observable.charAt(q) + ': tap to change it" data-obs="' + q + '" x="' + (g.W - 42) + '" y="' + (g.y(q) - ROW / 2) + '" width="42" height="' + ROW + '" rx="8"/>';
          }
          if (editing) {
            for (q = 0; q < state.n; q++) {
              for (var c = 0; c < COLS; c++) {
                if (!occupied(q, c)) s += '<circle class="qdx-slot" cx="' + g.x(c) + '" cy="' + g.y(q) + '" r="3"/>';
                s += '<rect class="qdx-hit" data-q="' + q + '" data-c="' + c + '" x="' + (g.x(c) - CELL / 2) + '" y="' + (g.y(q) - ROW / 2) + '" width="' + CELL + '" height="' + ROW + '"/>';
              }
            }
            return s;
          }
          /* Tap targets on the gates and on the lines, with a circle where a cut can go. */
          ops.forEach(function (op, i) {
            if (op.qubits.length !== 2) return;
            var a = g.y(Math.min.apply(null, op.qubits)), b = g.y(Math.max.apply(null, op.qubits)), cx = g.x(cols[i]);
            s += '<rect class="qdx-hit" role="button" tabindex="0" aria-label="Cut ' + shortName(op) + "_" + (i + 1) + '" data-op="' + i + '" x="' + (cx - 21) + '" y="' + (a - 19) + '" width="42" height="' + (b - a + 38) + '" rx="8"/>';
            if (!cutOps[i]) s += '<circle class="qdx-mark" cx="' + cx + '" cy="' + ((a + b) / 2) + '" r="5"/>';
          });
          edges.forEach(function (edge, i) {
            var cx = LEFT + ((cols[edge.from] + cols[edge.to] + 1) / 2) * CELL, y = g.y(edge.qubit);
            var on = cuts.some(function (k) { return k.kind === "wire" && k.qubit === edge.qubit && k.from === edge.from && k.to === edge.to; });
            /* A one-qubit gate does nothing when tapped, so the line may take its half too. */
            var reach = function (op) { return op.qubits.length === 1 ? 0 : 13; };
            var left = Math.min(g.x(cols[edge.from]) + reach(ops[edge.from]), cx - 13);
            var right = Math.max(g.x(cols[edge.to]) - reach(ops[edge.to]), cx + 13);
            s += '<rect class="qdx-hit" role="button" tabindex="0" aria-label="Cut the line of qubit ' + edge.qubit + " between gates " + (edge.from + 1) + " and " + (edge.to + 1) + '" data-edge="' + i + '" x="' + left + '" y="' + (y - 17) + '" width="' + (right - left) + '" height="34" rx="6"/>';
            if (!on) s += '<circle class="qdx-mark qdx-wiremark" cx="' + cx + '" cy="' + y + '" r="5"/>';
          });
          return s;
        }
      });
      board(boardHost, svg, true);
      /* The box follows the circuit: narrower or wider, shorter or taller, as it changes. */
      var drawn = LEFT + (state.editing ? COLS : Math.max(5, Math.max.apply(null, cols.concat([0])) + 1)) * CELL + 50;
      boardHost.style.width = "min(100%, " + Math.round(drawn * 1.35) + "px)";
      if (keyed) {  /* with a keyboard, stay on the same place */
        Array.prototype.some.call(boardHost.querySelectorAll(".qdx-hit"), function (node) {
          if (node.getAttribute("data-op") === keyed.op && node.getAttribute("data-edge") === keyed.edge && node.getAttribute("data-obs") === keyed.obs) { node.focus(); return true; }
        });
        keyed = null;
      }
      paint();

      /* What just happened, or what to do */
      hint.textContent = editing
        ? (state.tool === "erase" ? "Tap a gate to remove it." : TWO[state.tool] ? "Tap one qubit, then another in the same column." : "Tap an empty place on the circuit.")
        : state.said || (state.editing ? "Pick a gate here, then tap where it goes on the lines. Tap Edit the circuit again when you are done." : ops.length ? "Tap a two-qubit gate to cut it, or a line to cut the wire. Or try Auto cut." : "An empty circuit.");
      angleRow.style.display = ANGLED[state.tool] ? "" : "none";

      /* The result */
      var observable = state.observable, exact = E.exact(circuit, observable), bad = plan.invalid.length > 0;
      warn.hidden = !bad;
      payoff.textContent = "";
      strip.textContent = "";
      shotStats.textContent = "";
      termsBody.textContent = "";
      if (bad) {
        warn.textContent = "This cut leaves its two sides joined by another path, so it separates nothing. Qdislib refuses it with InvalidCutError. Cut that path as well, or tap the scissors to remove it.";
        payoff.appendChild(payoffValues(exact, null, observable));
      } else {
        var prepared = E.prepare(circuit, cuts, observable);
        payoff.appendChild(payoffValues(exact, cuts.length ? prepared.value : null, observable));
        payoff.appendChild(payoffCost(plan, cuts, prepared));
        drawShots(prepared, plan, exact);
        drawTerms(plan, cuts, circuit.ops);
      }
      clear.disabled = !state.cuts.length;
      codeBody.textContent = pythonFor(ops, cuts, observable);
    }

    function payoffValues(exact, value, observable) {
      var box = el("div", { class: "qdx-values" });
      var W = 340, x = function (v) { return 14 + ((v + 1) / 2) * (W - 28); };
      var numbers = box.appendChild(el("div", { class: "qdx-numbers" }));
      var first = numbers.appendChild(el("div"));
      first.appendChild(el("span", { class: "qdx-key qdx-ring" }, "Uncut circuit"));
      first.appendChild(el("b", {}, num(exact, 4)));
      var second = numbers.appendChild(el("div"));
      second.appendChild(el("span", { class: "qdx-key qdx-dot" }, "Reconstructed"));
      second.appendChild(el("b", {}, value === null ? "—" : num(value, 4)));
      var line = '<svg viewBox="0 0 ' + W + ' 44" style="max-width:560px" role="img" aria-label="Both values on a line from minus one to plus one">' +
        '<line x1="14" x2="' + (W - 14) + '" y1="16" y2="16" stroke="var(--qd-line)" stroke-width="4" stroke-linecap="round"/>' +
        '<text x="14" y="38" text-anchor="middle">−1</text><text x="' + x(0) + '" y="38" text-anchor="middle">0</text><text x="' + (W - 14) + '" y="38" text-anchor="middle">+1</text>' +
        '<circle cx="' + x(exact) + '" cy="16" r="10" fill="none" stroke="var(--qd-b)" stroke-width="3"/>';
      if (value !== null) line += '<circle cx="' + x(value) + '" cy="16" r="5.5" fill="var(--qd-a)"/>';
      box.appendChild(el("div", { class: "qdx-line" })).innerHTML = line + "</svg>";
      box.appendChild(el("p", { class: "qdx-verdict" }, value === null
        ? "<" + observable + "> of the whole circuit. No cuts yet."
        : Math.abs(value - exact) < 1e-9
          ? "The same value, from the pieces alone: <" + observable + ">, difference " + small(Math.abs(value - exact)) + "."
          : "These should agree: difference " + Math.abs(value - exact).toExponential(1)));
      return box;
    }
    function payoffCost(plan, cuts, prepared) {
      var box = el("div", { class: "qdx-cost" });
      var pieces = box.appendChild(el("div", { class: "qdx-pieces" }));
      plan.fragments.forEach(function (fragment, f) {
        var piece = pieces.appendChild(el("span", { class: "qdx-piece", title: "Fragment " + (f + 1) }));
        for (var i = 0; i < plan.widths[f]; i++) piece.appendChild(el("i", { style: "background:" + fragColor(f) }));
      });
      var largest = Math.max.apply(null, plan.widths.concat([0]));
      box.appendChild(el("p", { class: "qd-note" }, !cuts.length
        ? "The whole circuit needs " + state.n + " qubits at once" + (plan.fragments.length > 1 ? ", in pieces that already do not touch." : ".")
        : "Largest piece: " + largest + " qubit" + (largest === 1 ? "" : "s") + ", from a circuit of " + state.n + "."));
      var stats = box.appendChild(el("div", { class: "qd-stats" }));
      stat(stats, count(plan.subcircuits), plan.subcircuits === 1 ? "circuit to run" : "subcircuits to run");
      if (plan.optimized !== null && cuts.length) stat(stats, count(plan.optimized), "with optimized=True");
      stat(stats, gammaText(plan.gamma), "gamma");
      if (cuts.length) stat(stats, "×" + count(plan.gamma * plan.gamma), "shots, in the worst case");
      var notes = [], goldenCount = Object.keys(plan.golden).length;
      if (goldenCount) notes.push((goldenCount > 1 ? goldenCount + " wire cuts are golden" : "This wire cut is golden") + ": the state before it is real, so its Y terms are zero and are left out, with 6 terms and gamma 3.");
      if (plan.joint && plan.joint.shared) notes.push("Qdislib's default, joint=True, shares ancillas between the " + plan.joint.shared + " Hadamard cuts across one boundary: " + count(plan.joint.subcircuits) + " subcircuits, gamma " + gammaText(plan.joint.gamma) + ". The value is the same.");
      if (notes.length) box.appendChild(el("p", { class: "qd-note" }, notes.join(" ")));
      return box;
    }
    /* The axis stays put along the slider, set by the fewest shots, so the spread is seen to
       shrink. It belongs to the circuit and its cuts, not to one draw. */
    var axis = { key: null, extent: 1 };
    function axisFor(prepared, exact, key) {
      if (axis.key === key) return axis.extent;
      var rng = E.mulberry(7919), most = 1e-3, fewest = Math.pow(2, Number(shotSlider.min));
      for (var i = 0; i < 60; i++) most = Math.max(most, Math.abs(prepared.draw(fewest, rng) - exact));
      axis = { key: key, extent: most };
      return most;
    }
    function drawShots(prepared, plan, exact) {
      var shots = Math.pow(2, Number(shotSlider.value)), rng = E.mulberry(state.seed * 7919 + shots), draws = [];
      for (var i = 0; i < 30; i++) draws.push(prepared.draw(shots, rng));
      var extent = axisFor(prepared, exact, JSON.stringify([state.ops, state.cuts, state.observable]));
      var mean = draws.reduce(function (a, b) { return a + b; }, 0) / draws.length;
      var spread = Math.sqrt(draws.reduce(function (a, b) { return a + (b - mean) * (b - mean); }, 0) / (draws.length - 1));
      stat(shotStats, count(shots), "shots per subcircuit");
      stat(shotStats, num(draws[0], 4), "one run, exact " + num(exact, 4));
      stat(shotStats, spread.toFixed(4), "spread of 30 runs");
      stat(shotStats, count(shots * plan.subcircuits), "shots in total, per run");
      strip.innerHTML = stripPlot(draws, exact, extent);
    }
    function drawTerms(plan, cuts, ops) {
      if (!cuts.length) { termsBody.appendChild(el("p", { class: "qd-note" }, "No cuts yet.")); return; }
      cuts.forEach(function (cut, index) {
        var title = cut.kind === "wire"
          ? "Wire cut on q" + cut.qubit + (plan.golden[index] ? " (golden)" : "")
          : (cut.kind === "gate" ? "Gate cut of " : "Hadamard cut of ") + shortName(ops[cut.op]) + "_" + (cut.op + 1);
        termsBody.appendChild(el("p", {}, title));
        var isHadamard = cut.kind === "hadamard", list = E.terms(plan, cut, index).slice();
        /* In the order of Qdislib's describe(), whose (i,j) these are. */
        if (isHadamard) list.sort(function (x, y) {
          var p = x.a.had, q = y.a.had;
          return p.i - q.i || p.j - q.j || p.g - q.g;
        });
        var rows = list.map(function (t) {
          var a, b, term = "";
          if (cut.kind === "wire") { a = t.read === "I" ? "read nothing" : "measure " + t.read; b = "prepare |" + t.prep + "\u27e9"; }
          else if (cut.kind === "gate") { a = termName(t.a); b = termName(t.b); }
          else {
            a = hadName(t.a.had, 0);
            b = hadName(t.b.had, 1);
            var h = t.a.had;
            term = '<td class="qdx-term">(' + h.i + "," + h.j + ")" + (h.i === h.j ? "" : "<br>" + (h.g ? "imaginary" : "real")) + "</td>";
          }
          return "<tr>" + term + "<td>" + (t.w >= 0 ? "+" : "\u2212") + Math.abs(t.w).toFixed(4) + "</td><td>" + esc(a) + "</td><td>" + esc(b) + "</td></tr>";
        }).join("");
        termsBody.appendChild(el("table", { class: "qdx-t" })).innerHTML = "<thead><tr>" + (isHadamard ? "<th>Term</th>" : "") + "<th>Weight</th><th>Side 1</th><th>Side 2</th></tr></thead><tbody>" + rows + "</tbody>";
        if (isHadamard) termsBody.appendChild(el("p", { class: "qd-note" }, "(i,j) as subs.describe() prints it, with parity 0 for the real part and 1 for the imaginary one. |0\u27e9 and |1\u27e9 are the ancilla's two branches. Swapping the operators conjugates the term: the real parts match and the imaginary parts flip sign, which cancels in the product, so each pair adds the same amount."));
      });
    }
    function termName(a) {
      return { i: "nothing", z: "Z", m: "measure, sign", "r+": "RZ(π/2)", "r-": "RZ(−π/2)" }[a];
    }
    /* What the engine applies: U[i] on the ancilla's |0>, U[j] on its |1>. Side 1 carries the
       phase of the gate's -iZ. With i = j there is no ancilla, and that phase changes nothing. */
    function hadName(h, side) {
      if (h.i === h.j) return h.i ? "Z" : "nothing";
      var U = ["I", side ? "Z" : (h.phase < 0 ? "\u2212iZ" : "iZ")];
      return U[h.i] + " on |0\u27e9, " + U[h.j] + " on |1\u27e9";
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
    function stripPlot(draws, exact, extent) {
      var W = 340, H = 64, L = 8, R = 8, span = extent * 1.15;
      var x = function (v) { return L + ((Math.max(-span, Math.min(span, v - exact)) + span) / (2 * span)) * (W - L - R); };
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" style="max-width:560px" role="img" aria-label="30 runs against the exact value, on an axis that stays the same at every number of shots">';
      s += '<line x1="' + x(exact) + '" x2="' + x(exact) + '" y1="6" y2="42" stroke="var(--qd-b)" stroke-width="2"/>';
      draws.forEach(function (d, i) { s += '<circle cx="' + x(d) + '" cy="' + (12 + (i % 5) * 6) + '" r="3.2" fill="var(--qd-a)" opacity="0.8"/>'; });
      s += '<text x="' + x(exact) + '" y="56" text-anchor="middle">exact ' + num(exact, 3) + "</text>";
      s += '<text x="' + L + '" y="56">' + num(exact - span, 3) + "</text>";
      s += '<text x="' + (W - R) + '" y="56" text-anchor="end">' + num(exact + span, 3) + "</text>";
      return s + "</svg>";
    }

    /* Tapping */
    function activate(t) {
      if (t.dataset.obs !== undefined) {
        var q = Number(t.dataset.obs), order = "ZXYI", now = state.observable.charAt(q);
        state.observable = state.observable.slice(0, q) + order.charAt((order.indexOf(now) + 1) % 4) + state.observable.slice(q + 1);
        state.said = order.charAt((order.indexOf(now) + 1) % 4) === "I" ? "Nothing measured on q" + q + ". Tap the letter again for Z." : "Measuring " + order.charAt((order.indexOf(now) + 1) % 4) + " on q" + q + ". Tap the letter again to change it.";
        render();
      } else if (t.dataset.edge !== undefined) {
        var ops = ordered();
        tapWire(wireEdges(ops)[Number(t.dataset.edge)], ops);
        render();
      } else if (t.dataset.op !== undefined) {
        var list = ordered();
        tapGate(list[Number(t.dataset.op)], list);
        render();
      } else if (t.dataset.q !== undefined) {
        tapCell(Number(t.dataset.q), Number(t.dataset.c));
      }
    }
    boardHost.addEventListener("click", function (event) { activate(event.target); });
    boardHost.addEventListener("keydown", function (event) {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        keyed = { op: event.target.getAttribute("data-op"), edge: event.target.getAttribute("data-edge"), obs: event.target.getAttribute("data-obs") };
        activate(event.target);
      }
    });
    auto.addEventListener("click", function () {
      var ops = ordered(), circuit = buildCircuit(ops), found = E.autoCut(circuit, state.observable);
      state.cuts = found.map(function (c) { return { kind: "gate", opId: ops[c.op].id }; });
      state.auto = found.length > 0;
      if (found.length) {
        state.said = "Auto cut: every split of the qubits in two, keeping the lowest gamma. A simple rule; qd.find_cut searches further.";
      } else {
        state.said = !ops.length ? "Nothing to cut yet: this circuit is empty." : E.plan(circuit, [], state.observable).fragments.length > 1
          ? "Nothing to cut: this circuit is already in pieces."
          : "Auto cut found nothing: splitting these qubits in two takes more than " + MAX_CUTS + " gate cuts.";
      }
      render();
    });
    clear.addEventListener("click", function () { state.cuts = []; state.said = ""; state.auto = false; render(); });
    shotSlider.addEventListener("input", render);
    again.addEventListener("click", function () { state.seed += 1; render(); });
    EXAMPLES.chain[1]();
    render();
  }

  /* ---- What lowers the cost ---- */
  function lowerSection(root) {
    var data = D.lower, k = 3, angle = 1;
    root.textContent = "";
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-low-k" }));
    var slider = row.appendChild(el("input", { id: "qdx-low-k", type: "range", min: "1", max: String(data.cuts.length), step: "1", value: String(k), "aria-label": "Number of cuts" }));
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
    var data = D.find_cut, family = data[0].key, kind = "gate", max = 4;
    root.textContent = "";
    picker(root, "qdx-find-circuit", "Circuit", data.map(function (x) { return { value: x.key, label: x.label }; }), family, function (v) { family = v; update(); });
    var kindBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-max" }));
    var slider = row.appendChild(el("input", { id: "qdx-max", type: "range", step: "1", "aria-label": "max_qubits" }));
    var host = root.appendChild(el("div", { id: "qdx-find-board" }));
    root.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendHtml();
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var code = codeCell(root);
    var note = root.appendChild(el("p", { class: "qd-note" }));

    function update() {
      var f = data.filter(function (x) { return x.key === family; })[0];
      toggles(kindBox, [{ value: "gate", label: "Gate" }, { value: "wire", label: "Wire" }, { value: "mixed", label: "Any, by cost" }], kind, function (v) { kind = v; update(); });
      slider.min = "2"; slider.max = String(f.n - 1);
      max = Math.min(max, f.n - 1);
      slider.value = String(max);
      label.textContent = "max_qubits = " + max;
      var run = f.runs[kind + ":" + max], cuts = run.cuts || [];
      var plan = E.plan({ n: f.n, ops: f.ops }, cuts.map(function (c) { return Object.assign({}, c); }), new Array(f.n + 1).join("Z"));
      var cols = columns(f.n, f.ops), columnsUsed = Math.max.apply(null, cols) + 1;
      board(host, drawCircuit({ n: f.n, ops: f.ops, cols: cols, cuts: cuts, plan: plan, width: columnsUsed, aria: f.label }), true);
      /* The box follows the circuit, as in the playground. */
      host.style.width = "min(100%, " + (LEFT + columnsUsed * CELL + 12) + "px)";
      stats.textContent = "";
      if (cuts.length) {
        stat(stats, count(run.subcircuits), "subcircuits");
        stat(stats, gammaText(run.gamma), "gamma");
        stat(stats, String(Math.max.apply(null, plan.widths)), "qubits in the largest piece");
      } else {
        stat(stats, "none", "find_cut found no cut within this limit");
      }
      var options = kind === "mixed" ? "mixed=True" : kind === "gate" ? "wire_cut=False" : "gate_cut=False";
      code.textContent = "qd.find_cut(circuit, max_qubits=" + max + ", " + options + ", seed=0)\n# " + (cuts.length ? JSON.stringify(run.spec) : "[]");
      note.textContent = (run.difference !== undefined ? "Reconstructed with shots=None, this cut differs from the uncut value by " + small(run.difference) + ". " : "") + "Stored from a real run; the largest piece is computed here from the cut.";
    }
    slider.addEventListener("input", function () { max = Number(slider.value); update(); });
    update();
  }

  /* ---- 3. Planner ---- */
  function plannerSection(root) {
    var data = D.planner, k = 3;
    root.textContent = "";
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-k" }));
    var slider = row.appendChild(el("input", { id: "qdx-k", type: "range", min: "1", max: String(data.cuts.length), step: "1", value: String(k), "aria-label": "Number of cuts" }));
    var table = root.appendChild(el("table", { class: "qdx-t" }));
    root.appendChild(el("p", { class: "qd-note" },
      "An almost exact result: within ±" + data.epsilon + " of the true value, with probability " + (1 - data.delta) + ". The shots are what qd.shots_for(circuit, cuts, epsilon=" + data.epsilon + ") returns, 2 (gamma / epsilon)² ln(2 / delta), about " + (2 * Math.log(2 / data.delta)).toFixed(1) + " (gamma / epsilon)²: so gamma drives them. It holds for any circuit with this cut, so a real run can need fewer. Stored from qd.cutting_overhead and qd.shots_for."));
    function update() {
      label.textContent = "Cuts: " + k;
      var html = "<thead><tr><th>Cut</th><th>Subcircuits</th><th>gamma</th><th>Shots for ±" + data.epsilon + "</th></tr></thead><tbody>";
      Object.keys(data.kinds).forEach(function (key) {
        var kind = data.kinds[key], r = kind.rows[k - 1];
        html += "<tr><td>" + esc(kind.label) + "</td><td>" + count(r.subcircuits) + "</td><td>" + gammaText(r.gamma) + "</td><td>" + count(r.shots) + "</td></tr>";
      });
      table.innerHTML = html + "</tbody>";
    }
    slider.addEventListener("input", function () { k = Number(slider.value); update(); });
    update();
  }

  /* ---- 4. Distributions ---- */
  /* Bars for one run, ticks for the exact values. A reconstruction can dip below zero. */
  function hist(keys, exact, values, o) {
    var W = o.width || 170, H = o.height || 118, T = 20, B = keys.length <= 4 ? 20 : 4, L = 4, R = 4;
    var top = o.top, low = o.low, y = function (v) { return T + ((top - v) / (top - low)) * (H - T - B); };
    var step = (W - L - R) / keys.length, bar = Math.min(Math.max(step - 2, 2), 34), s = '<svg viewBox="0 0 ' + W + " " + H + '" role="img" aria-label="' + esc(o.title || "Distribution") + '">';
    if (o.title) s += '<text class="t" x="' + L + '" y="12">' + esc(o.title) + "</text>";
    keys.forEach(function (key, i) {
      var v = values[i], x = L + i * step + (step - bar) / 2;
      var note = key + ": " + (exact ? "exact " + exact[i].toFixed(4) + ", " : "") + "this run " + v.toFixed(4);
      s += "<g><title>" + esc(note) + "</title>";
      s += '<rect x="' + (L + i * step) + '" y="' + T + '" width="' + step + '" height="' + (H - T - B) + '" fill="transparent"/>';
      s += '<rect x="' + x + '" y="' + Math.min(y(0), y(v)) + '" width="' + bar + '" height="' + Math.max(1, Math.abs(y(v) - y(0))) + '" rx="1.5" fill="var(--qd-a)"/>';
      if (exact) s += '<line x1="' + (x - 0.5) + '" x2="' + (x + bar + 0.5) + '" y1="' + y(exact[i]) + '" y2="' + y(exact[i]) + '" stroke="var(--qd-b)" stroke-width="2.5"/>';
      s += "</g>";
      if (keys.length <= 4) {
        s += '<text x="' + (x + bar / 2) + '" y="' + (H - 5) + '" text-anchor="middle">' + esc(key) + "</text>";
      }
    });
    s += '<line x1="' + L + '" x2="' + (W - R) + '" y1="' + y(0) + '" y2="' + y(0) + '" stroke="var(--color-foreground-secondary, #555)" stroke-width="1"/>';
    return s + "</svg>";
  }
  /* About three plots at different shots, on one scale, a little larger than the full ones. */
  function shotPlots(host, shots, runs, exact) {
    var keys = ["00", "01", "10", "11"], all = [0].concat(exact || []), width = window.matchMedia("(max-width: 600px)").matches ? 120 : 170;
    runs.forEach(function (run) { all = all.concat(run); });
    var top = Math.max.apply(null, all) * 1.05, low = Math.min.apply(null, all);
    var grid = host.appendChild(el("div", { class: "qdx-grid qdx-grid3" }));
    shots.forEach(function (count_, i) {
      grid.appendChild(el("div")).innerHTML = hist(keys, exact, runs[i], { width: width, top: top, low: low, title: count(count_) + " shots" });
    });
  }
  function legendOf(items) {
    return items.map(function (x) {
      var mark = x[2] === "tick"
        ? '<span class="qd-swatch" style="height:3px;vertical-align:4px;background:' + x[0] + '"></span>'
        : '<span class="qd-swatch" style="background:' + x[0] + '"></span>';
      return "<span>" + mark + esc(x[1]) + "</span>";
    }).join("");
  }
  function distributionSection(root) {
    var data = D.distributions, tab = "full", seed = 0, a = 0, b = 7, pair = 10;
    root.textContent = "";
    var tabs = root.appendChild(el("div", { class: "qdx-tabs" }));
    var body = root.appendChild(el("div"));
    function update() {
      toggles(tabs, [{ value: "full", label: "Full distribution" }, { value: "marginal", label: "qubits= marginal" }, { value: "big", label: "40 qubits" }], tab, function (v) { tab = v; update(); });
      body.textContent = "";
      if (tab === "full") {
        var f = data.full, top = 0, low = 0;
        f.exact.forEach(function (v) { top = Math.max(top, v); });
        f.runs.forEach(function (level) { level.forEach(function (run) { run.top.forEach(function (v) { top = Math.max(top, v); low = Math.min(low, v); }); }); });
        body.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendOf([["var(--qd-a)", "circuit_cutting_sampling, one run"], ["var(--qd-b)", "exact", "tick"]]);
        var grid = body.appendChild(el("div", { class: "qdx-grid" }));
        f.shots.forEach(function (shots, level) {
          var run = f.runs[level][seed];
          grid.appendChild(el("div")).innerHTML = hist(f.keys, f.exact, run.top, { top: top * 1.05, low: low, title: count(shots) + " shots" }) +
            '<p class="qd-note">distance from exact ' + run.tvd.toFixed(3) + "</p>";
        });
        var again = body.appendChild(el("button", { type: "button" }, "Another run (" + (seed + 1) + " of 3)"));
        again.addEventListener("click", function () { seed = (seed + 1) % 3; update(); });
        body.appendChild(el("p", { class: "qd-note" }, "Shots per subcircuit. Each bar is one of the 12 likeliest of 256 outcomes, from " + f.keys[0] + " on, 8 qubits with one gate cut, " + data.cut + ". All plots share one scale; the distance is total variation over all 256. Stored from real runs."));
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
        var code = codeCell(body);
        code.textContent = 'qd.circuit_cutting_sampling(circuit, ["' + data.cut + '"], shots=4096, seed=1, qubits=[' + lo + ", " + hi + "])";
        body.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendOf([["var(--qd-a)", "qubits=[" + lo + ", " + hi + "]"], ["var(--qd-b)", "exact", "tick"]]);
        shotPlots(body, data.small_shots, m.runs, m.exact);
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
        body.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendOf([["var(--qd-a)", "one run"], ["var(--qd-b)", "exact, with shots=None", "tick"]]);
        shotPlots(body, data.small_shots, data.small_shots.map(function (n) { return big.sampled[n][pair + "," + (pair + 1)]; }), values);
        sl.addEventListener("input", function () { pair = Number(sl.value); update(); });
        body.appendChild(el("span", { class: "qd-label" }, "d.most_likely(5)"));
        var pre = codeCell(body);
        pre.textContent = Object.keys(big.most_likely).map(function (k) { return k + "  " + big.most_likely[k].toExponential(2); }).join("\n");
        body.appendChild(el("p", { class: "qd-note" }, big.text + ". Exact, in " + big.seconds + " s on one machine. Stored from a real run; no circuit runs when it is read. On qubits 2–7, its marginals differ from a 10-qubit statevector of the same gates by " + small(big.check) + "."));
      }
    }
    update();
  }

  /* ---- 5. Error bars ---- */
  function whisker(x, y1, y2, colour) {
    return '<line x1="' + x + '" x2="' + x + '" y1="' + y1 + '" y2="' + y2 + '" stroke="' + colour + '" stroke-width="1.5"/>' +
      '<line x1="' + (x - 3) + '" x2="' + (x + 3) + '" y1="' + y1 + '" y2="' + y1 + '" stroke="' + colour + '" stroke-width="1.5"/>' +
      '<line x1="' + (x - 3) + '" x2="' + (x + 3) + '" y1="' + y2 + '" y2="' + y2 + '" stroke="' + colour + '" stroke-width="1.5"/>';
  }
  /* A few round values on the y axis, faint. */
  function yGrid(low, high, y, x1, x2) {
    var span = high - low, unit = Math.pow(10, Math.floor(Math.log10(span)) - 1), step = [1, 2, 5, 10].map(function (m) { return m * unit; }).filter(function (d) { return span / d <= 5; })[0];
    var s = "", digits = step < 0.1 ? 2 : 1;
    if (step < 0.01) digits = 3;
    for (var v = Math.ceil(low / step) * step; v <= high + 1e-12; v += step) {
      s += '<line x1="' + x1 + '" x2="' + x2 + '" y1="' + y(v) + '" y2="' + y(v) + '" stroke="var(--qd-line)" stroke-width="1"/>';
      s += '<text x="' + (x1 - 4) + '" y="' + (y(v) + 4) + '" text-anchor="end">' + v.toFixed(digits) + "</text>";
    }
    return s;
  }
  function errorSection(root) {
    var data = D.errors, key = "gate";
    root.textContent = "";
    var tabs = root.appendChild(el("div", { class: "qdx-tabs" }));
    root.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendOf([["var(--qd-a)", "one run, with its error=True bar"], ["var(--qd-b)", "24 runs: each a dot, their mean and spread"]]);
    var plot = root.appendChild(el("div", { class: "qdx-plot qdx-plot-small" }));
    var table = root.appendChild(el("div"));
    root.appendChild(el("p", { class: "qd-note" }, "The purple bar is the ±1 standard error one run reported, from a bootstrap over its shots, at no device cost. The orange bar is how far 24 separate runs really scatter, each with its own seed: when the two match, the reported bar can be trusted. Stored from real runs."));
    function update() {
      toggles(tabs, Object.keys(data.cases).map(function (k) { return { value: k, label: data.cases[k].label }; }), key, function (v) { key = v; update(); });
      var rows = data.cases[key].rows, W = 360, H = 230, L = 44, R = 8, T = 22, B = 40;
      var points = [data.exact], means = rows.map(function (r) {
        return r.values.reduce(function (a, v) { return a + v; }, 0) / r.values.length;
      });
      rows.forEach(function (r, i) {
        points = points.concat(r.values, [r.value - r.reported, r.value + r.reported, means[i] - r.spread, means[i] + r.spread]);
      });
      var low = Math.min.apply(null, points), high = Math.max.apply(null, points), pad = (high - low) * 0.05;
      low -= pad; high += pad;
      var y = function (v) { return T + ((high - v) / (high - low)) * (H - T - B); };
      var band = (W - L - R) / rows.length;
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" role="img" aria-label="Reported error bars next to the spread of 24 runs, at each number of shots">';
      s += yGrid(low, high, y, L, W - R);
      s += '<line x1="' + L + '" x2="' + (W - R) + '" y1="' + y(data.exact) + '" y2="' + y(data.exact) + '" stroke="var(--color-foreground-primary, #111)" stroke-width="1.5" stroke-dasharray="5 3"/>';
      s += '<text class="t" x="' + (W - R) + '" y="14" text-anchor="end">dashed line: exact ' + num(data.exact, 3) + "</text>";
      rows.forEach(function (r, i) {
        var c = L + band * (i + 0.5), xa = c - 15, xb = c + 7;
        r.values.forEach(function (v, j) {
          s += '<circle cx="' + (c + 19 + (j % 3) * 4) + '" cy="' + y(v) + '" r="1.8" fill="var(--qd-b)" opacity="0.45"/>';
        });
        s += "<g><title>" + count(data.shots[i]) + " shots: one run " + num(r.value, 4) + " ±" + r.reported.toFixed(4) + "</title>" + whisker(xa, y(r.value + r.reported), y(r.value - r.reported), "var(--qd-a)") + '<circle cx="' + xa + '" cy="' + y(r.value) + '" r="4" fill="var(--qd-a)"/></g>';
        s += "<g><title>" + count(data.shots[i]) + " shots: 24 runs, mean " + num(means[i], 4) + ", spread " + r.spread.toFixed(4) + "</title>" + whisker(xb, y(means[i] + r.spread), y(means[i] - r.spread), "var(--qd-b)") + '<circle cx="' + xb + '" cy="' + y(means[i]) + '" r="4" fill="var(--qd-b)"/></g>';
        s += '<text x="' + c + '" y="' + (H - B + 16) + '" text-anchor="middle">' + count(data.shots[i]) + "</text>";
      });
      s += '<text x="' + (L + (W - L - R) / 2) + '" y="' + (H - 6) + '" text-anchor="middle">shots per subcircuit</text>';
      plot.innerHTML = s + "</svg>";
      var html = '<table class="qdx-t"><thead><tr><th>Shots</th><th>error=True, one run</th><th>Spread of 24</th><th>Shots in one run</th></tr></thead><tbody>';
      rows.forEach(function (r, i) {
        html += "<tr><td>" + count(data.shots[i]) + "</td><td>±" + r.reported.toFixed(4) + "</td><td>" + r.spread.toFixed(4) + "</td><td>" + count(r.total) + "</td></tr>";
      });
      table.innerHTML = html + "</tbody></table>";
    }
    update();
  }

  /* ---- 6. Noise ---- */
  function noiseSection(root) {
    var data = D.noise, W = 380, H = 240, L = 44, R = 106, T = 14, B = 38;
    var SERIES = [
      { name: "Uncut", colour: "var(--qdx-plain)", dx: 0, get: function (r) { return r.uncut; } },
      { name: "Gate cut", colour: "var(--qdx-gate)", dx: -3, get: function (r) { return r.gate; } },
      { name: "Hadamard cut", colour: "var(--qdx-had)", dx: 3, get: function (r) { return r.hadamard; } },
      { name: "Wire cut", colour: "var(--qdx-wire)", dx: 6, get: function (r) { return r.wire; } }
    ];
    root.textContent = "";
    root.appendChild(el("div", { class: "qd-legend" })).innerHTML = legendOf(SERIES.map(function (x) { return [x.colour, x.name]; }));
    var points = [data.exact];
    data.rows.forEach(function (r) {
      SERIES.forEach(function (series) {
        var v = series.get(r);
        points.push(v[0] - v[1], v[0] + v[1]);
      });
    });
    var low = Math.min(0, Math.min.apply(null, points)), high = Math.max.apply(null, points) * 1.05;
    var pMax = data.rows[data.rows.length - 1].p;
    var x = function (p) { return L + 8 + (p / pMax) * (W - L - R - 16); }, y = function (v) { return T + ((high - v) / (high - low)) * (H - T - B); };
    var s = '<svg viewBox="0 0 ' + W + " " + H + '" role="img" aria-label="The value against the noise strength, uncut and with each kind of cut">';
    s += yGrid(low, high, y, L, W - R);
    s += '<line x1="' + L + '" x2="' + (W - R) + '" y1="' + y(data.exact) + '" y2="' + y(data.exact) + '" stroke="var(--color-foreground-primary, #111)" stroke-width="1.5" stroke-dasharray="5 3"/>';
    data.rows.forEach(function (r) {
      if (r.p !== 0.005) s += '<text x="' + x(r.p) + '" y="' + (H - B + 16) + '" text-anchor="middle">' + r.p + "</text>";
    });
    s += '<text x="' + (L + (W - L - R) / 2) + '" y="' + (H - 6) + '" text-anchor="middle">noise strength p</text>';
    s += '<text transform="translate(11 ' + ((T + H - B) / 2) + ') rotate(-90)" text-anchor="middle">expectation value \u27e8ZZZZZZ\u27e9</text>';
    var ends = [{ name: "exact", value: data.exact }];
    SERIES.forEach(function (series) {
      var line = data.rows.map(function (r) { return (x(r.p) + series.dx) + "," + y(series.get(r)[0]); }).join(" ");
      s += '<polyline points="' + line + '" fill="none" stroke="' + series.colour + '" stroke-width="2"/>';
      data.rows.forEach(function (r) {
        var v = series.get(r), cx = x(r.p) + series.dx;
        s += "<g><title>" + series.name + ", p = " + r.p + ": " + num(v[0], 4) + " ±" + v[1].toFixed(4) + "</title>";
        s += whisker(cx, y(v[0] + v[1]), y(v[0] - v[1]), series.colour);
        s += '<circle cx="' + cx + '" cy="' + y(v[0]) + '" r="3" fill="' + series.colour + '"/></g>';
      });
      ends.push({ name: series.name, value: series.get(data.rows[data.rows.length - 1])[0] });
    });
    /* Labels at the right of each line, pushed apart when they would touch. */
    ends.sort(function (a, b) { return y(a.value) - y(b.value); });
    var at = -Infinity;
    ends.forEach(function (end) {
      at = Math.max(y(end.value) + 4, at + 13);
      s += '<text class="t" x="' + (W - R + 8) + '" y="' + at + '">' + end.name + "</text>";
    });
    root.appendChild(el("div", { class: "qdx-plot qdx-plot-small" })).innerHTML = s + "</svg>";
    root.appendChild(el("p", { class: "qdx-verdict" }, "As noise grows, every line falls away from the exact value. Here the gate and wire cuts stay close to the uncut circuit, within about their error bars. The Hadamard cut falls fastest: its ancillas add two-qubit gates of their own."));
    root.appendChild(el("p", { class: "qd-note" }, "<ZZZZZZ> on 6 qubits, cut once at " + data.gate + " into two pieces, exact " + num(data.exact, 4) + ", the dashed line. The model is on device gates: a depolarizing error p on each cz, p / 10 on each sx and x, none on rz, and a readout error p / 2. The backend compiles every circuit to these gates, so the noise reaches them all. All three cuts have gamma 3 here, as this wire is a golden cut point. " + count(data.shots) + " shots per circuit; bars are \u00b11 standard error. A simulation with a simple model, stored from real Aer runs: it says nothing about any device."));
    var all = root.appendChild(el("details"));
    all.appendChild(el("summary", {}, "All the values"));
    var html = '<table class="qdx-t"><thead><tr><th>p</th><th>Uncut</th><th>Gate cut</th><th>Hadamard cut</th><th>Wire cut</th></tr></thead><tbody>';
    var cell = function (v) { return num(v[0], 3) + " ±" + v[1].toFixed(3); };
    data.rows.forEach(function (r) {
      html += "<tr><td>" + r.p + "</td><td>" + cell(r.uncut) + "</td><td>" + cell(r.gate) + "</td><td>" + cell(r.hadamard) + "</td><td>" + cell(r.wire) + "</td></tr>";
    });
    all.appendChild(el("div")).innerHTML = html + "</tbody></table>";
    var code = root.appendChild(el("details"));
    code.appendChild(el("summary", {}, "The code"));
    codeCell(code, 'noise = NoiseModel(basis_gates=["cz", "sx", "x", "rz"])\nnoise.add_all_qubit_quantum_error(depolarizing_error(p, 2), ["cz"])\nnoise.add_all_qubit_quantum_error(depolarizing_error(p / 10, 1), ["sx", "x"])\nnoise.add_all_qubit_readout_error(ReadoutError([[1 - p / 2, p / 2], [p / 2, 1 - p / 2]]))\nbackend = qd.build_backend("AER_BACKEND", qd.BackendConfig(extra={"noise_model": noise}))\nqd.gate_cutting(circuit, ["' + data.gate + '"], "ZZZZZZ", backend=backend, shots=' + data.shots + ', seed=1, error=True)\nqd.gate_cutting_hadamard(circuit, ["' + data.gate + '"], "ZZZZZZ", backend=backend, shots=' + data.shots + ', seed=1, error=True)\nqd.wire_cutting(circuit, [("' + data.wire[0] + '", "' + data.wire[1] + '")], "ZZZZZZ", backend=backend, shots=' + data.shots + ", seed=1, error=True)");
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
      codeCell(body, s.install);
      codeCell(body, s.code);
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
    var all = D.pack, data = all.devices.filter(function (d) { return d.key === "montreal"; })[0], buffer = "1", job = 0;
    root.textContent = "";
    picker(root, "qdx-pack-device", "Layout", all.devices.map(function (d) {
      return { value: d.key, label: d.name + (d.family ? " (" + d.family + ")" : "") + ", " + d.qubits + " qubits" };
    }), data.key, function (v) {
      data = all.devices.filter(function (d) { return d.key === v; })[0];
      job = 0;
      update();
    });
    root.appendChild(el("span", { class: "qd-label" }, "Idle qubits between groups: buffer="));
    var bufferBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
    var row = root.appendChild(el("div", { class: "qdx-row" }));
    var label = row.appendChild(el("label", { for: "qdx-job" }));
    var slider = row.appendChild(el("input", { id: "qdx-job", type: "range", min: "0", step: "1", "aria-label": "Which job" }));
    var stats = root.appendChild(el("div", { class: "qd-stats" }));
    var map = root.appendChild(el("div", { class: "qdx-map" }));
    var chipBox = root.appendChild(el("div", { class: "qd-legend" }));
    var packNote = root.appendChild(el("p", { class: "qd-note" }));
    var COLOURS = ["#3f6fa8", "#c8640a", "#2e8b57", "#a8326a", "#7a6a1e", "#5a4fb0", "#2b8a99", "#8a4b2a", "#b0413e", "#4a7c1f", "#8e44ad", "#546e7a"];
    function update() {
      var run = data.runs[buffer];
      toggles(bufferBox, Object.keys(data.runs).map(function (b) { return { value: b, label: b }; }), buffer, function (v) { buffer = v; job = 0; update(); });
      row.style.display = run.jobs.length > 1 ? "" : "none";
      packNote.textContent = "Real device layouts, with calibration snapshots shipped in Qiskit: " + data.name + "'s is from " + data.snapshot + ". There is no live connection, and nothing ran on any device. These are the " + all.circuits + " subcircuits of a 6-qubit circuit with one cut, and each colour is one subcircuit's group. Placement only, stored from a real call to qd.pack.";
      slider.max = String(run.jobs.length - 1);
      job = Math.min(job, run.jobs.length - 1);
      slider.value = String(job);
      label.textContent = "Job " + (job + 1) + " of " + run.jobs.length;
      stats.textContent = "";
      stat(stats, String(run.jobs.length), "jobs for " + all.circuits + " subcircuits");
      stat(stats, String(run.regions.length), "groups of qubits on the device");
      stat(stats, String(run.jobs[job].length), "subcircuits in this job");
      var xs = [], ys = [];
      Object.keys(data.points).forEach(function (q) { xs.push(data.points[q][0]); ys.push(data.points[q][1]); });
      var minX = Math.min.apply(null, xs), maxX = Math.max.apply(null, xs), minY = Math.min.apply(null, ys), maxY = Math.max.apply(null, ys);
      var large = data.qubits > 40, W = large ? 640 : 340, H = large ? 580 : data.qubits < 10 ? 170 : 300, P = 16;
      var px = function (q) { var p = data.points[q]; return [P + (p[0] - minX) / ((maxX - minX) || 1) * (W - 2 * P), P + (p[1] - minY) / ((maxY - minY) || 1) * (H - 2 * P)]; };
      var owner = {};
      run.jobs[job].forEach(function (slot, i) { slot.qubits.forEach(function (q) { owner[q] = i; }); });
      var s = '<svg viewBox="0 0 ' + W + " " + H + '" style="max-width:' + (large ? 700 : 520) + 'px;margin:auto" role="img" aria-label="The device, with this job\'s groups of qubits">';
      data.edges.forEach(function (e) {
        var a = px(e[0]), b = px(e[1]);
        s += '<line x1="' + a[0] + '" y1="' + a[1] + '" x2="' + b[0] + '" y2="' + b[1] + '" stroke="var(--qd-line)" stroke-width="2"/>';
      });
      Object.keys(data.points).forEach(function (q) {
        var p = px(q), on = owner[q] !== undefined;
        s += '<circle cx="' + p[0] + '" cy="' + p[1] + '" r="9.5" fill="' + (on ? COLOURS[owner[q] % COLOURS.length] : "var(--color-background-primary, #fff)") + '" stroke="' + (on ? "none" : "var(--color-foreground-secondary, #777)") + '"/>';
        s += '<text x="' + p[0] + '" y="' + (p[1] + 3) + '" text-anchor="middle" style="fill:' + (on ? "#fff" : "var(--color-foreground-secondary, #555)") + '">' + q + "</text>";
      });
      map.innerHTML = s + "</svg>";
      chipBox.textContent = "";
      run.jobs[job].forEach(function (slot, i) {
        var item = chipBox.appendChild(el("span"));
        item.appendChild(el("span", { class: "qd-swatch", style: "background:" + COLOURS[i % COLOURS.length] }));
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
