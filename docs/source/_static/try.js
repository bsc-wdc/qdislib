/* The in-browser demo on the Try it page. It runs no quantum code: it shows runs of real
   Qdislib stored in try_data.js by docs/build_try_data.py. No other network request. */
(function () {
  "use strict";
  var data = window.QD_TRY;
  var root = document.getElementById("qd-try");
  if (!data || !root) return;

  var state = { method: "d", cuts: 2, shots: 2, seed: 0, joint: "on" };
  var extras = {
    a: "Mid-circuit measure: yes. Extra qubits: none. Best when no ancilla is available.",
    b: "Mid-circuit measure: no. Extra qubits: +2 per cut (one ancilla per side).",
    c: "Mid-circuit measure: no. Extra qubits: +1 per cut. The general wire cut.",
    d: "Mid-circuit measure: no. Extra qubits: +1 per cut. Each side runs once per setting of its own cuts."
  };

  function el(tag, attrs, text) {
    var node = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (key) { node.setAttribute(key, attrs[key]); });
    if (text !== undefined) node.textContent = text;
    return node;
  }
  function count(n) {
    if (n < 1e6) return n.toLocaleString("en-US");
    var e = Math.floor(Math.log10(n));
    return (n / Math.pow(10, e)).toFixed(1) + "e" + e;
  }
  function small(x) { return x < 1e-12 ? "below 1e-12" : x.toExponential(0); }
  function num(x, digits) { return (x >= 0 ? "+" : "−") + Math.abs(x).toFixed(digits); }
  function entry() { return data.runs[state.method][String(state.cuts)]; }
  /* Where joint cuts change nothing, one set of runs serves both. */
  function run() { var e = entry(); return e.states[state.joint] || e.states.on; }
  function jointMatters() { var from = data.joint_from[state.method]; return !!from && state.cuts >= from; }
  function shots() { return data.shots[state.shots]; }

  function choices(box, items, active, onPick) {
    items.forEach(function (item) {
      var button = el("button", { type: "button" }, item.label);
      button.dataset.value = item.value;
      button.addEventListener("click", function () { onPick(item.value); update(); });
      box.appendChild(button);
    });
    return function () {
      Array.prototype.forEach.call(box.children, function (b) {
        b.setAttribute("aria-pressed", String(b.dataset.value === String(active())));
      });
    };
  }

  /* Static skeleton */
  root.textContent = "";
  root.appendChild(el("span", { class: "qd-label" }, "1. Pick a method"));
  var methodBox = root.appendChild(el("div", { class: "qd-choices qd-methods" }));
  var markMethod = choices(methodBox, Object.keys(data.methods).map(function (m) {
    return { value: m, label: "(" + m + ") " + data.methods[m].label };
  }), function () { return state.method; }, function (v) { state.method = v; });
  var extra = root.appendChild(el("p", { class: "qd-note" }));

  root.appendChild(el("span", { class: "qd-label" }, "2. Pick the number of cuts, k"));
  var cutBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
  var markCuts = choices(cutBox, data.cuts.map(function (k) {
    return { value: String(k), label: "k = " + k };
  }), function () { return state.cuts; }, function (v) { state.cuts = Number(v); });

  var jointLabel = root.appendChild(el("span", { class: "qd-label" }, "3. Joint cuts"));
  var jointBox = root.appendChild(el("div", { class: "qd-choices qd-cuts" }));
  var markJoint = choices(jointBox, [{ value: "on", label: "On" }, { value: "off", label: "Off" }],
    function () { return state.joint; }, function (v) { state.joint = v; });
  var note = root.appendChild(el("p", { class: "qd-note" }));
  var stats = root.appendChild(el("div", { class: "qd-stats" }));
  var figure = root.appendChild(el("figure"));
  var image = figure.appendChild(el("img", { alt: "The circuit with its cuts and its fragments, drawn by qd.draw_cut" }));
  figure.appendChild(el("figcaption", { class: "qd-note" }, "Drawn by qd.draw_cut: each colour is a fragment, each pair of scissors a cut."));

  root.appendChild(el("span", { class: "qd-label" }, "4. The reconstruction equals the uncut result"));
  var check = root.appendChild(el("div", { class: "qd-check" }));

  root.appendChild(el("span", { class: "qd-label" }, "5. Now with shots"));
  var shotRow = root.appendChild(el("div", { class: "qd-shots" }));
  var slider = shotRow.appendChild(el("input", {
    type: "range", min: "0", max: String(data.shots.length - 1), step: "1", "aria-label": "Shots per subcircuit"
  }));
  slider.value = String(state.shots);
  slider.addEventListener("input", function () { state.shots = Number(slider.value); state.seed = 0; update(); });
  var again = shotRow.appendChild(el("button", { type: "button" }, "Another run"));
  again.addEventListener("click", function () { state.seed = (state.seed + 1) % data.seeds; update(); });
  var shotStats = root.appendChild(el("div", { class: "qd-stats" }));
  var plots = root.appendChild(el("div", { class: "qd-plots" }));
  var bars = root.appendChild(el("div"));

  root.appendChild(el("span", { class: "qd-label" }, "How the subcircuits grow with k"));
  var table = root.appendChild(el("table"));
  root.appendChild(el("p", { class: "qd-note" },
    "From qd.cutting_overhead, which runs nothing. Each cut also multiplies the shots needed in the worst case by gamma squared."));
  root.appendChild(el("p", { class: "qd-note" },
    "Nothing quantum runs in your browser. These are stored results of real runs: Qdislib " + data.version +
    " on Qiskit Aer, " + data.seeds + " runs per point, generated by docs/build_try_data.py. " +
    "The circuit has k + 1 qubits on each side, joined by k CZ gates. Run the real code below."));

  function stat(box, value, label) {
    var item = box.appendChild(el("div", { class: "qd-stat" }));
    item.appendChild(el("b", {}, value));
    item.appendChild(el("span", {}, label));
  }

  function plot(title, index) {
    var e = entry(), r = run(), W = 300, H = 170, L = 44, R = 18, T = 22, B = 30;
    var means = data.shots.map(function (s) {
      var runs = r.runs[String(s)];
      return runs.reduce(function (sum, x) {
        return sum + (index === 0 ? Math.abs(x[0] - e.uncut) : x[1]);
      }, 0) / runs.length;
    });
    var current = r.runs[String(shots())][state.seed];
    var now = index === 0 ? Math.abs(current[0] - e.uncut) : current[1];
    var low = Math.min.apply(null, means.concat([now]).map(function (v) { return Math.max(v, 1e-5); }));
    var high = Math.max.apply(null, means.concat([now]));
    var lo = Math.floor(Math.log10(low)), hi = Math.ceil(Math.log10(high));
    if (hi === lo) hi += 1;
    var x = function (i) { return L + i * (W - L - R) / (data.shots.length - 1); };
    var y = function (v) { return T + (hi - Math.log10(Math.max(v, 1e-5))) * (H - T - B) / (hi - lo); };
    var s = '<svg viewBox="0 0 ' + W + " " + H + '" role="img" aria-label="' + title + '">';
    s += '<text class="t" x="0" y="13">' + title + "</text>";
    for (var e = lo; e <= hi; e++) {
      s += '<line x1="' + L + '" x2="' + (W - R) + '" y1="' + y(Math.pow(10, e)) + '" y2="' + y(Math.pow(10, e)) + '" stroke="var(--qd-line)"/>';
      s += '<text x="' + (L - 5) + '" y="' + (y(Math.pow(10, e)) + 4) + '" text-anchor="end">1e' + e + "</text>";
    }
    data.shots.forEach(function (v, i) {
      s += '<text x="' + x(i) + '" y="' + (H - B + 15) + '" text-anchor="middle">2^' + Math.log2(v) + "</text>";
    });
    s += '<polyline fill="none" stroke="var(--qd-a)" stroke-width="2" points="' +
      means.map(function (v, i) { return x(i) + "," + y(v); }).join(" ") + '"/>';
    means.forEach(function (v, i) {
      s += '<circle cx="' + x(i) + '" cy="' + y(v) + '" r="4" fill="var(--qd-a)"/>';
    });
    s += '<rect x="' + (x(state.shots) - 6) + '" y="' + (y(now) - 6) + '" width="12" height="12" rx="2" fill="var(--qd-b)" stroke="var(--color-background-primary, #fff)" stroke-width="2"/>';
    s += '<text x="' + ((W + L - R) / 2) + '" y="' + (H - 2) + '" text-anchor="middle">shots per subcircuit</text></svg>';
    return s;
  }

  function update() {
    var e = entry(), r = run(), runs = r.runs[String(shots())], current = runs[state.seed];
    var form = jointMatters() ? state.joint : "on";
    markMethod();
    markCuts();
    markJoint();
    extra.textContent = extras[state.method];
    var joinable = !!data.joint_from[state.method];
    jointLabel.style.display = jointMatters() ? "" : "none";
    jointBox.style.display = jointMatters() ? "" : "none";
    stats.textContent = "";
    var formula = data.methods[state.method].formula[form].replace(/k/g, state.cuts);
    stat(stats, count(r.subcircuits), "subcircuits to run: " + formula);
    stat(stats, r.fragments.join(" + "), "qubits of the 2 fragments, from a " + e.qubits + "-qubit circuit");
    stat(stats, String(r.gamma), "gamma, the sampling overhead");
    var texts = [];
    if (jointMatters()) {
      var on = e.states.on, off = e.states.off;
      var line = "Joint cuts: " + count(on.subcircuits) + " subcircuits, without: " + count(off.subcircuits) + ", the form on the poster.";
      if (on.gamma !== off.gamma) line += " Gamma: " + on.gamma + " with them, " + off.gamma + " without.";
      texts.push(line + (state.method === "b"
        ? " The gates crossing one boundary share one ancilla per side."
        : " The wires a fragment sends are measured together."));
    } else if (joinable) {
      texts.push("With one cut there is nothing to join: joint cuts give the same circuits.");
    } else {
      texts.push("Joint cuts do not apply to this method: " + (state.method === "a" ? "each gate is cut on its own." : "they change the optimized wire cut, not this one."));
    }
    if ((state.method === "c" || state.method === "d") && state.cuts === 1) {
      texts.push("This circuit has complex gates. With only real ones before the cut, such as ry and cz, this wire is a golden cut point: its Y terms are zero and are left out, so " +
        data.golden.subcircuits + " subcircuits (" + data.golden.optimized + " optimized) and gamma " + data.golden.gamma + ".");
    }
    note.textContent = texts.join(" ");
    image.src = data.figures[(state.method === "d" ? "c" : state.method) + ":" + state.cuts];

    check.textContent = "";
    check.appendChild(el("p", {}, "Uncut circuit, exact statevector: <Z\u2026Z> = " + num(e.uncut, 6)));
    check.appendChild(el("p", {}, "Reconstructed from the " + count(r.subcircuits) + " subcircuits" + (jointMatters() ? (state.joint === "on" ? " with joint cuts" : " without joint cuts") : "") + ", exact: " + num(r.reconstructed, 6)));
    check.appendChild(el("p", {}, "Difference " + small(Math.abs(r.reconstructed - e.uncut)) +
      ", and total variation distance from the exact distribution " + small(r.tvd_exact) + ": rounding only."));

    slider.value = String(state.shots);
    again.textContent = "Another run (" + (state.seed + 1) + " of " + data.seeds + ")";
    shotStats.textContent = "";
    stat(shotStats, count(shots()), "shots per subcircuit");
    stat(shotStats, num(current[0], 4), "<Z…Z> reconstructed, exact " + num(e.uncut, 4));
    stat(shotStats, Math.abs(current[0] - e.uncut).toFixed(4), "absolute error of this run");
    stat(shotStats, current[1].toFixed(4), "TVD of this run from the exact distribution");

    plots.innerHTML = plot("Absolute error of <Z…Z>".replace("<", "&lt;").replace(">", "&gt;"), 0) + plot("Total variation distance", 1);

    bars.textContent = "";
    var legend = bars.appendChild(el("div", { class: "qd-legend" }));
    [["var(--qd-a)", "circle: mean of " + data.seeds + " runs"], ["var(--qd-b)", "square: this run"]].forEach(function (x) {
      var item = legend.appendChild(el("span"));
      var swatch = item.appendChild(el("span", { class: "qd-swatch" }));
      swatch.style.background = x[0];
      item.appendChild(document.createTextNode(x[1]));
    });
    bars.appendChild(el("span", { class: "qd-label" }, "The eight likeliest outcomes"));
    var legend2 = bars.appendChild(el("div", { class: "qd-legend" }));
    [["var(--qd-a)", "exact"], ["var(--qd-b)", "reconstructed, this run"]].forEach(function (x) {
      var item = legend2.appendChild(el("span"));
      var swatch = item.appendChild(el("span", { class: "qd-swatch" }));
      swatch.style.background = x[0];
      item.appendChild(document.createTextNode(x[1]));
    });
    var scale = Math.max.apply(null, e.exact_top.concat(current[2]));
    e.keys.forEach(function (key, i) {
      var row = bars.appendChild(el("div", { class: "qd-row" }));
      row.appendChild(el("code", {}, key));
      var cell = row.appendChild(el("div"));
      [["var(--qd-a)", e.exact_top[i]], ["var(--qd-b)", current[2][i]]].forEach(function (x) {
        var bar = cell.appendChild(el("div", { class: "qd-bar" }));
        bar.style.background = x[0];
        bar.style.width = Math.max(0, 100 * x[1] / scale) + "%";
      });
      row.appendChild(el("span", {}, current[2][i].toFixed(3)));
    });

    var html = "<thead><tr><th>k</th>" + Object.keys(data.methods).map(function (m) {
      return "<th>(" + m + ")<br>" + data.methods[m].formula[state.joint] + "</th>";
    }).join("") + "</tr></thead><tbody>";
    data.table_cuts.forEach(function (k, i) {
      html += "<tr><td>" + k + "</td>" + Object.keys(data.methods).map(function (m) {
        var on = m === state.method && k === state.cuts;
        return "<td" + (on ? ' class="on"' : "") + ">" + count(data.table[state.joint][m][i]) + "</td>";
      }).join("") + "</tr>";
    });
    table.innerHTML = html + "</tbody>";
  }

  update();
})();
