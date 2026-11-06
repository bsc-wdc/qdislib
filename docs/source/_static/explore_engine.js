/* The circuit playground's engine on the Explore page: a statevector simulator and the
   decompositions of a gate cut, a Hadamard gate cut and a wire cut. No dependencies; it
   runs in the browser and, for docs/check_explore_engine.py, in Node. That script
   compares what it returns with real Qdislib on random circuits.

   A circuit is { n, ops: [{ name, qubits, theta }] } in the order Qdislib numbers its
   gates from 1. A cut is { kind: "gate" | "hadamard", op } or
   { kind: "wire", qubit, from, to }, the indices of the two gates around the wire. */
(function (root) {
  "use strict";

  var S2 = Math.SQRT1_2;
  var TWO_QUBIT = { cz: 1, cx: 1, rzz: 1 };

  function matrix(name, theta) {
    var c, s;
    switch (name) {
      case "h": return [S2, 0, S2, 0, S2, 0, -S2, 0];
      case "x": return [0, 0, 1, 0, 1, 0, 0, 0];
      case "z": return [1, 0, 0, 0, 0, 0, -1, 0];
      case "s": return [1, 0, 0, 0, 0, 0, 0, 1];
      case "sdg": return [1, 0, 0, 0, 0, 0, 0, -1];
      case "t": return [1, 0, 0, 0, 0, 0, S2, S2];
      case "rx":
        c = Math.cos(theta / 2); s = Math.sin(theta / 2);
        return [c, 0, 0, -s, 0, -s, c, 0];
      case "ry":
        c = Math.cos(theta / 2); s = Math.sin(theta / 2);
        return [c, 0, -s, 0, s, 0, c, 0];
      case "rz":
        c = Math.cos(theta / 2); s = Math.sin(theta / 2);
        return [c, -s, 0, 0, 0, 0, c, s];
    }
    throw new Error("Unknown gate " + name);
  }

  function State(width) {
    this.width = width;
    this.re = new Float64Array(1 << width);
    this.im = new Float64Array(1 << width);
    this.re[0] = 1;
  }
  State.prototype.one = function (m, q) {
    var bit = 1 << q, re = this.re, im = this.im, size = re.length;
    for (var i = 0; i < size; i++) {
      if (i & bit) continue;
      var j = i | bit, ar = re[i], ai = im[i], br = re[j], bi = im[j];
      re[i] = m[0] * ar - m[1] * ai + m[2] * br - m[3] * bi;
      im[i] = m[0] * ai + m[1] * ar + m[2] * bi + m[3] * br;
      re[j] = m[4] * ar - m[5] * ai + m[6] * br - m[7] * bi;
      im[j] = m[4] * ai + m[5] * ar + m[6] * bi + m[7] * br;
    }
  };
  State.prototype.gate = function (name, q, theta) { this.one(matrix(name, theta), q); };
  State.prototype.cz = function (a, b) {
    var mask = (1 << a) | (1 << b), re = this.re, im = this.im;
    for (var i = 0; i < re.length; i++) if ((i & mask) === mask) { re[i] = -re[i]; im[i] = -im[i]; }
  };
  State.prototype.cx = function (c, t) {
    var cb = 1 << c, tb = 1 << t, re = this.re, im = this.im;
    for (var i = 0; i < re.length; i++) {
      if ((i & cb) && !(i & tb)) {
        var j = i | tb, r = re[i], m = im[i];
        re[i] = re[j]; im[i] = im[j]; re[j] = r; im[j] = m;
      }
    }
  };
  State.prototype.rzz = function (a, b, theta) {
    var c = Math.cos(theta / 2), s = Math.sin(theta / 2), re = this.re, im = this.im;
    for (var i = 0; i < re.length; i++) {
      var odd = ((i >> a) ^ (i >> b)) & 1, sign = odd ? 1 : -1, r = re[i], m = im[i];
      re[i] = c * r - sign * s * m;
      im[i] = c * m + sign * s * r;
    }
  };
  /* A Hadamard test's select: |0><0| (x) U_a + |1><1| (x) U_b, with each U = factor * Z^power. */
  State.prototype.select = function (anc, q, u0, u1) {
    var re = this.re, im = this.im;
    for (var i = 0; i < re.length; i++) {
      var u = (i >> anc) & 1 ? u1 : u0, sign = u.z && ((i >> q) & 1) ? -1 : 1;
      var fr = sign * u.re, fi = sign * u.im, r = re[i], m = im[i];
      re[i] = fr * r - fi * m;
      im[i] = fr * m + fi * r;
    }
  };
  /* <Z on the bits of mask> */
  State.prototype.parity = function (mask) {
    var total = 0, re = this.re, im = this.im;
    for (var i = 0; i < re.length; i++) {
      var x = i & mask, odd = 0;
      while (x) { odd ^= 1; x &= x - 1; }
      var p = re[i] * re[i] + im[i] * im[i];
      total += odd ? -p : p;
    }
    return total;
  };

  function applyOp(state, op, a, b) {
    if (op.name === "cz") state.cz(a, b);
    else if (op.name === "cx") state.cx(a, b);
    else if (op.name === "rzz") state.rzz(a, b, op.theta);
    else state.gate(op.name, a, op.theta);
  }

  /* The basis changes that read an X or a Y as a Z, as gates after the circuit's own. */
  function basisOps(observable) {
    var out = [];
    observable.split("").forEach(function (p, q) {
      if (p === "X") out.push({ name: "h", qubits: [q], basis: true });
      if (p === "Y") {
        out.push({ name: "sdg", qubits: [q], basis: true });
        out.push({ name: "h", qubits: [q], basis: true });
      }
    });
    return out;
  }

  function exact(circuit, observable) {
    var state = new State(circuit.n), mask = 0;
    circuit.ops.forEach(function (op) { applyOp(state, op, op.qubits[0], op.qubits[1]); });
    basisOps(observable).forEach(function (op) { applyOp(state, op, op.qubits[0]); });
    observable.split("").forEach(function (p, q) { if (p !== "I") mask |= 1 << q; });
    return state.parity(mask);
  }

  function isReal(op) {
    switch (op.name) {
      case "h": case "x": case "z": case "cz": case "cx": case "ry": return true;
      case "s": case "t": case "sdg": return false;
      default: return Math.abs(Math.sin(op.theta / 2)) < 1e-12;
    }
  }

  /* The rzz a cut gate is, with what surrounds it: CZ = (S (x) S) rzz(-pi/2) up to a
     phase, and CX is a CZ between Hadamards on its target. */
  function cutAngle(op) { return op.name === "rzz" ? op.theta : -Math.PI / 2; }

  function gateGamma(op) { return 1 + 2 * Math.abs(Math.sin(cutAngle(op))); }

  /* Fragments: each qubit's line is split at its wire cuts into segments, and segments
     joined by an uncut two-qubit gate are one fragment. */
  function plan(circuit, cuts, observable) {
    var n = circuit.n, ops = circuit.ops.concat(basisOps(observable));
    var gateCut = {}, wireAfter = {};
    cuts.forEach(function (cut) {
      if (cut.kind === "wire") wireAfter[cut.qubit + ":" + cut.from] = cut;
      else gateCut[cut.op] = cut;
    });

    var segments = [], current = [], opSegments = ops.map(function () { return []; });
    function newSegment(q, ordinal) {
      segments.push({ id: segments.length, qubit: q, ordinal: ordinal, ops: [], incoming: null, outgoing: null });
      return segments[segments.length - 1];
    }
    for (var q = 0; q < n; q++) current[q] = newSegment(q, 0);
    ops.forEach(function (op, i) {
      op.qubits.forEach(function (qubit, side) {
        var segment = current[qubit];
        segment.ops.push(i);
        opSegments[i][side] = segment.id;
        var cut = wireAfter[qubit + ":" + i];
        if (cut) {
          var next = newSegment(qubit, segment.ordinal + 1);
          segment.outgoing = cut;
          next.incoming = cut;
          cut.source = segment.id;
          cut.target = next.id;
          current[qubit] = next;
        }
      });
    });

    var parent = segments.map(function (s) { return s.id; });
    function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
    ops.forEach(function (op, i) {
      if (op.qubits.length === 2 && !gateCut[i]) parent[find(opSegments[i][0])] = find(opSegments[i][1]);
    });

    var byRoot = {}, fragments = [];
    segments.forEach(function (s) {
      if (!s.ops.length && !s.incoming && !s.outgoing) return; // an idle line
      var r = find(s.id);
      if (!byRoot[r]) {
        byRoot[r] = { id: fragments.length, segments: [], sides: [], measures: [], preps: [] };
        fragments.push(byRoot[r]);
      }
      byRoot[r].segments.push(s);
      s.fragment = byRoot[r].id;
    });

    var invalid = [];
    cuts.forEach(function (cut) {
      var a = cut.kind === "wire" ? cut.source : opSegments[cut.op][0];
      var b = cut.kind === "wire" ? cut.target : opSegments[cut.op][1];
      if (find(a) === find(b)) invalid.push(cut);
    });

    var wires = cuts.filter(function (c) { return c.kind === "wire"; });
    if (!invalid.length) {
      cuts.forEach(function (cut, index) {
        if (cut.kind === "wire") {
          fragments[segments[cut.source].fragment].measures.push(index);
          fragments[segments[cut.target].fragment].preps.push(index);
        } else {
          [0, 1].forEach(function (side) {
            fragments[segments[opSegments[cut.op][side]].fragment].sides.push({ cut: index, side: side });
          });
        }
      });
    }

    /* Golden wires: as Qdislib's _golden_wires. */
    var golden = {};
    if (!invalid.length) {
      var real = fragments.map(function (f) {
        return !f.sides.length && f.measures.length === 1 && f.segments.every(function (s) {
          return s.ops.every(function (i) { return isReal(ops[i]); });
        });
      });
      var changed = true;
      while (changed) {
        changed = false;
        fragments.forEach(function (f, i) {
          if (real[i] && f.preps.every(function (c) { return golden[c]; }) && !golden[f.measures[0]]) {
            golden[f.measures[0]] = true;
            changed = true;
          }
        });
      }
    }

    var gamma = 1, configurations = 1, hasHadamard = false;
    cuts.forEach(function (cut, index) {
      if (cut.kind === "wire") {
        gamma *= golden[index] ? 3 : 4;
        configurations *= golden[index] ? 6 : 8;
      } else {
        gamma *= gateGamma(ops[cut.op]);
        configurations *= 6;
        if (cut.kind === "hadamard") hasHadamard = true;
      }
    });
    if (!fragments.length) fragments.push({ id: 0, segments: [], sides: [], measures: [], preps: [] });

    /* Optimized counts, with joint wire measurement: Qdislib's num_units. */
    var optimized = 0;
    fragments.forEach(function (f) {
      var count = 1;
      f.sides.forEach(function () { count *= 5; });
      f.preps.forEach(function (c) { count *= golden[c] ? 3 : 4; });
      if (f.measures.length === 1 && golden[f.measures[0]]) count *= 2;
      else if (f.measures.length) count *= Math.pow(2, f.measures.length) + 1;
      optimized += count;
    });

    /* Hadamard cuts across one boundary, as joint=True runs them: for k rzz with
       coefficients |c|, |s|, 2^k product terms, 2 r^2 - r settings for r terms, and
       gamma = 2 (sum |c_i|)^2 - sum c_i^2 = 2 prod(1 + |sin theta|) - 1. */
    var joint = null;
    if (!invalid.length && hasHadamard) {
      var groups = {};
      cuts.forEach(function (cut) {
        if (cut.kind !== "hadamard") return;
        var a = segments[opSegments[cut.op][0]].fragment, b = segments[opSegments[cut.op][1]].fragment;
        var key = Math.min(a, b) + ":" + Math.max(a, b);
        (groups[key] = groups[key] || []).push(gateGamma(ops[cut.op]));
      });
      var jointGamma = 1, jointConfigurations = 1, shared = 0;
      cuts.forEach(function (cut, index) {
        if (cut.kind === "gate") { jointGamma *= gateGamma(ops[cut.op]); jointConfigurations *= 6; }
        if (cut.kind === "wire") { jointGamma *= golden[index] ? 3 : 4; jointConfigurations *= golden[index] ? 6 : 8; }
      });
      Object.keys(groups).forEach(function (key) {
        var product = 1, r = Math.pow(2, groups[key].length);
        groups[key].forEach(function (g) { product *= 1 + (g - 1) / 2; });
        jointGamma *= 2 * product - 1;
        jointConfigurations *= 2 * r * r - r;
        if (groups[key].length > 1) shared += groups[key].length;
      });
      joint = { shared: shared, gamma: jointGamma, configurations: jointConfigurations, subcircuits: jointConfigurations * fragments.length };
    }

    var widths = fragments.map(function (f) {
      var ancillas = 0;
      f.sides.forEach(function (side) { if (cuts[side.cut].kind === "hadamard") ancillas += 1; });
      return f.segments.length + ancillas;
    });

    return {
      ops: ops, cuts: cuts, segments: segments, opSegments: opSegments, fragments: fragments,
      invalid: invalid, golden: golden, gamma: gamma, configurations: configurations,
      subcircuits: configurations * fragments.length, optimized: hasHadamard ? null : optimized,
      widths: widths, wires: wires.length, joint: joint
    };
  }

  /* One cut's terms: the action on each side, and the weight. */
  function gateTerms(op) {
    var theta = cutAngle(op), c = Math.cos(theta / 2), s = Math.sin(theta / 2);
    return [
      { w: c * c, a: "i", b: "i" },
      { w: s * s, a: "z", b: "z" },
      { w: c * s, a: "m", b: "r+" },
      { w: -c * s, a: "m", b: "r-" },
      { w: c * s, a: "r+", b: "m" },
      { w: -c * s, a: "r-", b: "m" }
    ];
  }
  /* The Hadamard test's settings (i, j, g), as Qdislib lists them, of
     rzz(theta) = |c| I (x) I + |s| (sign(c s) (-i Z)) (x) Z times the sign of c. */
  function hadamardTerms(op) {
    var theta = cutAngle(op), c = Math.cos(theta / 2), s = Math.sin(theta / 2);
    var w = [Math.abs(c), Math.abs(s)];
    var out = [];
    [[0, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1]].forEach(function (t) {
      var setting = { i: t[0], j: t[1], g: t[2], phase: (c * s >= 0 ? -1 : 1) };
      out.push({ w: w[t[0]] * w[t[1]] * (t[2] ? -1 : 1), a: { had: setting, side: 0 }, b: { had: setting, side: 1 } });
    });
    return out;
  }
  var WIRE_TERMS = [
    { w: 0.5, read: "I", prep: "0" }, { w: 0.5, read: "I", prep: "1" },
    { w: 0.5, read: "X", prep: "+" }, { w: -0.5, read: "X", prep: "-" },
    { w: 0.5, read: "Y", prep: "+i" }, { w: -0.5, read: "Y", prep: "-i" },
    { w: 0.5, read: "Z", prep: "0" }, { w: -0.5, read: "Z", prep: "1" }
  ];

  function termsOf(p, cut, index) {
    if (cut.kind === "gate") return gateTerms(p.ops[cut.op]);
    if (cut.kind === "hadamard") return hadamardTerms(p.ops[cut.op]);
    return p.golden[index] ? WIRE_TERMS.filter(function (t) { return t.read !== "Y"; }) : WIRE_TERMS;
  }

  /* <Pauli read> of one fragment with one term chosen for each cut that touches it. */
  function fragmentValue(p, fragment, observable, choice) {
    var cuts = p.cuts, ops = p.ops, local = {}, width = 0, ordered = [];
    fragment.segments.forEach(function (s) { local[s.id] = width++; });
    fragment.sides.forEach(function (side) {
      var cut = cuts[side.cut], term = choice[side.cut], action = side.side ? term.b : term.a;
      if (cut.kind === "gate" ? action === "m" : action.had.i !== action.had.j) side.ancilla = width++;
      else side.ancilla = null;
    });
    if (width > 24) throw new Error("Fragment too wide");
    var state = new State(width), mask = 0;

    function sideOfOp(i, s) {
      for (var k = 0; k < fragment.sides.length; k++) {
        var side = fragment.sides[k];
        if (cuts[side.cut].op === i && side.side === s) return side;
      }
      return null;
    }

    /* Wire states prepared at the start of a segment. */
    fragment.segments.forEach(function (s) {
      if (!s.incoming) return;
      var term = choice[cuts.indexOf(s.incoming)], q = local[s.id];
      if (term.prep === "1") state.gate("x", q);
      if (term.prep === "+" || term.prep === "+i") state.gate("h", q);
      if (term.prep === "-" || term.prep === "-i") { state.gate("x", q); state.gate("h", q); }
      if (term.prep === "+i" || term.prep === "-i") state.gate("s", q);
    });
    fragment.sides.forEach(function (side) {
      if (side.ancilla === null || cuts[side.cut].kind !== "hadamard") return;
      var action = side.side ? choice[side.cut].b : choice[side.cut].a;
      state.gate("h", side.ancilla);
      if (action.had.g) state.gate("s", side.ancilla);
    });

    /* Gates in time order. Two segments of one fragment cannot be out of step, as a
       fragment's gates keep the circuit's order. */
    var touched = {};
    fragment.segments.forEach(function (s) { s.ops.forEach(function (i) { touched[i] = true; }); });
    Object.keys(touched).map(Number).sort(function (x, y) { return x - y; }).forEach(function (i) {
      var op = ops[i], cut = p.gateCutOf ? p.gateCutOf[i] : null;
      var index = -1;
      for (var k = 0; k < cuts.length; k++) if (cuts[k].kind !== "wire" && cuts[k].op === i) index = k;
      if (index < 0) {
        var segA = p.opSegments[i][0], segB = p.opSegments[i][1];
        applyOp(state, op, local[segA], op.qubits.length === 2 ? local[segB] : undefined);
        return;
      }
      [0, 1].forEach(function (s) {
        var seg = p.opSegments[i][s];
        if (!(seg in local)) return;
        var q = local[seg], side = sideOfOp(i, s), action = s ? choice[index].b : choice[index].a;
        /* CX is a CZ between Hadamards on its target. */
        if (op.name === "cx" && s === 1) state.gate("h", q);
        if (cuts[index].kind === "gate") {
          if (action === "z") state.gate("z", q);
          if (action === "r+") state.gate("rz", q, Math.PI / 2);
          if (action === "r-") state.gate("rz", q, -Math.PI / 2);
          if (action === "m") { state.cx(q, side.ancilla); mask |= 1 << side.ancilla; }
        } else {
          var had = action.had;
          /* The two unitaries of this side: side 0 carries the phase of -i. */
          var unitaries = [{ re: 1, im: 0, z: false }, s === 0 ? { re: 0, im: had.phase, z: true } : { re: 1, im: 0, z: true }];
          if (had.i === had.j) {
            if (had.i === 1) state.select(q, q, unitaries[1], unitaries[1]);
          } else {
            state.select(side.ancilla, q, unitaries[had.i], unitaries[had.j]);
          }
        }
        /* The local gates that make a CZ of the rzz(-pi/2). */
        if (op.name !== "rzz") state.gate("s", q);
        if (op.name === "cx" && s === 1) state.gate("h", q);
      });
    });

    fragment.sides.forEach(function (side) {
      if (side.ancilla === null || cuts[side.cut].kind !== "hadamard") return;
      state.gate("h", side.ancilla);
      mask |= 1 << side.ancilla;
    });
    /* What is read: the observable on a line's last segment, a wire's Pauli on the others. */
    fragment.segments.forEach(function (s) {
      var q = local[s.id];
      if (s.outgoing) {
        var read = choice[cuts.indexOf(s.outgoing)].read;
        if (read === "X") state.gate("h", q);
        if (read === "Y") { state.gate("sdg", q); state.gate("h", q); }
        if (read !== "I") mask |= 1 << q;
      } else if (observable.charAt(s.qubit) !== "I") {
        mask |= 1 << q;
      }
    });
    return state.parity(mask);
  }

  /* The reconstruction: the sum over every term of every cut of the weight and of the
     product of the fragments' values, each fragment run once per setting of its own cuts.
     An idle qubit (no gates, no basis change) is |0>: Z and I both read +1. */
  function prepare(circuit, cuts, observable) {
    var p = plan(circuit, cuts, observable);
    if (p.invalid.length) throw new Error("A cut does not separate its two sides");
    var terms = cuts.map(function (cut, i) { return termsOf(p, cut, i); });
    var memo = p.fragments.map(function () { return {}; });
    var touching = p.fragments.map(function (f) {
      var list = [];
      f.sides.forEach(function (s) { if (list.indexOf(s.cut) < 0) list.push(s.cut); });
      f.measures.concat(f.preps).forEach(function (c) { if (list.indexOf(c) < 0) list.push(c); });
      return list;
    });
    var digits = cuts.map(function () { return 0; }), choice = [], runs = 0;

    function key(f) {
      var out = "";
      for (var k = 0; k < touching[f].length; k++) out += digits[touching[f][k]] + ",";
      return out;
    }
    function exactValue(f) {
      var k = key(f);
      if (k in memo[f]) return memo[f][k];
      runs++;
      return (memo[f][k] = p.fragments[f].segments.length ? fragmentValue(p, p.fragments[f], observable, choice) : 1);
    }
    /* draws: null for the exact value, or per fragment a map of this draw's values. */
    function sum(leaf) {
      var total = 0;
      (function walk(i, weight) {
        if (i === cuts.length) {
          var product = weight;
          for (var f = 0; f < p.fragments.length; f++) product *= leaf(f);
          total += product;
          return;
        }
        for (var d = 0; d < terms[i].length; d++) {
          digits[i] = d;
          choice[i] = terms[i][d];
          walk(i + 1, weight * terms[i][d].w);
        }
      })(0, 1);
      return total;
    }
    var exactResult = sum(exactValue);
    return {
      plan: p,
      value: exactResult,
      runs: runs,
      /* One run with shots: every distinct fragment circuit gets `shots` shots. */
      draw: function (shots, rng) {
        var draws = p.fragments.map(function () { return {}; });
        return sum(function (f) {
          var k = key(f);
          if (k in draws[f]) return draws[f][k];
          return (draws[f][k] = sampled(memo[f][k], shots, rng));
        });
      }
    };
  }

  function reconstruct(circuit, cuts, observable) {
    return prepare(circuit, cuts, observable);
  }

  /* A run of `shots` shots of a circuit whose per-shot result is +1 or -1, with mean v. */
  function sampled(v, shots, rng) {
    var prob = (1 + v) / 2, ones = 0;
    if (prob <= 0) return -1;
    if (prob >= 1) return 1;
    if (shots <= 400) {
      for (var i = 0; i < shots; i++) if (rng() < prob) ones++;
    } else {
      var u = 1 - rng(), w = rng();
      var normal = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * w);
      ones = Math.max(0, Math.min(shots, Math.round(shots * prob + normal * Math.sqrt(shots * prob * (1 - prob)))));
    }
    return (2 * ones) / shots - 1;
  }

  function mulberry(seed) {
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  var api = { exact: exact, plan: plan, prepare: prepare, reconstruct: reconstruct, mulberry: mulberry, TWO_QUBIT: TWO_QUBIT, terms: termsOf };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.QD_ENGINE = api;
})(typeof window !== "undefined" ? window : this);
