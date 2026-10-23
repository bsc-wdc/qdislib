# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- **Many cuts reconstruct in seconds: the reconstruction is a contraction.** Each
  cut (a measure-and-prepare gate cut, a Hadamard group, a wire) is a bond and each
  fragment a tensor over the bonds it touches; contracting them replaces the loop
  over every configuration. Ten cuts in a line now reconstruct exactly in about
  0.2-9 s (`6**10` to `8**10` configurations, which did not finish before), and the
  cost follows the fragment graph, not the number of cuts. An expectation value is
  one task; a distribution reduces each fragment's results in its own task first.
- **`num_samples` on every cut type.** `gate_cutting`, `wire_cutting` and
  `circuit_cutting` (with their `_sampling` and three-step forms) can draw `N`
  configurations with probability `|weight| / gamma` instead of running them all,
  and run only the fragment circuits the draws use, with `shots` per draw -- the
  quasiprobability Monte-Carlo estimate, so the shots follow the weights. Unbiased
  (error ~ `gamma / sqrt(N * shots)`); the way to run cuts with too many
  configurations to list (`6**n`, `8**n`), and much cheaper than equal shots when
  weights differ (Hadamard cuts of rotations: about 1,600x fewer shots on a QAOA-100
  cut). With equal weights and few configurations, running them all with the same
  total shots stays slightly more accurate (about 15% lower error). Three-step lists
  give each subcircuit a `draws` count to size its shots. (`gate_cutting_hadamard`
  already had `num_samples`; this brings it to every cut type.)
- **`circuit_cutting`: any combination of cuts on one circuit.** Measure-and-prepare
  gate cuts, Hadamard gate cuts and wire cuts, in any number, with the same four forms
  as every other cut: `circuit_cutting` (expectation value, X/Y and weighted Pauli
  sums), `circuit_cutting_subcircuits`, `circuit_cutting_sampling` (output
  distribution) and `circuit_cutting_sampling_subcircuits`. `optimized=True` runs each
  fragment once per setting of the cuts that touch it, so the circuit count is a sum
  over fragments instead of a product over cuts (e.g. 166 instead of 6,912 for two
  gate cuts and two wire cuts across three fragments); this also brings each side of a
  measure-and-prepare gate cut from 6 circuits to 5. `joint` covers the joint wire
  measurement and the joint Hadamard cut. It replaces `gate_wire_cutting` (deprecated).
- **`draw_cut` draws the circuit and what the cut does to it.** A circuit diagram
  (qubit lines, gates in time order, each with its index -- the name to cut) with
  every fragment the cut leaves in its own colour and every cut marked: measure-and-
  prepare gate cuts, Hadamard gate cuts and wire cuts, in any mix, as a list or a
  `{"gate", "hadamard", "wire"}` dict. A cut that does not separate its two sides is
  drawn and named instead of raising; the title gives the fragments, the cuts and
  `gamma`. Takes an `ax=` to draw into and returns the `Axes`. It used to draw the
  gate DAG with a spring layout, which placed gates at random and only knew a gate
  list or a wire list.
- **`find_cut(mixed=True)` searches every kind of cut, priced as it runs.** It
  minimises the sampling overhead `gamma` (shots grow as `gamma**2`) over gate,
  Hadamard and wire cuts together:
  - each gate cut takes the cheaper of measure-and-prepare (3 per CZ of the gate's
    decomposition) and the Hadamard test, whose `gamma` follows the gate's angle;
  - the Hadamard cuts between two fragments are priced jointly;
  - wire cuts split a qubit's timeline where that is cheaper.

  Every fragment fits `max_qubits` with its Hadamard ancillas, every cut separates its
  two sides, and the search (a local search over qubit timelines, seeded by the
  classic partitioners and, up to 40 qubits, an exact gate-partition MILP with
  SciPy/HiGHS) returns the best cut found within `time_limit` (10 s by default; about
  1-8 s at 100 qubits). It returns `{"gate", "hadamard", "wire"}` for
  `circuit_cutting`. On QAOA circuits it lowers `gamma` from 729 to about 32 (10
  qubits) and from 10^13.4 to 10^4.4 (100 qubits, two halves); on circuits where
  measure-and-prepare is already optimal it returns the same cut. New `find_cut`
  options: `hadamard`, `joint`, `exact`, `time_limit`. Research notes and benchmarks
  behind it: Fröhler et al. 2608.05287, Brandhofer et al. 2308.09567.
- **`cutting_overhead(circuit, cuts, ...)`** reports the circuits to run, the
  configurations, the fragments and `gamma` of any cut or mix, without running it.
- **Exact results on every cut type**: `shots=None` is a statevector reconstruction
  (on qibo or Qiskit) for the expectation value and the distribution, all-in-one and
  three-step; mid-circuit measurements are deferred to an extra qubit.
- **Output distributions with the Hadamard cut (`gate_cutting_hadamard_sampling`,
  `gate_cutting_hadamard_sampling_subcircuits`).** Every fragment is read in the Z
  basis, so one run per fragment gives its correlators for every output subset; their
  weighted products give every `<Z^A>` and the inverse Walsh-Hadamard transform the
  distribution, at the cost of one expectation value. Exact with `shots=None`;
  works with `joint` and `num_samples`. Every cut type now has both an
  expectation value and an output distribution, each in all-in-one and three-step form.
- **Joint Hadamard gate cut (`joint=True`).** `gate_cutting_hadamard` and its
  `_subcircuits` / `_sampling` variants can cut all the gates crossing one boundary
  together, with one shared ancilla per fragment: `gamma = 2**(k+1) - 1` instead of
  `3**k` for k CZ/CX, `2**k (2**(k+1) - 1)` settings instead of `6**k`, and one
  ancilla per side instead of k. Works for gates at different times and for general
  (KAK) gates, needs no communication. On by default (`joint=False` cuts each gate on
  its own).
- **The Hadamard cut measures X and Y and takes weighted Pauli sums**, like
  `gate_cutting`: X/Y qubits are rotated into the Z basis before cutting, and a sum
  is cut term by term (`*_subcircuits` return a `CombinedCutSubcircuits`).
- **Optimized wire cutting measures cut wires jointly.** With `optimized=True`, the
  wires that leave the same fragment are measured together: `2**k + 1` measurement
  circuits for k wires instead of `3**k`, each reading a commuting group of Paulis
  (mutually unbiased bases from GF(2**k)), every Pauli once. Applies to
  `wire_cutting`, `wire_cutting_subcircuits` and `wire_cutting_sampling`, for mixed
  (chain) fragments too. Some of these circuits add CZ gates between the cut qubits;
  `joint=False` keeps the previous per-wire measurement, whose results are unchanged.
- **`CombinedCutSubcircuits`**: the subcircuit container `*_subcircuits` return
  when a circuit has several independent pieces or the observable is a weighted
  Pauli sum. It is still one flat list to evaluate; `reconstruct` multiplies the
  pieces' results or sums the terms with their coefficients.
- **The IQM backend can run a cut circuit on IQM hardware.** `IQMBackend` now
  implements the cutting hooks (`evaluate_expectation` for the all-in-one path,
  `evaluate_subcircuit` for the three-step one), reconstructing each subcircuit's
  expectation with the same Z/I + cut-point-sign rule the simulators use, so a cut
  evaluated on the device recombines exactly as it does on Aer. Step 2 of the
  three-step workflow takes the backend directly::

      iqm = build_backend("IQM_BACKEND", BackendConfig(quantum_computer="garnet"))
      subs = qd.gate_cutting_subcircuits(circuit, ["CZ_3"], "ZZZZ")
      values = [qd.expectation_value(s, shots=2000, backend=iqm) for s in subs]
      result = subs.reconstruct(values)

  Verified end to end on IQM Resonance (Garnet, 20 qubits): a 4-qubit circuit cut
  into 12 subcircuits reconstructed `<ZZZZ> = 0.852` against an exact 1.0, the gap
  being hardware noise. Garnet accepts the mid-circuit measurement the gate-cut
  decomposition needs, so no fallback is required. Expectations remain Z/I-only
  (an X or Y term raises), and `gate_cutting(qpu=True, ...)` still routes only to
  `MN_Ona` / `IBM_Quantum` -- the three-step workflow is the IQM route.

### Fixed
- **`sync=False` returns a future from every cutting function.** Under PyCOMPSs,
  several paths reconstructed on the master from futures or waited anyway:
  `gate_cutting` / `wire_cutting` with shots (`TypeError` summing futures),
  `wire_cutting(optimized=True)`, the exact Hadamard cut, weighted Pauli sums in
  `wire_cutting` (one wait per term), the three-step exact runs (task converters
  called on the master) and the Hadamard three-step `reconstruct`. A strict fake
  PyCOMPSs in the tests now catches any of these.
- **Exact runs are parallel.** `shots=None` fragment circuits ran on the master one
  after the other; each is now a task, in the all-in-one and three-step forms.
- **`gate_cutting_hadamard` on Qiskit circuits** with `u`, `p`, `cp` or `cx` gates
  failed (`qibo.gates has no gate class 'u'`).
- **Large joint Hadamard groups no longer hang.** A joint group of k gates has
  `2**k` product terms and `2 * 4**k - 2**k` settings; they were written out as Python
  lists, so `cutting_overhead`, `draw_cut` and `num_samples` hung on cuts like the
  13-gate group `find_cut` proposes for a 100-qubit QAOA circuit. The terms are now
  built on demand from their index, `gamma` and the setting counts come in closed
  form, and settings are drawn gate by gate (also for
  `gate_cutting_hadamard(..., num_samples=N)`): `cutting_overhead` on that cut takes
  0.06 s.
- **METIS partitioning in `find_cut` ignored repeated gates.** Two gates on the same
  qubit pair (or a gate-graph edge severing two wires) counted as one edge, so METIS
  minimised the wrong cut; on a ring with heavy and light links it cut the heavy ones.
  Edges now carry their multiplicity, and METIS and spectral clustering check
  `max_cuts` against the weighted count.
- **`cutting_overhead` took minutes for about a dozen cuts**: it listed every
  configuration (6**g * 8**w). Its counts now come from closed forms (0.04 s for 12
  CZ cuts, against minutes).
- **Gate cutting on IBM Quantum hardware returned about 0 for any measured cut.**
  The IBM reader looked only at the final measurements, so it dropped every
  cut's decomposition sign, and it ran the sampler with its default shot count
  while dividing by `shots`. It now reads the cut-point registers like the
  simulators and passes `shots`; checked through SamplerV2 in local testing mode.
- **Several gate cuts on one qubit lost all but one sign.** Two cuts on the same
  qubit pair, or a general gate decomposed into several CZs (SWAP, RXX, ...),
  kept only the last cut-point outcome on the Qiskit path (two CZ cuts: +0.07
  instead of -0.95), and on the sampling paths of every backend. Each cut now
  keeps its own outcome. Reported by David Ryan (bsc-wdc/qdislib#3).
- **A circuit's final measurements broke gate cutting.** Qiskit's `measure_all()`
  was read as the gate cut's own cut-point measurement (ZZ gave 1.0 instead of
  -0.009), and a final Qibo `M(0, 1)` joined all its qubits so a valid cut was
  rejected. Final measurements are now left out of the DAG; measurements followed
  by more gates are kept. Reported by David Ryan (bsc-wdc/qdislib#4).
- **The three-step workflow ignored a weighted Pauli sum** and gave the all-Z
  value; `gate_cutting_hadamard` silently read a list, or an X/Y, as identity. It
  now measures X/Y (rotating them into Z) and splits a Pauli sum into terms.
  Reported by David Ryan (bsc-wdc/qdislib#6).
- **Two wire cuts in series on one wire gave a wrong value** (0.000 instead of
  -0.226): the second cut used the qubit index from before the first one moved
  the wire. Reported by David Ryan (bsc-wdc/qdislib#7).
- **`*_subcircuits` and `cache=` were wrong when part of the circuit was not
  cut**, or when several pieces were cut (+0.75 instead of +0.17): the pieces'
  subcircuits were grouped by guessing from the list length. Reported by David
  Ryan (bsc-wdc/qdislib#9).
- **`find_cut` wire mode ignored `max_qubits` and `forbidden_qubits`.** It
  counted gates instead of qubits, and never matched a protected wire. Every wire
  candidate is now checked as `wire_cutting` will apply it. Reported by David Ryan
  (bsc-wdc/qdislib#10).
- **The circuit cache served one circuit's result for the same circuit on other
  qubits** (X on qubit 0 vs qubit 1): the ZX key now tells qubits apart. Entries
  stored under the old keys are no longer hit. Reported by David Ryan
  (bsc-wdc/qdislib#11).
- **A Qibo CU3 changed meaning on the Aer backend** (0.999 instead of 0.696), and
  Qibo circuits with SWAP, SX, CRX, RXX, iSWAP, CCZ and similar gates could not be
  converted for Aer at all. Reported by David Ryan (bsc-wdc/qdislib#13).
- **`find_cut` could return a cut on a Qiskit barrier**, which `wire_cutting` then
  rejected with `KeyError`. Reported by David Ryan (bsc-wdc/qdislib#16).
- **Qiskit circuits with more than one quantum register were cut on the wrong
  qubits.** `circuit_qiskit_to_dag` and `iter_gates` (which `find_cut` reads) used
  each qubit's index within its own register, so a gate between `a[0]` and `b[0]`
  landed on `(0, 0)` and cutting returned a wrong value with no error. They now use
  the global index (`circuit.find_bit(q).index`). Contributed by David Ryan
  (bsc-wdc/qdislib#5, fixes bsc-wdc/qdislib#2).
- **Sampling on `QIBO_BACKEND` changed the caller's circuit.** `sample_circuit`
  added an `M` gate to the Qibo circuit it was given, so every later run of that
  circuit read the wrong result. It now samples a copy. Reported by David Ryan
  (bsc-wdc/qdislib#12).
- **`find_cut(max_cuts=0)` still returned a cut.** The budget checks treated `0` as
  "no limit". `max_cuts=0` now allows no cut in every mode, including `mixed=True`
  (which otherwise still returns an over-budget cut with a warning), and
  `max_qubits < 1` or `max_cuts < 0` raise `ValueError`. Reported by David Ryan
  (bsc-wdc/qdislib#15).
- **`wire_cutting(optimized=True)` was ignored for a weighted sum of Paulis.** Each
  term of the Hamiltonian ran the standard 8**n reconstruction instead of the
  optimized one.
- **A Qiskit subcircuit could only be evaluated once.** `_counts_qiskit` (and the
  GPU/QPU evaluators) appended `measure_all` **in place**, so evaluating the same
  cut twice -- simulate first, then submit to hardware -- silently read the wrong
  classical register and returned a plausible but wrong number. They now measure a
  copy. The IQM path additionally refuses a subcircuit that already carries a
  `meas` register rather than paying for the shots to get a wrong answer.
- **The IQM backend could not authenticate against IQM Resonance.** It required a
  Cortex-CLI tokens *file* (`tokens_path` -> `IQM_TOKENS_FILE`) and hardcoded
  `quantum_computer="default"`, so a plain cloud token raised
  `KeyError: 'tokens_path'`. `token`, `tokens_path`, `url` and `quantum_computer`
  are now all configurable, each falling back to its `IQM_*` env var. The URL is
  the server base URL -- the current client rejects one with the device name in it.
- **`QdBackEnd` assumed a simulator-sized width for hardware.** A `qpu_dict` entry
  without `num_qubits` was recorded as `DEFAULT_NUM_QUBITS` (30), so a 20-qubit
  device would accept 25-qubit circuits. `num_qubits` is now optional everywhere: a
  backend that talks to a device reports the device's own width, and only a
  backend that can resolve nothing falls back to the default.

### Changed
- **Fragment circuits are built where they run.** `circuit_cutting` (and every cut
  type's exact, optimized and `num_samples` runs, which go through it) ships each
  fragment once and a small spec per circuit; the circuit is built in its task, not
  on the master. For an expectation value a circuit's task returns only the few Z
  correlators the reconstruction reads, not all `2**width`.
- **`gate_cutting_hadamard(_sampling)` exact and `num_samples` runs go through
  `circuit_cutting`**, as `gate_cutting`'s already did.
- **`reduce_chunk` is deprecated** (removed in the next release): the
  reconstruction is one cheap task. It still works, for shot runs of the Hadamard
  cut.
- **Breaking: one calling convention for every cut type.** Gate, wire, Hadamard and
  `circuit_cutting` take the same arguments, with the same names and defaults:

      X(circuit, cuts, observable=None, *, backend=None, shots=1024, seed=None, sync=True, ...)
      X_sampling(circuit, cuts, *, backend=None, shots=1024, seed=None, sync=True, ...)
      X_subcircuits(circuit, cuts, observable=None, *, software=None, ...)
      X_sampling_subcircuits(circuit, cuts, *, software=None, ...)

  - `backend` is where the subcircuits run: an SDK name (its simulator), a registered
    backend name or a `build_backend(...)` object; `None` is the circuit's own SDK.
    It replaces `software` on the calls that run circuits, and `qpu`, `qpu_dict`,
    `qpu_config`, `gpu`, `gpu_min_qubits`, `method` and wire cutting's
    `backend="numpy"`. Choosing a device per subcircuit by size (what `qpu_dict` /
    `gpu_min_qubits` did) is done in the three-step loop with `s.num_qubits`.
  - `software` is, on the three-step generators only, the SDK the subcircuits are
    written in (`None`: the circuit's own).
  - `observable` (was `observables`) defaults to all `Z` everywhere (the Hadamard cut
    required it); `cuts` (was `gates_cut` / `cut`); `circuit` (was `dag` /
    `rand_qc`). Only `circuit`, `cuts` and `observable` are positional; a call in the
    old positional order (e.g. `wire_cutting(qc, cut, True, "ZZ")`) raises a
    `TypeError` that says so.
  - `shots=1024` everywhere, `shots=None` is exact everywhere. The Hadamard cut used
    to default to exact.
  - `seed` on every call that samples; each subcircuit gets its own stream.
  - `sync` on every all-in-one call (`gate_cutting_hadamard_sampling` lacked it).
  - `joint=True` everywhere: the Hadamard cut now cuts the gates crossing one
    boundary together by default (`joint=False` for one ancilla per gate).
  - `optimized` on gate cutting too (each side of a cut gate: 5 circuits instead of
    6), and on `wire_cutting_sampling_subcircuits`.
  - Every three-step item is a `Subcircuit`, whatever the cut type: `s.circuit` (the
    native circuit in `software`), `s.software`, `s.num_qubits`, `s.observable`
    (`None` for distribution lists). `expectation_value` and `sample_subcircuit`
    take `shots`, `backend` and `seed`; `s[0]` / `s[1]` still give the circuit and
    the observable, as the 1.x pairs did.
  - `estimate` follows the same convention.

  The 1.x keywords keep working in this release with a `DeprecationWarning` naming
  the replacement, and are removed in the next one. Deprecated the same way:
  `gate_wire_cutting` (use `circuit_cutting`), `cut_cost` (use
  `cutting_overhead`) and the standalone `gate_cutting_subcircuit_reconstruction`,
  `gate_cutting_sampling_reconstruction`, `wire_cutting_subcircuit_reconstruction`
  and `gate_cutting_hadamard_subcircuit_reconstruction` (use
  `subcircuits.reconstruct(results)`).
- **`find_cut(mixed=True)` replaced.** The old mixed mode bisected a qubit-timeline
  graph by node count with 6/8 term weights: it ignored `max_qubits` (a fragment of 6
  qubits came back for `max_qubits=4`) and returned cuts with up to 65,000x the
  `gamma` of a gate-only cut. See "Added" for the cost-aware search that replaces it.
- **SparseCut (`find_cut(coupling_map=...)`) returns a separating cut.** It used to
  return the longest-routing gates whatever they did to the circuit, which the cutting
  functions now refuse when they do not separate. It now returns the cheapest
  separating gate cut made of non-adjacent gates, preferring the longest-routing
  ones, within `max_qubits` when given (`top_n` gives that one cut, scored by its
  total routing distance).
- **Removed `find_cut(implementation="ibm-ckt")`** and the `qiskit-addon-cutting`
  dependency. The IBM cut finder was slow on larger circuits (11.6 s without a proof
  of optimality on a 16-qubit grid, where an exact MILP takes 0.06 s) and Qdislib's
  wrapper dropped the wire cuts it found. Use `mixed=True`.
- **Every cut has to separate its two sides.** All cutting algorithms (gate,
  Hadamard, wire, `gate_wire_cutting`, `circuit_cutting`, with their three-step and
  sampling forms) now check each cut on the circuit with all cuts applied: its two
  sides must land in different fragments. A cut that another, uncut path bypasses
  raises `InvalidCutError` naming it. Before, only the cut set as a whole had to add
  a fragment, so such a cut ran and silently multiplied the circuits and shots for
  nothing. Cuts that separate only in combination (several gates across one
  boundary, a gate cut plus a wire cut) stay valid.
- **Breaking: `_sampling` means the output distribution for every cut type.**
  `gate_cutting_hadamard_sampling` and `gate_cutting_hadamard_sampling_subcircuits`
  used to estimate an expectation value by drawing configurations at random; they now
  return the output distribution, like `gate_cutting_sampling` and
  `wire_cutting_sampling`. The random-configuration estimate is now
  `gate_cutting_hadamard(..., num_samples=N)` (and `gate_cutting_hadamard_subcircuits(...,
  num_samples=N)`). Calling the old way, with an observable, raises a `TypeError` that
  says so.
- **`optimized=True` wire cutting runs fewer measurement circuits by default** (see
  Added: joint measurement). Pass `joint=False` for the previous circuits.
- **`wire_cutting` runs the subcircuits in the circuit's own SDK by default.**
  `software` defaulted to `"qibo"` whatever the input, so a Qiskit circuit
  containing `cx` failed unless `software="qiskit"` was passed. The default is now
  taken from the circuit (Qiskit Aer for a DAG). Reported by David Ryan
  (bsc-wdc/qdislib#14).
- **Backend options are named directly on `BackendConfig`** --
  `BackendConfig(quantum_computer="garnet")` rather than
  `BackendConfig(20, extra={"quantum_computer": "garnet"})`. They are still
  collected into `.extra` (which every backend reads), and `extra={...}` is still
  accepted for callers holding a ready-made dict.
- **The `iqm` extra installs `iqm-client[qiskit]>=35`**, not the older standalone
  `qiskit-iqm` distribution, whose latest release pins `qiskit<1.3` and would
  downgrade the Qiskit stack on install.

## [1.0.1] - 2026-07-20
### Added
- **Tutorial 09** (`09_scaling_up`): scale cutting past H₂. Reconstruct the energy of
  a real **N₂ molecule at 12 qubits** (active space CAS(6,6), 383-term Jordan–Wigner
  Hamiltonian, precomputed in `notebooks/data/n2_cas66_hamiltonian.json`; generation
  cell shown) cut into two 6-qubit halves, and a **25-qubit Ising lattice** cut into
  four blocks with **3 cuts** — both in parallel under PyCOMPSs. Illustrates that the
  cutting cost is `#terms × 6**(#cuts)`, so lattices chase qubit count while molecules
  (O(N⁴) terms) chase accuracy at a modest size.
- **Tutorial 08** (`08_vqe_h2`): estimate H₂'s ground-state energy with VQE by
  cutting a 4-qubit **CUDA-Q** ansatz into two 2-qubit halves. The energy
  `⟨H⟩ = Σ_i c_i ⟨P_i⟩` (15-term Jordan–Wigner Hamiltonian) is rebuilt via the
  three-step workflow, so the ~500 subcircuits (15 terms × 36 cut configurations)
  evaluate **in parallel under PyCOMPSs**; the exact landscape uses `cudaq.observe`.
  Notebooks 07 and 08 are now in the docs learning-path toctree.
- **Parallel Hamiltonian evaluation in `gate_cutting`.** A weighted Pauli-sum
  observable now submits every term's cut without waiting (`sync=False`) and reduces
  with a single COMPSs task, so all `(term × cut-configuration)` subcircuits schedule
  together under PyCOMPSs instead of term-by-term.
- **QPU routing for sampling reconstruction.** `gate_cutting_sampling`,
  `wire_cutting_sampling` and their three-step APIs now accept `qpu` / `qpu_dict`
  (as `gate_cutting` already did), routing each subcircuit's sampling to a QPU. The
  per-subcircuit sampling task resolves its executor through `build_execution_backend`
  (the same seam as the expectation path) instead of an SDK `if`-ladder, via a new
  `BaseBackend.evaluate_joint_counts` hook that captures the joint final/cut-point
  measurements. The local backends (Aer/Qibo/CUDA-Q) implement it; `BSCQpuBackend`
  samples through the shared-filesystem bridge (exact for terminal cuts -- real QPUs
  report only aggregated final counts, so the cut-point is taken to equal the final
  measurement). `_SampleUnit` carries the routing so it survives the three-step API.
- **Tutorial 07** (image reconstruction) additions: a **geometric-transform** step
  (`encode_image(..., transform=...)`) showing that a mirror reflection -- an `X` on
  each address qubit, which stays inside a register and never crosses the cut -- is
  reconstructed faithfully through the cut (scaling is *not* a clean gate: it is a
  resampling, hence non-unitary or qubit-adding); a **"send the pieces to a real BSC
  QPU"** cell (`qpu=True, qpu_dict={"MN_Ona": ...}`); and a rendered comparison of an
  ideal reconstruction versus a noisy (hardware-like) one.
- **Wire-cutting sampling** (`wire_cutting_sampling`, `wire_cutting_sampling_subcircuits`):
  reconstruct the output probability **distribution** of a *wire*-cut (space-like)
  circuit -- the counterpart of `gate_cutting_sampling`. One-pass fragment
  tomography (Perlin et al., arXiv:2005.12702): each fragment is sampled **once**
  (parallel under PyCOMPSs), its Pauli-`Z` correlators are combined via the base-8
  wire-cut decomposition into `<Z^A>` (a measured cut-source qubit is a `Z` factor
  for the X/Y/Z terms and observed with `I` for the identity terms; a cut wire's
  output lives on its fresh qubit line), and `P(x)=2**-N * sum_A (-1)**(A.x) <Z^A>`
  is the inverse Walsh-Hadamard transform, abs-normalised. Full parity with gate
  sampling: all-in-one **and** a three-step API (`wire_cutting_sampling_subcircuits`
  -> `WireSampledCutSubcircuits` -> `sample_subcircuit` -> `reconstruct`) with the
  same create/execute task split. Validated against the exact distribution (TVD at
  shot-noise level) for one and two cuts on Aer/Qibo.
- **Three-step sampling workflow**, mirroring the expectation one:
  `gate_cutting_sampling_subcircuits(circuit, cut, software)` returns a
  `SampledCutSubcircuits` (a `list` of sampling subcircuits + cut metadata); sample
  each with `sample_subcircuit(s, shots)` in your own loop (real QPUs, batching,
  caching, ...), then `subs.reconstruct(sampled)` recombines them into the output
  distribution. The all-in-one `gate_cutting_sampling` and this three-step API share
  one subcircuit-generation path, and the distribution reconstruction is now a task
  fed directly by the sampling futures (parallel under PyCOMPSs, no master barrier).
- **Sampling-based gate cutting** (`gate_cutting_sampling`): reconstruct the full
  output **probability distribution** over bitstrings from a gate-cut circuit,
  not just an expectation value. Each subcircuit is sampled independently
  (embarrassingly parallel -- no classical dependency between subcircuits),
  combined via the quasiprobability decomposition with correct original-qubit
  mapping, then abs-normalised into a valid distribution. Works on Aer / Qibo /
  CUDA-Q via the backend layer; validated against the exact statevector
  distribution (total-variation distance at shot-noise level) for one and two
  cuts.
- **Pluggable backend layer** (`Qdislib.core.backend`): a user-swappable
  execution abstraction so backend choice lives outside the cutting algorithms.
  `QdBackEnd` routes each `sample_circuit` / `get_expectation_value` job to the
  most suitable registered backend (shortest queue first, qubit-count
  compatible, with fallback). Backends implement `BaseBackend` and register via
  `register_backend(...)`; shipped: `AER_BACKEND`, `QIBO_BACKEND`,
  `CUDAQ_BACKEND`, `IQM_BACKEND` (lazy-imported, so only the extras you use are
  required). New `iqm` extra. (Foundation for extracting backend selection out
  of `gate_cutting`/`wire_cutting`.)
- CUDA-Q backend support for gate cutting, wire cutting and `find_cut`
  (`software="cudaq"`), validated against analytical expectation values.
- First-class subcircuit workflow: `expectation_value()` and the
  `CutSubcircuits` object returned by `gate_cutting_subcircuits` /
  `wire_cutting_subcircuits`, with a `.reconstruct()` method that needs no
  manual cut bookkeeping.
- Optional-dependency extras for the quantum SDKs: `qdislib[qibo]`,
  `qdislib[qiskit]`, `qdislib[cudaq]` and `qdislib[all]`.
- Observables can now be a **weighted sum of Pauli operators** (a Hamiltonian),
  not just a single Pauli string: pass a list of `(coefficient, pauli_string)`
  pairs or a `{pauli_string: coefficient}` mapping to `gate_cutting`,
  `wire_cutting` and `analytical_solution`. New `normalize_observables` helper.
- **Semantic circuit cache** (`Qdislib.utils.circuit_cache`): reuse results of
  functionally-equivalent subcircuits, keyed by a ZX-reduced Weisfeiler–Lehman
  hash (`semantic_key`, backend-agnostic across Qibo/Qiskit/CUDA-Q). Pluggable
  stores (`InMemoryCircuitCache`, `LMDBCircuitCache`, `RedisCircuitCache`,
  `RedisPSCOCircuitCache`) via `default_circuit_cache(backend=...)`; `"redis"`
  resolves to the BSC COMPSs PSCO store under PyCOMPSs. Enable with `cache=` on
  `gate_cutting`/`wire_cutting`/`expectation_value`; choose the granularity with
  `cache_mode` (`"expectation"`, `"counts"`, or a custom `CacheStrategy`). Under
  PyCOMPSs a placeholder-reservation + reconstruction-time refresh keeps hit
  rates up under high concurrency. New `cache` extra (`pip install Qdislib[cache]`).
- `set_log_level(...)` (and the `QDISLIB_LOG_LEVEL` environment variable) to turn
  on Qdislib's diagnostic logging with one switch.
- Cross-cluster **QPU bridge** (`Qdislib.utils.qpu_bridge`): `run_qpu_listener()`
  runs on the quantum cluster and executes submitted subcircuits on a BSC QPU
  (Blue or Red) via qililab, communicating with the classical cut job through a
  shared filesystem directory (`subcircuits_dir()` / `$QDISLIB_SUBCIRCUITS_DIR`,
  the single source of truth used by both sides). See `examples/qpu_listener.py`.
- `estimate(circuit, observable, max_qubits=...)`: a one-call helper that
  find_cuts, executes and reconstructs in one step (`return_metrics=True` also
  reports the chosen cut, timings and cache stats).
- `from_qasm(qasm, software=...)`: build a circuit from an OpenQASM 2.0 string,
  so any framework that exports QASM 2.0 (Cirq, PennyLane, Braket, ...) can feed
  Qdislib's cutting workflows.
- `draw_cut(circuit, cut)`: visualise a circuit's DAG with the cut highlighted.
- `CircuitCache.stats()`: hit/miss counts, hit rate and size for a cache.

### Changed
- **Subcircuit creation and execution are now separate PyCOMPSs tasks.** Each
  per-subcircuit fan-out (`_build_and_run_gate_subcircuit`,
  `_build_and_run_subcircuit`, `_sample_gate_subcircuit`) is split into a *create*
  task (prepare the subcircuit at DAG level -> serializable spec) feeding a
  *execute* task (build the native circuit + run/sample it), task -> task. This
  gives finer scheduling granularity -- the execution task is the one to route
  (e.g. with `@constraint` or a QPU backend) to control **where** subcircuits run.
  The native circuit (Qibo collapse circuit / CUDA-Q kernel) is still built and run
  inside a single task, so it never has to cross a task boundary. The three-step
  workflow (`gate_cutting_subcircuits` / `wire_cutting_subcircuits`) likewise now
  creates each subcircuit in its own `_make_subcircuit_task` feeding
  `_evaluate_subcircuit`; under PyCOMPSs those functions therefore return
  `Subcircuit` *futures* (`compss_wait_on` one to inspect it before evaluating).
- **Gate and wire cutting now delegate subcircuit execution to the backend
  layer** instead of branching on `software` internally. The per-backend
  `if software == "qibo"/"qiskit"/"cudaq"` dispatch is replaced by
  `backend.evaluate_expectation(dag, ...)`, and `software=` accepts the name of
  any backend registered via `register_backend` -- so a user can plug in their
  own execution backend without editing the algorithms. GPU is now an Aer
  backend configuration. Reconstruction (including the identity-observable
  sign-decoupling) is unchanged, so results are identical.
- Backend extraction completed for **all** execution targets: the GPU and
  hardware-QPU paths are backends too -- BSC Blue/Red via the shared-filesystem
  bridge (`BSC_QPU_BACKEND`) and IBM Quantum (`IBM_BACKEND`) -- selected by
  `build_execution_backend`, so `gate_cutting`/`wire_cutting` carry no per-device
  branches. Backend estimator observables (`get_expectation_value`) now use a
  consistent qubit-0-first convention across Aer / Qibo / CUDA-Q / IQM.
- `gate_cutting_sampling` runs its subcircuits as PyCOMPSs tasks, so the
  distribution reconstruction is parallel under PyCOMPSs (serial otherwise).
- **IBM Quantum configured explicitly**, not only via environment variables. The
  new `IBMBackend` takes `token` / `channel` / `instance` / `backend_name` /
  `max_time` / `proxies` (or a pre-built `service`) in its config and runs the
  connection + backend selection once in its constructor; `gate_cutting` accepts
  a `qpu_config` dict for the same. Every key falls back to its `IBM_QUANTUM_*`
  env var, so existing setups keep working. The formerly hardcoded BSC proxy
  (`localhost:44433`) is now a config option rather than baked in.
- The quantum SDKs (Qibo / Qiskit / CUDA-Q) are now **optional**: the library
  imports and runs with any single SDK installed, and imports with none.
- `find_cut` is fully backend-agnostic across Qibo, Qiskit and CUDA-Q.
- PyCOMPSs is optional; without it the library runs serially.
- Replaced most `print()` calls with module-level logging.
- **Removed the `verbose=` parameters** from `gate_cutting`, `find_cut` and
  `analytical_solution` (and internal helpers); their diagnostics now log at
  DEBUG level. Use `qd.set_log_level("DEBUG")` instead of `verbose=True`.

### Fixed
- **Semantic cache collided gate-cut subcircuits with opposite cut-point signs.** A
  gate cut's mid-circuit *collapse* measurement carries the quasiprobability
  reconstruction sign, but `semantic_key` hashes only the unitary (QASM → ZX
  `full_reduce`), which drops the collapse -- so two subcircuits with the same unitary
  but opposite sign shared a cache entry, and a wrong cache hit flipped the
  reconstructed expectation (`test_circuit_cache` saw +0.72 vs the exact -0.72). The
  cache key now folds in a fingerprint of the collapse measurements. Wire cutting and
  uncut circuits have none, so their keys are unchanged. (Pre-existing; surfaced once
  the lint fix let CI's test stage run with the `cache` extra.)
- **Three-step gate cutting now honours an identity observable on a non-terminal
  cut qubit** (CPU backends). The subcircuit workflow's guard was called without
  `software`, so it defaulted to "unsupported" and rejected the CPU 3-step even
  though its evaluator recovers the cut-point sign (prep qubits auto-detected from
  the subcircuit's mid-circuit measurements) exactly as the all-in-one path does.
  Passing `software` lets the guard fire only for the paths that genuinely can't
  (GPU / QPU). This is what lets the H₂ VQE tutorial (notebook 08) reconstruct
  `⟨H⟩` with the three-step workflow.
- **Three-step reconstruction no longer forces a master barrier.** The 3-step
  reconstructors used to `compss_wait_on` every subcircuit result on the master
  before applying the per-configuration sign; now each configuration's subcircuit
  futures feed straight into `_change_sign` / `_change_sign_gate_cutting`
  (task -> task), with the semantic-cache placeholder resolve moved into the sign
  task. Under PyCOMPSs the `_evaluate_subcircuit -> _change_sign` dependency is now
  direct (no sync node between them; confirmed in the `-g` task graph).
- **No nested PyCOMPSs tasks in the cutting paths.** The backend-abstraction
  refactor left the per-subcircuit leaf helpers (`_counts_qiskit`, `_counts_qibo`,
  `_expec_value_cudaq`, and the QPU/GPU `_expec_value_*`) decorated with `@task`
  even though they only ever run *inside* a per-subcircuit `@task`
  (`_build_and_run_*`, `_sample_gate_subcircuit`, `_evaluate_subcircuit`) -- a
  `@task` calling a `@task`, which under PyCOMPSs can hand the worker an
  unresolved future. Those leaves are now plain functions; every remaining
  `@task` is either a per-subcircuit fan-out or a master-side reduction over
  futures (`_change_sign*`, `_sum_results`, `COLLECTION_IN`). Verified under
  `runcompss` (expectation, and qiskit + CUDA-Q sampling, all reconstruct
  correctly with a flat fan-out and no task-to-task edges in the `-g` graph).
- **CUDA-Q sampling/expectation under `cudaq` >= 0.12.** `cudaq.run` no longer
  supports returning a `list` from entry-point kernels, so the gate-cut record
  kernel (`_build_cudaq_run_kernel`) now returns the per-shot outcomes packed into
  an int (bit `k` = outcome `k` of `[cut-points..., finals...]`), unpacked on the
  host with `(record >> k) & 1`. Reconstruction is unchanged.
- **Gate-cut reconstruction corrected for cut qubits that keep evolving after the
  cut.** A gate cut's measure term is the projective `(I ± Z)/2` (Mitarai & Fujii,
  arXiv:1909.07534): the cut qubit must be *measured and left collapsed to its
  outcome*, then continue, contributing the decomposition sign `-(-1)**outcome`
  for **any** observable. The converters instead **forced** the qubit to a fixed
  state after measuring (`|0>` on Qiskit, `|1>` on Qibo/CUDA-Q), which is harmless
  only when nothing follows the cut on that qubit -- but once a cut qubit had
  further gates (e.g. two cuts sharing an internally-entangled block), both the
  expectation value and the sampled distribution came out systematically wrong (a
  cross-cut `<Z Z>` could flip sign; sampling TVD ~0.13-0.22, independent of
  shots). The converters now **collapse-and-continue** (no reset) and the
  reconstruction applies the cut-point decomposition sign for every observable
  (`Z` multiplies it by the final measurement, identity is just the sign). Fixed
  for the expectation and `gate_cutting_sampling` paths, the all-in-one and 3-step
  `gate_cutting_subcircuits` workflows, and the GPU path, across **Aer, Qibo and
  CUDA-Q** -- validated against the exact statevector for terminal *and*
  non-terminal cuts, one and two cuts, product and entangled blocks (TVD
  ~0.001-0.02). A subcircuit that is a single bare qubit (touched only by a removed
  cut node) no longer crashes the sampler.
- **CUDA-Q gate cutting now executes via `cudaq.run`**, not `cudaq.sample`. The
  sampler does not collapse mid-circuit measurements (and rejects conditional
  feedback), so the measure-and-prepare gate-cut path was silently wrong whenever a
  cut qubit kept evolving. The subcircuit is now emitted as a `@cudaq.kernel` that
  returns its per-shot measurement record and is run with `cudaq.run` (which
  collapses), fixing the CUDA-Q backend for non-terminal and multi-cut/entangled
  cases to match Aer/Qibo.
- **Sampling gate cutting** (`gate_cutting_sampling`) reconstructs the full output
  distribution from the (now-correct) Pauli-Z correlators as
  `P(x) = 2**-n * sum_A (-1)**(A·x) <Z^A>` (a Walsh-Hadamard transform), so it
  inherits the expectation path's correctness rather than needing a separate
  derivation.
- **CX/CNOT gate cutting.** `gate_cutting` can now cut a `CX`/`CNOT` directly
  (previously only `CZ` was sound; cutting a bare `CNOT` silently returned a wrong
  value). A `CX` is cut in place via the identity `CX(c,t) = H(t)·CZ(c,t)·H(t)` --
  `find_cut` still names the cut after the gate the user wrote (no phantom `CZ`),
  and the Hadamard conjugation is applied around the target-qubit measure-and-
  prepare at the cut site. Cutting any other two-qubit gate now raises
  `InvalidCutError` instead of mis-cutting it as a `CZ`.
- Identity (`"I"`) observables in **gate cutting** now reconstruct correctly for
  any qubit on all three CPU backends. Previously an `I` on a qubit whose last
  gate is multi-qubit was silently dropped (turned into `Z`) in the DAG
  converters, and an `I` on a *cut* qubit that feeds a multi-qubit gate after the
  cut produced a sign-flipped result. The marker is now attached to the last gate
  acting on the qubit (multi-qubit terminals included), and the gate-cut
  reconstruction recovers the measure-and-prepare sign from the **cut-point
  measurement** (decoupled from the final observable) on `qiskit` (cut-measure
  classical register), `qibo` (collapse gate per-shot samples) and `cudaq`
  (`sample(explicit_measurements=True)`). The all-in-one, 3-step and GPU paths now
  all use this decoupled reconstruction (the 3-step and GPU auto-detect the cut
  qubits from the subcircuit); only the hardware-QPU path -- which sees aggregated
  counts and so cannot recover the per-shot cut-point correlation -- still raises
  `InvalidCutError` rather than returning a wrong value. Wire cutting was
  unaffected by the sign issue and also benefits from the converter fix.
- QPU hand-off no longer hangs forever: the wait for a result is now bounded
  (`QPUTimeoutError` after `$QDISLIB_QPU_TIMEOUT`, default 1 h) and files are
  written atomically, so a down listener or a half-written file no longer stalls
  or corrupts a run. Measurements are added once (on the listener side).
- Identity (`"I"`) observables on the Qiskit backend: `get_counts()` bitstrings
  are little-endian and multi-register, so they are now reduced to the
  measurement register and reversed before being paired with the observable.
- `X`/`Y` observables now work for any input backend: measurement-basis rotation
  is applied on the backend-neutral DAG (`change_basis_dag`) instead of the input
  circuit, which only worked for Qibo inputs.
- Y-basis change-of-basis in wire cutting; gate-cut sign convention for the
  Qibo and CUDA-Q backends.
- `find_cut` no longer crashes on small circuits (Kernighan-Lin / Girvan-Newman
  on graphs too small to partition).
- Multi-edge graph construction in `find_cut` that could return a wrong cut.
- `find_cut` Kernighan-Lin search no longer aborts after the first
  constraint-violating run.
- `gate_cutting` no longer raises `TypeError` when `qpu=True` and
  `qpu_dict=None`.

### Removed
- Dead/commented-out `qubit_mapping` module and unused converter duplicates.

## [1.0.0] - 2025-05-30
### Added
- Added support for sequential exectution

### Changed
- Nothing

### Fixed
- Fixed several bugs

### Removed
- Nothing



## [0.0.1] - 2024-03-01
### Added
- Initial code

### Changed
- Nothing

### Fixed
- Nothing

### Removed
- Nothing
