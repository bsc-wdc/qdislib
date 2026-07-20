# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
