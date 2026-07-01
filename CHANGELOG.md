# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
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
- The quantum SDKs (Qibo / Qiskit / CUDA-Q) are now **optional**: the library
  imports and runs with any single SDK installed, and imports with none.
- `find_cut` is fully backend-agnostic across Qibo, Qiskit and CUDA-Q.
- PyCOMPSs is optional; without it the library runs serially.
- Replaced most `print()` calls with module-level logging.
- **Removed the `verbose=` parameters** from `gate_cutting`, `find_cut` and
  `analytical_solution` (and internal helpers); their diagnostics now log at
  DEBUG level. Use `qd.set_log_level("DEBUG")` instead of `verbose=True`.

### Fixed
- Identity (`"I"`) observables in **gate cutting** now reconstruct correctly for
  any qubit on all three CPU backends. Previously an `I` on a qubit whose last
  gate is multi-qubit was silently dropped (turned into `Z`) in the DAG
  converters, and an `I` on a *cut* qubit that feeds a multi-qubit gate after the
  cut produced a sign-flipped result. The marker is now attached to the last gate
  acting on the qubit (multi-qubit terminals included), and the gate-cut
  reconstruction recovers the measure-and-prepare sign from the **cut-point
  measurement** (decoupled from the final observable) on `qiskit` (cut-measure
  classical register), `qibo` (collapse gate per-shot samples) and `cudaq`
  (`sample(explicit_measurements=True)`). The GPU/QPU paths and the 3-step
  `gate_cutting_subcircuits` workflow do not decouple this yet and raise
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
