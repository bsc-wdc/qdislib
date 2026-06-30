# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
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

### Changed
- The quantum SDKs (Qibo / Qiskit / CUDA-Q) are now **optional**: the library
  imports and runs with any single SDK installed, and imports with none.
- `find_cut` is fully backend-agnostic across Qibo, Qiskit and CUDA-Q.
- PyCOMPSs is optional; without it the library runs serially.
- Replaced most `print()` calls with module-level logging.

### Fixed
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
