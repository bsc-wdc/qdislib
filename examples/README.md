# Qdislib examples

Runnable examples for Qdislib's distributed (PyCOMPSs / `runcompss`) execution
path. The Jupyter notebooks under [`../notebooks/`](../notebooks) cover the
serial API; this directory covers running the same circuit-cutting workflows
**distributed across COMPSs tasks**.

## `compss_app.py`

A self-contained smoke test that builds a small 5-qubit Qibo circuit, finds a
cut, and reconstructs `<ZZZZZ>` four ways:

- all-in-one `gate_cutting` / `wire_cutting`, and
- the recommended three-step workflow
  (`gate_cutting_subcircuits` → `expectation_value` → `reconstruct`).

It compares every result against `analytical_solution` and prints
`RESULT: PASS` / `RESULT: FAIL` (exiting non-zero on failure), so it works as
both a usage example and a CI check of the parallel path.

### Serial (no PyCOMPSs)

```bash
python3 examples/compss_app.py
```

### Distributed (PyCOMPSs)

```bash
runcompss --pythonpath=$(pwd) --python_interpreter=python3 examples/compss_app.py
```

Select the backend with the `QDISLIB_SOFTWARE` environment variable
(`qibo` — the default, `qiskit`, or `cudaq`):

```bash
QDISLIB_SOFTWARE=qiskit runcompss --pythonpath=$(pwd) --python_interpreter=python3 examples/compss_app.py
```

## Serial vs. distributed: what changes

The public API is identical in both modes — the difference is internal:

- **Serial** (PyCOMPSs not installed): the `@task` / `compss_wait_on` decorators
  in [`Qdislib/utils/runtime.py`](../Qdislib/utils/runtime.py) degrade to
  no-ops, and `gate_cutting_subcircuits` / `wire_cutting_subcircuits` hand back
  **native circuits** (Qibo/Qiskit objects, or a CUDA-Q `[kernel, observables]`
  pair).
- **Distributed** (under `runcompss`): each subcircuit evaluation runs as its
  own COMPSs `@task`. Because native Qibo subcircuits and CUDA-Q kernels cannot
  be pickled to remote workers, `*_subcircuits` instead returns picklable
  `Subcircuit` wrappers (holding the circuit DAG); each worker materialises the
  native circuit locally before executing it. Qiskit circuits pickle directly
  but use the same wrapper for uniformity.

So you write the three-step workflow once and it runs serially or distributed
unchanged.

## Running everything under COMPSs in Docker

[`Dockerfile.compss`](../Dockerfile.compss) builds on the official COMPSs
tutorial image and runs both the example app and the full unit-test suite under
`runcompss` (via [`scripts/run_compss_tests.sh`](../scripts/run_compss_tests.sh)):

```bash
# qibo + qiskit (default); pass QDISLIB_EXTRAS=all to also install cudaq
docker build -f Dockerfile.compss -t qdislib-compss .
docker run --rm qdislib-compss
```

`runcompss` exits 0 even when an inner unittest fails, so the script inspects
the captured output for `RESULT: PASS` and the unittest `OK` summary line
instead of trusting the exit code.

## Semantic circuit cache

Cutting produces many *functionally equivalent* subcircuits (repeated fragments,
or the same ansatz across variational iterations). The semantic cache
([`Qdislib/utils/circuit_cache.py`](../Qdislib/utils/circuit_cache.py)) keys each
subcircuit by a ZX-reduced Weisfeiler–Lehman hash, so an equivalent fragment is
simulated once and reused. It never changes the result — an uncacheable circuit
just falls back to a normal simulation (requires the `cache` extra:
`pip install Qdislib[cache]`).

You pick the store; `"redis"` becomes the BSC COMPSs PSCO store under PyCOMPSs
(shared across tasks/nodes) and a plain redis-py store otherwise:

```python
import Qdislib as qd

cache = qd.default_circuit_cache("lmdb", path="/scratch/qcache")   # or "redis" / "memory"

# one-call API:
value = qd.gate_cutting(circuit, cut, software="qibo", cache=cache)

# or the three-step workflow, switching granularity / a custom payload:
subs   = qd.wire_cutting_subcircuits(circuit, cut, software="qibo")
values = [qd.expectation_value(s, cache=cache, cache_mode="counts") for s in subs]
result = subs.reconstruct(values)
```

Under PyCOMPSs, a subcircuit whose key another task is already computing is
deferred (a placeholder is reserved) and refreshed at reconstruction — so high
concurrency doesn't turn would-be hits into duplicate simulations.

## Diagnostic logging

The old per-function `verbose=` flags are gone; enable Qdislib's diagnostics
(candidate cut scores, chosen backend, reconstruction values, ...) with a single
switch, or the `QDISLIB_LOG_LEVEL` environment variable:

```python
import Qdislib as qd
qd.set_log_level("DEBUG")   # or "INFO"; default is quiet
```
