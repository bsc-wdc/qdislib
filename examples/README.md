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
