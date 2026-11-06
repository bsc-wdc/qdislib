# Qdislib examples

Scripts for Qdislib's distributed path. The [notebooks](../notebooks) cover the
API; these run the same cutting functions as COMPSs tasks, against the exact answer,
so they are usage examples and checks at once. They also run serially: Qdislib's tasks
then run in the same process.

Run them from the repository root, with Qdislib and its SDKs installed
(`pip install ".[qibo,qiskit]"`). Under COMPSs (`runcompss`, COMPSs 3.3 or later) add
`-g` to write the task graph under `$HOME/.COMPSs/<app>/monitor/`.

| Script | What it does | Prints |
|---|---|---|
| [`compss_app.py`](compss_app.py) | Gate, wire and Hadamard cuts, all-in-one and three-step, on 5 qubits | `RESULT: PASS` |
| [`mixed_cutting_compss.py`](mixed_cutting_compss.py) | A gate cut and a wire cut on one circuit | `RESULT: PASS` |
| [`hadamard_cut_compss.py`](hadamard_cut_compss.py) | The Hadamard cut: every configuration, drawn ones and three-step, overlapping with `sync=False` | `RESULT: PASS` |
| [`verify_sampling_pycompss.py`](verify_sampling_pycompss.py) | Distributions against the exact one, on each SDK in `QDISLIB_SOFTWARE` (default `qiskit`) | `RESULT: PASS` |
| [`verify_engine_pycompss.py`](verify_engine_pycompss.py) | Backend objects, circuit groups, marginals, error bars, free parameters and a 26-qubit distribution (needs Qiskit) | `RESULT: PASS` |
| [`verify_gpu_pycompss.py`](verify_gpu_pycompss.py) | GPU tasks: a GPU backend gets one GPU per task, a CPU one none (needs GPU nodes and CUDA-Q). `_gpu_probe.py` holds its backends | `RESULT: PASS` |
| [`qpu_listener.py`](qpu_listener.py) | The listener for BSC QPUs, on the quantum cluster (needs `qdislib[qpu]`) | |

```bash
python3 examples/compss_app.py                                                # serial
runcompss --pythonpath=$(pwd) --python_interpreter=python3 examples/compss_app.py   # with COMPSs

QDISLIB_SOFTWARE=qiskit python3 examples/compss_app.py                        # another backend
```

## In a container

[`Dockerfile.compss`](../Dockerfile.compss) builds on the COMPSs tutorial image and runs
the examples and the unit tests under `runcompss`, through
[`scripts/run_compss_tests.sh`](../scripts/run_compss_tests.sh):

```bash
docker build -f Dockerfile.compss -t qdislib-compss .   # QDISLIB_EXTRAS=all also installs CUDA-Q
docker run --rm qdislib-compss
```

`runcompss` exits 0 even when an inner test fails, so the script checks the output for
`RESULT: PASS` and the unittest `OK` line instead.

## On MareNostrum 5

[`scripts/mn5/check.sh`](../scripts/mn5/check.sh) submits `verify_engine_pycompss.py`
(`engine`, on a GPP login) or `verify_gpu_pycompss.py` (`gpu`, on an ACC login):

```bash
scripts/mn5/check.sh engine <account> <qos>
```

Each log ends with `RESULT: PASS` or `RESULT: FAIL ...`: grep for it in `compss-<job>.out`.

## The task graphs here

`hadamard_compss_graph.*` and `mixed_compss_graph.*` are the task graphs of
`hadamard_cut_compss.py` and `mixed_cutting_compss.py` from a real run
(`scripts/run_hadamard_graph.sh`, `scripts/run_mixed_graph.sh`, in the container).

## Running subcircuits on a BSC QPU

The classical (PyCOMPSs) cut job and the quantum-control job run on separate clusters
that cannot share a scheduler allocation, so they communicate through a shared directory,
`$QDISLIB_SUBCIRCUITS_DIR`. The cut job writes each subcircuit there and waits for its
result; a listener on the quantum cluster runs it on a QPU, Blue or Red, and writes the
result back.

On the quantum cluster:

```bash
export QDISLIB_SUBCIRCUITS_DIR=/gpfs/projects/.../subcircuits   # shared by both clusters
export RUNCARD_RED=/path/to/red_runcard.yml                     # or RUNCARD_BLUE
python examples/qpu_listener.py --qpu red --nshots 1000
```

On the classical cluster, with the same directory:

```python
qd.gate_cutting(circuit, cut, backend="BSC_QPU_BACKEND")
```

The listener runs until a `STOP` file appears in the directory, Ctrl-C, or
`--idle-timeout`. If it is down, submissions fail with a `QPUTimeoutError` after
`$QDISLIB_QPU_TIMEOUT` (default 1 h) instead of hanging.
