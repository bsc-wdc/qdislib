## Quickstart guide

### Dependencies

Qdislib needs Python 3.10 or later. Its core dependencies, installed with it, are `numpy`,
`networkx`, `matplotlib`, `scikit-learn` and, where it has a wheel, `pymetis`.

You also need **at least one quantum SDK**, installed as an optional extra. Use `pip` or
[uv](https://docs.astral.sh/uv/), which installs the same package and is typically much faster:
`uv pip install "qdislib[qiskit]"` in an environment, `uv add "qdislib[qiskit]"` in a uv project,
`uv run --with "qdislib[qiskit]" app.py` for a one-off.

* `pip install "qdislib[qiskit]"` : Qiskit, with Aer and the IBM runtime
* `pip install "qdislib[qibo]"` : Qibo
* `pip install "qdislib[cudaq]"` : CUDA-Q
* `pip install "qdislib[pennylane]"` : PennyLane
* `pip install "qdislib[qilisdk]"` : QiliSDK
* `pip install "qdislib[cache]"` : PyZX and LMDB, for the semantic cache keys and the on-disk cache
* `pip install "qdislib[all]"` : all of the above: about 3.3 GB, mostly NVIDIA libraries for CUDA-Q, so most people do not want it
* `pip install "qdislib[iqm]"` : IQM hardware (not part of `[all]`)
* `pip install "qdislib[qpu]"` : the BSC QPU listener, on the quantum cluster (not part of `[all]`)

A bare `pip install qdislib` has no SDK, so it can import Qdislib but cannot run a circuit.

**Linux arm64 and Windows:** `pymetis` has no wheel there, so it is not installed. `find_cut` then
warns and runs without METIS, and a wire cut can come out larger. `pip install "qdislib[metis]"`
adds it, and needs a C compiler. Windows has been checked for wheel availability only, not run.

A GPU may need `qiskit-aer-gpu` (with Qiskit), or a CUDA-Q with its `nvidia` target. Cutting
some gates with the Hadamard test also needs Qibo, for a circuit in any SDK.

**Troubleshooting:** on Linux with glibc older than 2.32 (e.g. Debian 11), the OpenMP runtime
that scikit-learn and qiskit-aer each bundle (`libgomp`) can fail to load with `cannot allocate
memory in static TLS block`, the more likely the more SDKs are imported. (`find_cut` then skips
its spectral strategy with a warning; Aer cannot run at all.) Preload every bundled copy before
starting Python:

```bash
export LD_PRELOAD="$(python -c 'import os, site; sp = site.getsitepackages()[0]; print(":".join(sorted(os.path.join(d, f) for d, _, fs in os.walk(sp) for f in fs if f.startswith("libgomp") and ".so" in f)))')"
```

### Installation

1. Install Qdislib with a backend extra: `pip install "qdislib[qiskit]"` or `uv pip install "qdislib[qiskit]"`.

2. Check that it works:

    ```python
    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(2)
    circuit.h(0)
    circuit.cz(0, 1)
    circuit.ry(0.8, 0)

    subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"])
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]
    print(subcircuits.reconstruct(values))          # about -0.717
    print(qd.analytical_solution(circuit, "ZZ"))    # -0.717, exact
    ```

### To run in parallel

Without PyCOMPSs, Qdislib runs serially, with the same results. To run in parallel:

1. Install **PyCOMPSs 3.3** or later following
   [its instructions](https://compss-doc.readthedocs.io/en/stable/Sections/01_Installation_Configuration.html).
   `pip` does not install the `pycompss` module with Qdislib.

2. Launch your script with `runcompss --pythonpath=$(pwd) your_app.py`.

`examples/compss_app.py` is a script to try this with, and it checks its results. To try COMPSs
without installing it, `Dockerfile.compss` builds an image with it (see
[`examples/README.md`](examples/README.md)). On MareNostrum 5, see
[the HPC guide](https://qdislib.readthedocs.io/en/latest/guide/hpc.html).
