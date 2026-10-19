
## Quickstart guide

There are two ways in which you can get started with Qdislib. You can perform
a manual installation, or you can download our ready-to-use docker image.

### Manual installation

#### Dependencies

Qdislib's core dependencies (installed automatically with `pip install qdislib`)
are: `numpy`, `networkx`, `matplotlib`, `pymetis` and `scikit-learn`.

In addition you need **at least one quantum SDK backend**, installed as an
optional extra:

* `pip install "qdislib[qibo]"`   — Qibo
* `pip install "qdislib[qiskit]"` — Qiskit (+ Aer, qiskit-addon-cutting, IBM runtime)
* `pip install "qdislib[cudaq]"`  — CUDA-Q
* `pip install "qdislib[pennylane]"` — PennyLane
* `pip install "qdislib[all]"`    — all of the above, plus the semantic cache (`[cache]`)
* `pip install "qdislib[iqm]"`    — IQM hardware (not part of `[all]`)

For distributed execution on HPC you also need **PyCOMPSs >= 3.3** (see the
installation steps below). GPU acceleration may additionally require
`qiskit-aer-gpu` (with Qiskit) or `qibojit` (with Qibo).

**Troubleshooting:** on Linux with glibc older than 2.32 (e.g. Debian 11),
scikit-learn's bundled OpenMP runtime can fail to load with `cannot allocate
memory in static TLS block`. `find_cut` then skips its spectral strategy with a
warning. To keep it, preload the library before starting Python:

```bash
export LD_PRELOAD="$(python -c 'import glob, site; print(":".join(glob.glob(site.getsitepackages()[0] + "/scikit_learn.libs/libgomp*.so*")))')"
```

#### Installation steps

1. Check which PyCOMPSs version to install.
    * Latest Qdislib release requires **PyCOMPSs 3.3** or greater (check [Qdislib releases](https://github.com/bsc-wdc/Qdislib/releases) for information about other releases).

2. Install PyCOMPSs following these [instructions](https://compss-doc.readthedocs.io/en/stable/Sections/01_Installation.html).

3. Install the latest Qdislib version together with a backend extra, e.g.
   ``pip3 install "Qdislib[qibo]"`` (see Dependencies above for the options).
   * **IMPORTANT:** Qdislib uses the ``pycompss`` Python module for distributed
     execution, but it is **not** installed automatically by pip. It becomes
     available after manually installing PyCOMPSs (step 2). Without PyCOMPSs,
     Qdislib still runs correctly in serial mode.

4. Check that everything works by running the three-step circuit-cutting
   workflow from the [Quickstart](https://qdislib.readthedocs.io/en/latest/quickstart.html):

    ```python
    import Qdislib as qd
    from qibo import models, gates

    circuit = models.Circuit(2)
    circuit.add(gates.H(0)); circuit.add(gates.CZ(0, 1)); circuit.add(gates.RY(0, 0.8))

    subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]
    print(subcircuits.reconstruct(values))  # ≈ -0.717
    ```

   To run distributed on PyCOMPSs, launch your script with ``runcompss your_app.py``.
