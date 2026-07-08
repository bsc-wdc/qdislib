<h1 align="center">
    <img src="https://github.com/bsc-wdc/qdislib/raw/master/docs/logos/dislib-logo-full.png" alt="Qdislib - The Quantum Distributed Computing Library" height="90px">
</h1>

<h3 align="center">Quantum distributed computing library implemented over PyCOMPSs programming model for HPC.</h3>

<p align="center">
  <a href="https://qdislib.readthedocs.io/?badge=latest">
    <img src="https://readthedocs.org/projects/qdislib/badge/?version=latest" alt="Documentation Status"/>
  </a>
  <a href="https://badge.fury.io/py/qdislib">
      <img src="https://badge.fury.io/py/qdislib.svg" alt="PyPI version" height="18">
  </a>
  <a href="https://badge.fury.io/py/qdislib">
      <img src="https://img.shields.io/badge/python-3.10-blue.svg" alt="Python version" height="18">
  </a>
</p>

<!-- <p align="center"><b>
    <a href="https://qdislib.bsc.es">Website</a> •
    <a href="https://qdislib.bsc.es/en/stable/api-reference.html">Documentation</a> •
    <a href="https://github.com/bsc-wdc/qdislib/releases">Releases</a> •
    <a href="https://bit.ly/bsc-wdc-community">Slack</a>
</b></p> -->

 **Table of content:**
- [Description](#description)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Availability](#availability)
- [URLs](#urls)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Description

**Qdislib** is a Python library designed for scalable quantum circuit execution using *circuit cutting* techniques. It enables the simulation of large quantum circuits by splitting them into smaller, manageable subcircuits that can be executed independently—either on classical simulators, GPUs, or quantum hardware.

Qdislib is built on top of the [PyCOMPSs](<https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar>) parallel runtime, allowing seamless distributed execution of quantum workloads across CPUs, GPUs, and QPUs.

With Qdislib, researchers and developers can:

- Perform **gate** and **wire cutting** to decompose complex quantum circuits — including an ancilla-based **Hadamard-test gate cut** (`gate_cutting_hadamard`) that cuts *any* two-qubit gate at optimal overhead with no mid-circuit measurement, and **mixed** wire + gate cutting on one circuit (`gate_wire_cutting`).
- Cut **any two-qubit gate** (CZ/CX natively; others auto-decomposed to CZ, or use the Hadamard cut).
- Leverage **GPU acceleration** using cuQuantum or Qibojit.
- Submit subcircuits to **remote QPUs** like IBM Quantum.
- Work with circuits defined in **Qibo**, **Qiskit**, **CUDA-Q**, and **PennyLane**.
- Automatically identify good cut points with `find_cut`.
- Extract and manipulate subcircuits independently.

Whether you're targeting HPC systems, hybrid quantum-classical setups, or constrained simulators, Qdislib is a flexible and modular tool to bridge the gap between current hardware limitations and large-scale quantum algorithm design.

Explore the sections below to get started with installation, quickstart examples, user guides, API references, and more.


Qdislib has been implemented on top of [PyCOMPSs](<https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>) programming model,
and it is being developed by the [Workflows and Distributed Computing](<https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>) group of the [Barcelona Supercomputing Center](<http://www.bsc.es>).



## Requirements

- Python >= 3.10
- At least one quantum SDK backend — **Qibo**, **Qiskit**, or **CUDA-Q** — installed as an extra (see below)
- *(optional)* COMPSs >= 3.3 for distributed execution on HPC; without it, Qdislib runs serially

The core Python dependencies are listed in `requirements.txt`; the quantum SDK
backends are optional extras declared in `pyproject.toml`.

## Installation

Install Qdislib together with at least one backend extra:

```bash
pip3 install "qdislib[qibo]"      # or [qiskit], [cudaq], or [all]
```

Installing the bare package (`pip3 install qdislib`) pulls only the core
dependencies — you need a backend extra to run circuits.

## Quickstart

The recommended way to use Qdislib is the **three-step circuit-cutting workflow**:
generate subcircuits, evaluate each one in your own loop, then recombine the
results. Evaluating the subcircuits yourself is the whole point — that loop is
where you dispatch to a simulator, a GPU, or a real QPU, batch them, or cache
results, all while Qdislib handles the cutting and the reconstruction.

```python
import Qdislib as qd
from qibo import models, gates

# A small example circuit (works the same with Qiskit or CUDA-Q circuits).
circuit = models.Circuit(2)
circuit.add(gates.H(0))
circuit.add(gates.CZ(0, 1))
circuit.add(gates.RY(0, 0.8))

# 1. Cut the circuit into independently executable subcircuits.
#    Use find_cut to discover good cut points automatically.
subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")

# 2. Evaluate each subcircuit — this loop is yours to control.
values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

# 3. Recombine into the original circuit's expectation value.
result = subcircuits.reconstruct(values)
print(result)  # ≈ -0.717
```

`qd.wire_cutting_subcircuits(circuit, [("H_1", "CZ_2")], software="qibo")` works
the same way for **wire cutting**. The `software` argument (`"qibo"`, `"qiskit"`
or `"cudaq"`) selects the backend the subcircuits are emitted for; it must match
the backend you evaluate them on.

If you'd rather not manage the loop, the all-in-one `qd.gate_cutting(...)` and
`qd.wire_cutting(...)` run the cut, execution and reconstruction in a single call
— convenient for quick experiments, but the three-step workflow above is the
flexible, production-oriented path.

## Availability

Currently, the following supercomputers have already PyCOMPSs and Qdislib installed and ready to use. If you need help configuring your own cluster or supercomputer, drop us an email and we will be pleased to help.

- Marenostrum 5 - Barcelona Supercomputing Center (BSC)


## Citing Qdislib

If you use Qdislib in a scientific publication, we would appreciate citations to
the following papers. Please cite the original paper (ACM, SC '25 Workshops)
first:

M. Tejedor, B. Casas, J. Conejero, A. Cervera-Lierta and R. M. Badia,
"Orchestrating Quantum-HPC Workflows with Distributed Quantum Circuit Cutting,"
in *Proceedings of the SC '25 Workshops of the International Conference for High
Performance Computing, Networking, Storage and Analysis*, ACM, 2025, pp.
1898–1906. https://doi.org/10.1145/3731599.3767547

M. Tejedor, J. Conejero and R. M. Badia, "A Semantic Quantum Circuit Cache for
Scalable and Distributed Quantum-Classical Workflows," arXiv:2604.26788, 2026.
https://arxiv.org/abs/2604.26788

### Bibtex:

```latex
@inproceedings{Qdislib,
    title     = {{Orchestrating Quantum-HPC Workflows with Distributed Quantum Circuit Cutting}},
    author    = {Mar Tejedor and Berta Casas and Javier Conejero and Alba Cervera-Lierta and Rosa M. Badia},
    booktitle = {Proceedings of the SC '25 Workshops of the International Conference for High Performance Computing, Networking, Storage and Analysis},
    publisher = {ACM},
    pages     = {1898--1906},
    year      = {2025},
    doi       = {10.1145/3731599.3767547},
}

@article{QdislibCache,
    title   = {{A Semantic Quantum Circuit Cache for Scalable and Distributed Quantum-Classical Workflows}},
    author  = {Mar Tejedor and Javier Conejero and Rosa M. Badia},
    journal = {arXiv preprint arXiv:2604.26788},
    year    = {2026},
}
```

## URLs

[Quantum Distributed computing libraries BSC](https://www.bsc.es/research-development/research-areas/distributed-computing/distributed-computing-libraries-and)


## Acknowledgements

The project acknowledges funding from the Spanish Ministry for Digital Transformation and of Civil Service of the Spanish Government through the QUANTUM ENIA project call - Quantum Spain, EU through the Recovery, Transformation and Resilience Plan – NextGenerationEU within the framework of the Digital Spain 2026. It acknowledges funding from Grant RYC2022-037769-I funded by MICIU/AEI/10.13039/501100011033 and by “ESF+”. It also acknowledges funding from projects  CEX2021-001148-S, and PID2023-147979NB-C21 from the  MCIN/AEI and MICIU/AEI /10.13039/501100011033 and by FEDER, UE, by the Departament de Recerca i Universitats de la Generalitat de Catalunya, research group MPiEDist (2021 SGR 00412).

## License

Apache License Version 2.0, see [LICENSE](LICENSE)
