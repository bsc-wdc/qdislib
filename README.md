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
- [Availability](#availability)
- [URLs](#urls)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Description

**Qdislib** is a Python library designed for scalable quantum circuit execution using *circuit cutting* techniques. It enables the simulation of large quantum circuits by splitting them into smaller, manageable subcircuits that can be executed independently—either on classical simulators, GPUs, or quantum hardware.

Qdislib is built on top of the [PyCOMPSs](<https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar>) parallel runtime, allowing seamless distributed execution of quantum workloads across CPUs, GPUs, and QPUs.

With Qdislib, researchers and developers can:

- Perform **gate** and **wire cutting** to decompose complex quantum circuits.
- Leverage **GPU acceleration** using cuQuantum or Qibojit.
- Submit subcircuits to **remote QPUs** like IBM Quantum.
- Work with circuits defined in both **Qibo** and **Qiskit**.
- Automatically identify good cut points with `find_cut`.
- Extract and manipulate subcircuits independently.

Whether you're targeting HPC systems, hybrid quantum-classical setups, or constrained simulators, Qdislib is a flexible and modular tool to bridge the gap between current hardware limitations and large-scale quantum algorithm design.

Explore the sections below to get started with installation, quickstart examples, user guides, API references, and more.


Qdislib has been implemented on top of [PyCOMPSs](<https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>) programming model,
and it is being developed by the [Workflows and Distributed Computing](<https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>) group of the [Barcelona Supercomputing Center](<http://www.bsc.es>).



## Requirements

- Python >= 3.10
- COMPSs >= 3.3

Python packages required are defined in `requirements.txt`:

## Installation

Qdislib can be installed with the following command:

```bash
pip3 install qdislib
```

## Availability

Currently, the following supercomputers have already PyCOMPSs and Qdislib installed and ready to use. If you need help configuring your own cluster or supercomputer, drop us an email and we will be pleased to help.

- Marenostrum 5 - Barcelona Supercomputing Center (BSC)


## Citing Qdislib


If you use Qdislib in a scientific publication, we would appreciate citations to the following paper:

\M. Tejedor, B. Casas, J. Conejero, A. Cervera-Lierta and R. M. Badia, "Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing" in *ArXiv*, 2025, pp. 1-12

### Bibtex:

```latex

   @inproceedings{Qdislib,
               title       = {{Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing}},
               author      = {Mar Tejedor and Berta Cervera and Javier Conejero and Alba Cervera-Lierta and Rosa M. Badia},
               booktitle   = {ArXiv},
               pages       = {1-12},
               year        = {2025},
    }
```

## URLs

[Quantum Distributed computing libraries BSC](https://www.bsc.es/research-development/research-areas/distributed-computing/distributed-computing-libraries-and)


## Acknowledgements

The project acknowledges funding from the Spanish Ministry for Digital Transformation and of Civil Service of the Spanish Government through the QUANTUM ENIA project call - Quantum Spain, EU through the Recovery, Transformation and Resilience Plan – NextGenerationEU within the framework of the Digital Spain 2026. It acknowledges funding from Grant RYC2022-037769-I funded by MICIU/AEI/10.13039/501100011033 and by “ESF+”. It also acknowledges funding from projects  CEX2021-001148-S, and PID2023-147979NB-C21 from the  MCIN/AEI and MICIU/AEI /10.13039/501100011033 and by FEDER, UE, by the Departament de Recerca i Universitats de la Generalitat de Catalunya, research group MPiEDist (2021 SGR 00412).

## License

Apache License Version 2.0, see [LICENSE](LICENSE)
