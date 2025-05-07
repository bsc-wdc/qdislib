.. Qdislib documentation master file, created by
   sphinx-quickstart on Tue Mar  5 09:41:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to Qdislib!
===================

.. image:: ../logo.png
   :alt: Qdislib Logo
   :align: left
   :scale: 100%

**Qdislib** is a Python library designed for scalable quantum circuit execution using *circuit cutting* techniques. It enables the simulation of large quantum circuits by splitting them into smaller, manageable subcircuits that can be executed independently—either on classical simulators, GPUs, or quantum hardware.

Qdislib is built on top of the `PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar>`_ parallel runtime, allowing seamless distributed execution of quantum workloads across CPUs, GPUs, and QPUs.

With Qdislib, researchers and developers can:

- Perform **gate** and **wire cutting** to decompose complex quantum circuits.
- Leverage **GPU acceleration** using cuQuantum or Qibojit.
- Submit subcircuits to **remote QPUs** like IBM Quantum.
- Work with circuits defined in both **Qibo** and **Qiskit**.
- Automatically identify good cut points with `find_cut`.
- Extract and manipulate subcircuits independently.

Whether you're targeting HPC systems, hybrid quantum-classical setups, or constrained simulators, Qdislib is a flexible and modular tool to bridge the gap between current hardware limitations and large-scale quantum algorithm design.

Explore the sections below to get started with installation, quickstart examples, user guides, API references, and more.


Qdislib has been implemented on top of
`PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`_ programming model,
and it is being developed by the
`Workflows and Distributed Computing <https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>`_
group of the `Barcelona Supercomputing Center <http://www.bsc.es>`_.


Documentation
-------------

* :doc:`Quickstart <quickstart>`
* :doc:`User guide <user_guide>`
* :doc:`API Reference <api_reference>`
* :doc:`Modules <modules>`
* :doc:`Examples <examples>`
* :doc:`Development <development>`
* :doc:`Glossary <glossary>`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api_reference
   development
   examples
   glossary
   modules
   quickstart
   user_guide

Source code
-----------

The source code of Qdislib is available online at `Github <https://github.com/bsc-wdc/Qdislib>`_.

Support
-------

If you have questions or issues about the Qdislib you can join us in `Slack <https://bit.ly/bsc-wdc-community>`_.

Alternatively, you can send us an e-mail to `support-compss@bsc.es <mailto:support-compss@bsc.es>`_.

Citing dislib
-------------

If you use Qdislib in a scientific publication, we would appreciate citations to the following paper:

\M. Tejedor, B. Casas, J. Conejero, A. Cervera-Lierta and R. M. Badia, "Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing" in *https://arxiv.org/abs/2505.01184*, 2025, pp. 1-12

Bibtex
......

.. code:: latex

   @inproceedings{Qdislib,
               title       = {{Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing}},
               author      = {Mar Tejedor and Berta Cervera and Javier Conejero and Alba Cervera-Lierta and Rosa M. Badia},
               booktitle   = {https://arxiv.org/abs/2505.01184},
               pages       = {1-12},
               year        = {2025},
    }


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
