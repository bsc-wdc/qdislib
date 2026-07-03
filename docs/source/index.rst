.. Qdislib documentation master file, created by
   sphinx-quickstart on Tue Mar  5 09:41:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to Qdislib!
===================

**Qdislib** is a Python library designed for scalable quantum circuit execution using *circuit cutting* techniques. It enables the simulation of large quantum circuits by splitting them into smaller, manageable subcircuits that can be executed independently—either on classical simulators, GPUs, or quantum hardware.

Qdislib is built on top of the `PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar>`__ parallel runtime, allowing seamless distributed execution of quantum workloads across CPUs, GPUs, and QPUs.

With Qdislib, researchers and developers can:

- Perform **gate** and **wire cutting** to decompose complex quantum circuits.
- Leverage **GPU acceleration** using Qiskit Aer, cuQuantum or Qibojit.
- Submit subcircuits to **remote QPUs** like IBM Quantum.
- Work with circuits defined in **Qibo**, **Qiskit**, and **CUDA-Q**.
- Automatically identify good cut points with `find_cut`.
- Extract and manipulate subcircuits independently.

Whether you're targeting HPC systems, hybrid quantum-classical setups, or constrained simulators, Qdislib is a flexible and modular tool to bridge the gap between current hardware limitations and large-scale quantum algorithm design.

Explore the sections below to get started with installation, quickstart examples, user guides, API references, and more.

Qdislib has been implemented on top of
`PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`__ programming model,
and it is being developed by the
`Workflows and Distributed Computing <https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>`__
group of the `Barcelona Supercomputing Center <http://www.bsc.es>`_.

.. toctree::
   :maxdepth: 2
   :hidden:

   quickstart
   user_guide
   examples
   api_reference
   modules
   development
   glossary

.. nbgallery::
    :caption: Start with our notebooks!
    :name: rst-gallery

    notebooks/01_quickstart
    notebooks/02_cutting_workflow
    notebooks/03_finding_cuts
    notebooks/04_backends
    notebooks/05_observables
    notebooks/06_sampling
    notebooks/07_caching
    notebooks/08_backend_internals
    notebooks/09_distributed_and_qpu
    notebooks/applications/01_image_reconstruction
    notebooks/applications/02_vqe_h2
    notebooks/applications/03_scaling

Documentation
-------------

.. grid:: 2 2 3 3
    :gutter: 3

    .. grid-item-card:: 🚀 Quickstart
        :link: quickstart
        :link-type: doc

        Install and cut your first circuit in three steps.

    .. grid-item-card:: 📖 User guide
        :link: user_guide
        :link-type: doc

        Gate vs wire cutting, observables, sampling, caching, backends & GPUs.

    .. grid-item-card:: 🧪 Examples
        :link: examples
        :link-type: doc

        The hands-on notebook learning path and end-to-end applications.

    .. grid-item-card:: 🧩 API reference
        :link: api_reference
        :link-type: doc

        Every public function, class and module.

    .. grid-item-card:: 🛠️ Development
        :link: development
        :link-type: doc

        Build, test and release Qdislib.

    .. grid-item-card:: 📚 Glossary
        :link: glossary
        :link-type: doc

        Circuit-cutting terms, defined.

Source code
-----------

The source code of Qdislib is available online at `Github <https://github.com/bsc-wdc/Qdislib>`__.

Support
-------

If you have questions or issues about the Qdislib you can join us in `Slack <https://bit.ly/bsc-wdc-community>`__.

Alternatively, you can send us an e-mail to `support-compss@bsc.es <mailto:support-compss@bsc.es>`__.

Citing Qdislib
--------------

If you use Qdislib in a scientific publication, we would appreciate citations to
the following papers. Please cite the original paper (ACM, SC '25 Workshops) first:

\M. Tejedor, B. Casas, J. Conejero, A. Cervera-Lierta and R. M. Badia, "Orchestrating Quantum-HPC Workflows with Distributed Quantum Circuit Cutting," in *Proceedings of the SC '25 Workshops of the International Conference for High Performance Computing, Networking, Storage and Analysis*, ACM, 2025, pp. 1898-1906. https://doi.org/10.1145/3731599.3767547

\M. Tejedor, J. Conejero and R. M. Badia, "A Semantic Quantum Circuit Cache for Scalable and Distributed Quantum-Classical Workflows," arXiv:2604.26788, 2026. https://arxiv.org/abs/2604.26788

Bibtex
......

.. code:: latex

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


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
