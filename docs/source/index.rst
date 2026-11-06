Qdislib
=======

**Qdislib** runs quantum circuits by **cutting** them into smaller subcircuits, running those independently, and reconstructing the result of the original circuit. It runs serially on a laptop or in parallel on HPC with `PyCOMPSs <https://compss-doc.readthedocs.io/>`__, on simulators, GPUs and QPUs.

.. code-block:: bash

    pip install "qdislib[qiskit]"        # or: uv pip install "qdislib[qiskit]"

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(10)
    circuit.ry(0.5, range(10))
    for q in range(9):
        circuit.cz(q, q + 1)

    cut = qd.find_cut(circuit, max_qubits=5, wire_cut=False)    # ['CZ_15']: two fragments

    # In one call
    print(qd.gate_cutting(circuit, cut, shots=8192, seed=1, error=True))

    # Or in three steps, to run the subcircuits where you want
    subcircuits = qd.gate_cutting_subcircuits(circuit, cut)
    values = [qd.expectation_value(s, shots=8192, seed=1) for s in subcircuits]
    print(subcircuits.reconstruct(values))

Use it to run a circuit too wide for your device or simulator, when a few cuts split it. Each cut multiplies the shots a given precision needs, by ``gamma**2``, so the number of cuts is what limits it: ``qd.cutting_overhead`` tells you the cost before anything runs.

.. tip::

   New to Qdislib? :doc:`Try it in your browser <try>`, with nothing to install, or :doc:`explore the whole library <explore>`.

.. grid:: 1 2 4 4
    :gutter: 3

    .. grid-item-card:: Quickstart
        :link: quickstart
        :link-type: doc

        Install, and cut your first circuit.

    .. grid-item-card:: Notebooks
        :link: examples
        :link-type: doc

        A learning path with results and figures. Open any of them in Colab or Binder.

    .. grid-item-card:: User guide
        :link: user_guide
        :link-type: doc

        What each part does and when to use it.

    .. grid-item-card:: PyCOMPSs and HPC
        :link: guide/hpc
        :link-type: doc

        Running in parallel, GPU tasks and MareNostrum 5.

    .. grid-item-card:: API reference
        :link: api_reference
        :link-type: doc

        Every public function and class.

    .. grid-item-card:: What is new in 2.0
        :link: whats_new
        :link-type: doc

        What changed, and how to move from 1.x.

    .. grid-item-card:: Cheat sheet
        :link: Qdislib_cheatsheet.pdf
        :link-type: url

        The main calls on two pages (PDF).

    .. grid-item-card:: Previous paper
        :link: https://dl.acm.org/doi/10.1145/3731599.3767547
        :link-type: url

        Distributed circuit cutting for quantum-HPC workflows (SC '25 Workshops).

What it does
------------

- **Three kinds of cut, in any mix** (``circuit_cutting``): measure-and-prepare gate cuts that follow the gate's angle, Hadamard gate cuts with no mid-circuit measurement, and wire cuts.
- **Expectation values and distributions.** Pauli strings and Hamiltonians; output distributions, marginals with ``qubits=``, and ``FactoredDistribution`` over more than 24 qubits.
- **Finding cuts** (``find_cut``): for a size, by cost over all kinds of cut (``mixed=True``), or for a device's connectivity (``coupling_map``).
- **Shots and error bars.** ``error=``, ``total_shots=``, ``num_samples=`` and ``shots_for``.
- **Any SDK.** Circuits from Qiskit, Qibo, CUDA-Q, PennyLane and QiliSDK, run on any of them. ``qd.BaseBackend`` adds your own in one method.
- **PyCOMPSs.** Every fragment circuit is a task, including GPU tasks. Without PyCOMPSs it runs serially, with the same results.
- **QPUs.** IBM, IQM and BSC devices, and ``qd.pack`` to run many subcircuits in one job.
- **A cache** of fragment results, across observables and runs.

Notebooks
---------

Every :doc:`notebook <examples>` runs top to bottom on a laptop with no credentials, stores its outputs, and has links to open it in Colab or Binder. Some have sliders that work on the page itself.

- **Start here:** :doc:`notebooks/Qdislib_try_it` follows the poster's four steps, then :doc:`notebooks/01_quickstart`, then :doc:`notebooks/02_cutting_workflow` and :doc:`notebooks/03_finding_cuts`.
- **Topics:** :doc:`observables <notebooks/05_observables>`, :doc:`distributions <notebooks/06_sampling>`, :doc:`caching <notebooks/07_caching>`, :doc:`cut types <notebooks/10_cut_types>`, :doc:`error bars <notebooks/11_shots_and_errors>`, :doc:`free parameters <notebooks/12_parameters>`.
- **Running it elsewhere:** :doc:`other SDKs <notebooks/04_backends>`, :doc:`your own backend <notebooks/08_backend_internals>`, :doc:`PyCOMPSs and QPUs <notebooks/09_distributed_and_qpu>`, :doc:`packing on a QPU <notebooks/13_packing>`.
- **Applications:** :doc:`a VQE <notebooks/applications/01_vqe_h2>`, :doc:`scaling up <notebooks/applications/02_scaling>`, :doc:`an image <notebooks/applications/03_image_reconstruction>`.

.. toctree::
   :maxdepth: 2
   :hidden:

   try
   explore
   quickstart
   examples
   user_guide
   poster
   api_reference
   whats_new
   modules
   development
   glossary

Source code
-----------

The source code of Qdislib is on `GitHub <https://github.com/bsc-wdc/qdislib>`__. Qdislib is developed by the `Workflows and Distributed Computing <https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>`__ group of the `Barcelona Supercomputing Center <https://www.bsc.es>`_.

Support
-------

For questions or issues, join us on `Slack <https://bit.ly/bsc-wdc-community>`__, or write to `support-compss@bsc.es <mailto:support-compss@bsc.es>`__.

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
