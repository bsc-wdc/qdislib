Notebooks and examples
======================

Running the notebooks
---------------------

Every notebook runs top to bottom on a laptop in a few minutes at most, with no credentials and no GPU, and its stored outputs came from a real run. Each page has the links to open it in **Colab**, whose first cell installs what it needs, or **Binder**, whose image has Qdislib and COMPSs ready.

Some have sliders. A slider on a page like this one moves through figures rendered when the notebook ran, so it works with no kernel; in Colab or Binder you can change what it renders, and cells with live widgets recompute at each move.

To run them locally:

.. code-block:: bash

    git clone https://github.com/bsc-wdc/qdislib.git
    cd qdislib
    pip install ".[qibo,qiskit,cache]" jupyterlab ipywidgets    # or: uv pip install ...
    jupyter lab notebooks

A first look
------------

:doc:`notebooks/Qdislib_try_it` goes through a cut in four steps and compares the cutting methods, in a few minutes. Its results come from a real run under COMPSs. See also :doc:`try`.

Learning path
-------------

Read these in order. Each one stands on its own, with its own imports and example circuit.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Notebook
     - What it covers
   * - :doc:`notebooks/01_quickstart`
     - A first cut, what it costs, and what to read next
   * - :doc:`notebooks/02_cutting_workflow`
     - The three steps; run the subcircuits yourself and read back your counts
   * - :doc:`notebooks/03_finding_cuts`
     - ``find_cut`` and its limits, ``mixed=True``, ``coupling_map``, ``draw_cut`` (with a slider)
   * - :doc:`notebooks/04_backends`
     - Qiskit, Qibo and the other SDKs; ``from_qasm``; dense unitary gates
   * - :doc:`notebooks/05_observables`
     - Pauli strings and Hamiltonians; terms in commuting groups
   * - :doc:`notebooks/06_sampling`
     - Output distributions, ``qubits=``, ``FactoredDistribution`` at 40 qubits (with a slider)
   * - :doc:`notebooks/07_caching`
     - The engine's cache, on a Hamiltonian and in a variational loop
   * - :doc:`notebooks/08_backend_internals`
     - Your own backend in one method, on a noisy simulator
   * - :doc:`notebooks/09_distributed_and_qpu`
     - PyCOMPSs tasks, ``sync=False``, ``circuits_per_task=``, the task graph; run under real COMPSs. QPU code, not run

Going further
-------------

Read any of these after notebook 02.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Notebook
     - What it covers
   * - :doc:`notebooks/10_cut_types`
     - Gate, Hadamard and wire cuts compared, and any mix with ``circuit_cutting`` (with a slider)
   * - :doc:`notebooks/11_shots_and_errors`
     - ``error=`` checked against repeated runs, ``total_shots=``, ``num_samples=``, ``shots_for`` (with a slider and live widgets)
   * - :doc:`notebooks/12_parameters`
     - ``parameter_values=``, ``subs.bind``, and an optimiser on one cut
   * - :doc:`notebooks/13_packing`
     - ``qd.pack`` on a fake backend, with no credentials

Applications
------------

Larger, self-contained studies.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Notebook
     - What it covers
   * - :doc:`notebooks/applications/01_vqe_h2`
     - A small VQE: the H2 energy from a cut ansatz, with a Hamiltonian and an optimiser
   * - :doc:`notebooks/applications/02_scaling`
     - N2 on 12 qubits, and a 25-qubit lattice in blocks, with the cost of each
   * - :doc:`notebooks/applications/03_image_reconstruction`
     - An image in 8 qubits, rebuilt from two halves (with a slider)

Scripts
-------

``examples/`` in the repository holds scripts for the distributed path. Each docstring says how to run it.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Script
     - What it does
   * - ``examples/compss_app.py``
     - Every cut type on a small circuit against the exact value, serially or under ``runcompss``. Prints ``RESULT: PASS``.
   * - ``examples/mixed_cutting_compss.py``
     - A gate cut and a wire cut on one circuit, for the task graph
   * - ``examples/hadamard_cut_compss.py``
     - The Hadamard cut: all configurations, drawn ones and three-step, overlapping with ``sync=False``
   * - ``examples/verify_sampling_pycompss.py``
     - Distributions under ``runcompss`` against the exact one
   * - ``examples/verify_engine_pycompss.py``
     - Backend objects, groups, marginals, error bars, parameters and a 26-qubit distribution under ``runcompss``
   * - ``examples/verify_gpu_pycompss.py``
     - GPU tasks under ``runcompss``, on GPU nodes
   * - ``examples/qpu_listener.py``
     - The listener for BSC QPUs, on the quantum cluster

.. toctree::
   :maxdepth: 1
   :caption: Learning path
   :hidden:

   notebooks/Qdislib_try_it
   notebooks/01_quickstart
   notebooks/02_cutting_workflow
   notebooks/03_finding_cuts
   notebooks/04_backends
   notebooks/05_observables
   notebooks/06_sampling
   notebooks/07_caching
   notebooks/08_backend_internals
   notebooks/09_distributed_and_qpu

.. toctree::
   :maxdepth: 1
   :caption: Going further
   :hidden:

   notebooks/10_cut_types
   notebooks/11_shots_and_errors
   notebooks/12_parameters
   notebooks/13_packing

.. toctree::
   :maxdepth: 1
   :caption: Applications
   :hidden:

   notebooks/applications/01_vqe_h2
   notebooks/applications/02_scaling
   notebooks/applications/03_image_reconstruction
