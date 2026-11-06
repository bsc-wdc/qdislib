User guide
==========

These pages explain what each part of Qdislib does and when to use it. Each code block ran against this version. For the same material with results and figures, see the :doc:`notebooks <examples>`.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - To...
     - Read
   * - understand what a cut costs and how every function is called
     - :doc:`guide/concepts`
   * - choose between gate, Hadamard and wire cuts, or mix them
     - :doc:`guide/cut_types`
   * - let Qdislib choose where to cut
     - :doc:`guide/finding_cuts`
   * - run the subcircuits yourself
     - :doc:`guide/three_step`
   * - measure X, Y and Hamiltonians
     - :doc:`guide/observables`
   * - rebuild distributions, marginals, and above 24 qubits
     - :doc:`guide/distributions`
   * - get error bars and spend shots well
     - :doc:`guide/shots`
   * - run one circuit for many values of its parameters
     - :doc:`guide/parameters`
   * - reuse results across runs
     - :doc:`guide/caching`
   * - use another SDK, a GPU or your own device
     - :doc:`guide/backends`
   * - run on IBM, IQM or BSC hardware, or pack many subcircuits in one job
     - :doc:`guide/qpus`
   * - run in parallel with PyCOMPSs, or on MareNostrum 5
     - :doc:`guide/hpc`

.. toctree::
   :maxdepth: 2
   :hidden:

   guide/concepts
   guide/cut_types
   guide/finding_cuts
   guide/three_step
   guide/observables
   guide/distributions
   guide/shots
   guide/parameters
   guide/caching
   guide/backends
   guide/qpus
   guide/hpc
