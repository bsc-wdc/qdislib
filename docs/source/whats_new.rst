What is new in 2.0
==================

This page summarises the release and how to move from 1.x. ``CHANGELOG.md`` in the repository has every change, with its reasons.

New
---

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - What
     - Where
   * - One engine and one calling convention for every cut type
     - :doc:`guide/concepts`
   * - ``circuit_cutting``: any mix of gate, Hadamard and wire cuts on one circuit
     - :doc:`guide/cut_types`
   * - The gate cut follows the gate's angle: ``gamma = 1 + 2|sin(phi)|``
     - :doc:`guide/cut_types`
   * - Joint Hadamard cut, joint wire measurement, and golden wire cut points
     - :doc:`guide/cut_types`
   * - A reconstruction that is a contraction, so many cuts reconstruct in seconds
     - :doc:`guide/cut_types`
   * - ``find_cut(mixed=True)``: a search over all cut types, priced as they run
     - :doc:`guide/finding_cuts`
   * - ``draw_cut`` draws the circuit and what the cut does to it
     - :doc:`guide/finding_cuts`
   * - ``cutting_overhead``: the cost of any cut, without running it
     - :doc:`guide/concepts`
   * - ``subs.circuits(software)``, ``reconstruct_counts``, ``result_from_counts`` and ``describe``
     - :doc:`guide/three_step`
   * - Pauli sums run once per commuting group, not once per term
     - :doc:`guide/observables`
   * - ``qubits=`` marginals, and ``FactoredDistribution`` over more than 24 qubits
     - :doc:`guide/distributions`
   * - ``error=`` bars, ``total_shots=``, ``num_samples=`` on every cut type, and ``shots_for``
     - :doc:`guide/shots`
   * - ``parameter_values=`` and ``subs.bind``: cut once, run for many values
     - :doc:`guide/parameters`
   * - The engine's cache: ``cache=``, ``use_cache``, ``cache_key``, ``cache_refresh``
     - :doc:`guide/caching`
   * - QiliSDK as a supported SDK, and dense ``unitary`` gates
     - :doc:`guide/backends`
   * - A backend needs one method, ``evaluate_samples``, and ``sample_with_qiskit`` writes it
     - :doc:`guide/backends`
   * - ``qd.pack`` and ``pack=True`` on the IBM backend
     - :doc:`guide/qpus`
   * - IBM and IQM give joint final and cut-point counts
     - :doc:`guide/qpus`
   * - ``circuits_per_task=`` and GPU tasks under PyCOMPSs, ``scripts/mn5/check.sh``
     - :doc:`guide/hpc`

Moving from 1.x
---------------

The 1.x keywords keep working in this release with a ``DeprecationWarning`` that names the replacement, and are removed in the next one.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - 1.x
     - 2.0
   * - ``observables=``
     - ``observable=`` (default: ``Z`` on every qubit, on every cut type)
   * - ``gates_cut=``, ``cut=``
     - ``cuts=``
   * - ``dag=``, ``rand_qc=``
     - ``circuit=``
   * - ``software=`` on a call that runs circuits
     - ``backend=``
   * - ``software=`` on ``*_subcircuits``
     - ``subs.circuits(software)``
   * - ``qpu=``, ``qpu_dict=``, ``qpu_config=``
     - ``backend=`` with a registered name or a ``build_backend(...)`` object
   * - ``gpu=``, ``gpu_min_qubits=``
     - ``backend=`` with ``extra={"use_gpu": True}``; by size, ``s.num_qubits`` in the three-step loop
   * - ``method=``, wire cutting's ``backend="numpy"``
     - ``backend=``
   * - ``gate_wire_cutting``
     - ``circuit_cutting``
   * - ``cut_cost``
     - ``cutting_overhead``
   * - ``gate_cutting_subcircuit_reconstruction``, ``gate_cutting_sampling_reconstruction``, ``wire_cutting_subcircuit_reconstruction``, ``gate_cutting_hadamard_subcircuit_reconstruction``
     - ``subs.reconstruct(results)``
   * - ``HadamardCutSubcircuits``, ``HadamardSampledCutSubcircuits``, ``SampledCutSubcircuits``, ``WireSampledCutSubcircuits``
     - ``CutCircuits``, which every ``*_subcircuits`` function returns. The old names warn.
   * - ``cache_mode=``
     - nothing: the cache keeps each circuit's whole result
   * - ``reduce_chunk=``
     - nothing: the reconstruction is one cheap task
   * - ``find_cut(implementation="ibm-ckt")``
     - ``find_cut(mixed=True)``. The old option and the ``qiskit-addon-cutting`` dependency are gone.

These break 1.x scripts:

- Only ``circuit``, ``cuts`` and ``observable`` are positional. A call in the old positional order, such as ``wire_cutting(qc, cut, True, "ZZ")``, raises a ``TypeError`` that says so.
- ``shots=1024`` is the default everywhere, and ``shots=None`` is exact. The Hadamard cut used to default to exact.
- ``joint=True`` is the default, so the Hadamard cut now cuts the gates crossing one boundary together.
- ``gate_cutting_hadamard_sampling`` means the output distribution, like every other ``_sampling`` function. The random-configuration estimate is now ``gate_cutting_hadamard(..., num_samples=N)``.
- The deprecated positional ``*_subcircuit_reconstruction`` functions cannot read the results of the new lists, nor do these have ``number_components``, ``weights`` or ``config_of``. They raise an error naming ``subs.reconstruct(results)``.
- A cut has to separate its two sides, with all the others applied, or it raises ``InvalidCutError``.
- A DAG's ``gate`` attribute uses OpenQASM 3 names, whatever SDK the circuit came from. Node names do not change.

The full change log
-------------------

.. mdinclude:: ../../CHANGELOG.md
