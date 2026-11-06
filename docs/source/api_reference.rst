API reference
=============

Everything here is available as ``qd.<name>`` after ``import Qdislib as qd``. The :doc:`module pages <modules>` hold the same objects by module, along with the internals.

.. currentmodule:: Qdislib

Cutting: all in one call
------------------------

Cut, run and reconstruct. :doc:`guide/concepts` has the arguments they share.

.. autosummary::

   gate_cutting
   wire_cutting
   gate_cutting_hadamard
   circuit_cutting
   gate_cutting_sampling
   wire_cutting_sampling
   gate_cutting_hadamard_sampling
   circuit_cutting_sampling
   estimate

.. autofunction:: Qdislib.circuit_cutting
.. autofunction:: Qdislib.circuit_cutting_sampling
.. autofunction:: Qdislib.gate_cutting
.. autofunction:: Qdislib.gate_cutting_sampling
.. autofunction:: Qdislib.wire_cutting
.. autofunction:: Qdislib.wire_cutting_sampling
.. autofunction:: Qdislib.gate_cutting_hadamard
.. autofunction:: Qdislib.gate_cutting_hadamard_sampling
.. autofunction:: Qdislib.estimate

Cutting: the three-step workflow
--------------------------------

See :doc:`guide/three_step`.

.. autosummary::

   circuit_cutting_subcircuits
   circuit_cutting_sampling_subcircuits
   gate_cutting_subcircuits
   gate_cutting_sampling_subcircuits
   wire_cutting_subcircuits
   wire_cutting_sampling_subcircuits
   gate_cutting_hadamard_subcircuits
   gate_cutting_hadamard_sampling_subcircuits
   expectation_value
   sample_subcircuit
   result_from_counts
   Subcircuit
   CutCircuits
   CombinedCutSubcircuits

.. autofunction:: Qdislib.circuit_cutting_subcircuits
.. autofunction:: Qdislib.circuit_cutting_sampling_subcircuits
.. autofunction:: Qdislib.gate_cutting_subcircuits
.. autofunction:: Qdislib.gate_cutting_sampling_subcircuits
.. autofunction:: Qdislib.wire_cutting_subcircuits
.. autofunction:: Qdislib.wire_cutting_sampling_subcircuits
.. autofunction:: Qdislib.gate_cutting_hadamard_subcircuits
.. autofunction:: Qdislib.gate_cutting_hadamard_sampling_subcircuits
.. autofunction:: Qdislib.expectation_value
.. autofunction:: Qdislib.sample_subcircuit
.. autofunction:: Qdislib.result_from_counts
.. autoclass:: Qdislib.Subcircuit
   :members:
.. autoclass:: Qdislib.CutCircuits
   :members:
.. autoclass:: Qdislib.CombinedCutSubcircuits
   :members:

Results
-------

.. autosummary::

   Estimate
   Distribution
   FactoredDistribution
   shots_for

.. autoclass:: Qdislib.Estimate
.. autoclass:: Qdislib.Distribution
.. autoclass:: Qdislib.FactoredDistribution
   :members:
.. autofunction:: Qdislib.shots_for

Finding and inspecting cuts
---------------------------

See :doc:`guide/finding_cuts`.

.. autosummary::

   find_cut
   draw_cut
   cutting_overhead
   analytical_solution
   exact_statevector

.. autofunction:: Qdislib.find_cut
.. autofunction:: Qdislib.draw_cut
.. autofunction:: Qdislib.cutting_overhead
.. autofunction:: Qdislib.analytical_solution
.. autofunction:: Qdislib.exact_statevector

Backends
--------

See :doc:`guide/backends` and :doc:`guide/qpus`.

.. autosummary::

   BaseBackend
   BackendConfig
   build_backend
   build_execution_backend
   register_backend
   available_backends
   sample_with_qiskit
   QdBackEnd

.. autoclass:: Qdislib.BaseBackend
   :members: evaluate_samples, evaluate_joint_counts, evaluate_joint_counts_many, evaluate_expectation, packs, uses_gpu
.. autoclass:: Qdislib.BackendConfig
.. autofunction:: Qdislib.build_backend
.. autofunction:: Qdislib.build_execution_backend
.. autofunction:: Qdislib.register_backend
.. autofunction:: Qdislib.available_backends
.. autofunction:: Qdislib.sample_with_qiskit
.. autoclass:: Qdislib.QdBackEnd
   :members:

Packing
-------

See :doc:`guide/qpus`.

.. autosummary::

   pack
   PackedCircuits

.. autofunction:: Qdislib.pack
.. autoclass:: Qdislib.PackedCircuits
   :members:

Cache
-----

See :doc:`guide/caching`.

.. autosummary::

   use_cache
   default_circuit_cache
   semantic_key
   CircuitCache
   InMemoryCircuitCache
   LMDBCircuitCache
   RedisCircuitCache

.. autofunction:: Qdislib.use_cache
.. autofunction:: Qdislib.default_circuit_cache
.. autofunction:: Qdislib.semantic_key
.. autoclass:: Qdislib.CircuitCache
   :members:
.. autoclass:: Qdislib.InMemoryCircuitCache
.. autoclass:: Qdislib.LMDBCircuitCache
.. autoclass:: Qdislib.RedisCircuitCache

Other
-----

.. autosummary::

   from_qasm
   generate_random_graph
   set_log_level

.. autofunction:: Qdislib.from_qasm
.. autofunction:: Qdislib.generate_random_graph
.. autofunction:: Qdislib.set_log_level

Exceptions
----------

- :class:`~Qdislib.utils.exceptions.QdislibException`
- :class:`~Qdislib.utils.exceptions.UnsupportedBackendError`
- :class:`~Qdislib.utils.exceptions.InvalidCutError`
- :class:`~Qdislib.utils.qpu_bridge.QPUExecutionError`
- :class:`~Qdislib.utils.qpu_bridge.QPUTimeoutError`

Deprecated
----------

These work in this release with a ``DeprecationWarning``, and are removed in the next one. :doc:`whats_new` has their replacements.

.. autosummary::

   gate_wire_cutting
   cut_cost
   gate_cutting_subcircuit_reconstruction
   gate_cutting_sampling_reconstruction
   wire_cutting_subcircuit_reconstruction
   gate_cutting_hadamard_subcircuit_reconstruction

.. autofunction:: Qdislib.gate_wire_cutting
.. autofunction:: Qdislib.cut_cost
.. autofunction:: Qdislib.gate_cutting_subcircuit_reconstruction
.. autofunction:: Qdislib.gate_cutting_sampling_reconstruction
.. autofunction:: Qdislib.wire_cutting_subcircuit_reconstruction
.. autofunction:: Qdislib.gate_cutting_hadamard_subcircuit_reconstruction
