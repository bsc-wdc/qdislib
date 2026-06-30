Glossary of terms
=================

This page lists alphabetically the meaning of various terms used in
Qdislib's documentation:

.. these terms can be referenced across the documentation using the :term: environment

.. glossary::

    circuit cutting
      decomposing a large quantum circuit into smaller subcircuits that can be
      executed independently, then recombining their results. Qdislib supports
      two strategies: :term:`gate cutting` and :term:`wire cutting`.

    expectation value
      the average value of an observable measured on a quantum state. Qdislib
      reconstructs the original circuit's expectation value from those of its
      subcircuits.

    gate cutting
      a circuit-cutting strategy that replaces a selected (two-qubit) gate with a
      quasiprobability decomposition over single-qubit operations, splitting the
      circuit into independent subcircuits.

    qubit
      Quantum bit.

    reconstruction
      the post-processing step that recombines the per-subcircuit
      :term:`expectation value`\ s (with the appropriate signs and normalization)
      into the expectation value of the original, uncut circuit.

    subcircuit
      one of the smaller circuits produced by cutting a larger circuit; it can be
      executed independently of the others.

    synchronization
      when the execution of a parallel application stalls until certain data is
      generated, or until all tasks have finished.

    task
      a unit of computation that can be executed in a remote computer.

    wire cutting
      a circuit-cutting strategy that cuts a "wire" (a qubit's worldline between
      two gates) using measurement and state preparation, disconnecting the
      circuit into independent subcircuits.