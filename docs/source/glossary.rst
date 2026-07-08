Glossary of terms
=================

This page lists alphabetically the meaning of various terms used in
Qdislib's documentation:

.. these terms can be referenced across the documentation using the :term: environment

.. glossary::

    backend
      the execution target Qdislib runs a :term:`subcircuit` on -- a CPU simulator
      (Qibo / Qiskit / CUDA-Q / PennyLane), a GPU, or a real :term:`quantum processing unit`.
      Selected per call with the ``software`` / ``gpu`` / ``qpu`` arguments.

    circuit cutting
      decomposing a large quantum circuit into smaller subcircuits that can be
      executed independently, then recombining their results. Qdislib supports
      :term:`gate cutting`, :term:`Hadamard-test gate cutting`, :term:`wire cutting`,
      and :term:`mixed cutting` (gate + wire on one circuit).

    expectation value
      the average value of an observable measured on a quantum state. Qdislib
      reconstructs the original circuit's expectation value from those of its
      subcircuits.

    gate cutting
      a circuit-cutting strategy that replaces a selected (two-qubit) gate with a
      measure-and-prepare quasiprobability decomposition over single-qubit
      operations, splitting the circuit into independent subcircuits. Cuts any
      two-qubit gate (``CZ``/``CX`` natively; others decomposed into ``CZ`` +
      single-qubit gates). See also :term:`Hadamard-test gate cutting`.

    Hadamard-test gate cutting
      an alternative :term:`gate cutting` method (Harrow & Lowe,
      `arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) that cuts a two-qubit
      gate with **one local ancilla** per fragment running a Hadamard test instead of
      a mid-circuit measurement. It cuts **any** gate as a single cut at the optimal
      quasiprobability overhead, and (having no feed-forward) keeps the fragments
      fully independent. Exposed as ``gate_cutting_hadamard``.

    mixed cutting
      cutting some gates *and* some wires on the same circuit in one workflow
      (``gate_wire_cutting``); the two quasiprobability decompositions tensor and all
      fragments run in parallel.

    output distribution
      the probability of each measurement bitstring, ``{bitstring: probability}``.
      Qdislib's ``*_sampling`` functions reconstruct it from the subcircuits, as an
      alternative to a single :term:`expectation value`.

    quantum processing unit
      (QPU) real quantum hardware. Qdislib can route each :term:`subcircuit` to a
      QPU with ``qpu=True`` and a ``qpu_dict``.

    quasiprobability overhead
      the extra sampling cost of :term:`circuit cutting`: each cut expands into
      several terms (6 for a :term:`gate cutting` cut, 8 for a :term:`wire cutting`
      cut), so the number of subcircuit evaluations -- and the shots needed for a
      given accuracy -- grows exponentially in the number of cuts.

    qubit
      Quantum bit.

    reconstruction
      the post-processing step that recombines the per-subcircuit
      :term:`expectation value`\ s (with the appropriate signs and normalization)
      into the expectation value of the original, uncut circuit.

    semantic cache
      a content-addressed cache that keys each :term:`subcircuit` by a ZX-reduced
      graph hash, so *functionally equivalent* subcircuits are simulated only once.

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