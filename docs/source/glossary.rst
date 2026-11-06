Glossary of terms
=================

This page lists alphabetically the meaning of various terms used in
Qdislib's documentation:

.. these terms can be referenced across the documentation using the :term: environment

.. glossary::

    backend
      where Qdislib runs a :term:`subcircuit`: a simulator (Qiskit Aer, Qibo, CUDA-Q,
      PennyLane, QiliSDK), a GPU, or a real :term:`quantum processing unit`. Chosen per
      call with ``backend=``.

    circuit cutting
      decomposing a large quantum circuit into smaller subcircuits that can be
      executed independently, then recombining their results. Qdislib supports
      :term:`gate cutting`, :term:`Hadamard-test gate cutting`, :term:`wire cutting`,
      and any mix of them on one circuit (:term:`mixed cutting`).

    configuration
      one term of the combined quasiprobability decomposition of a set of cuts: a choice
      of one term per cut.

    distribution
      see :term:`output distribution`.

    expectation value
      the average value of an observable measured on a quantum state. Qdislib
      reconstructs the original circuit's expectation value from those of its
      subcircuits.

    fragment
      a part of a circuit that a cut leaves disconnected from the others. A fragment runs
      as several :term:`subcircuit`\ s, one per setting of the cuts that touch it.

    gamma
      the :term:`sampling overhead` of a cut.

    gate cutting
      a circuit-cutting strategy that replaces a two-qubit gate with a measure-and-prepare
      quasiprobability decomposition over local operations. A gate other than CZ or CX is cut
      as the ``rzz`` rotations it decomposes into, each at its own angle. See also
      :term:`Hadamard-test gate cutting`.

    golden cut point
      a wire leaving a fragment whose state is real. Its ``Y`` terms are zero and are left
      out: ``gamma`` 3 instead of 4.

    Hadamard-test gate cutting
      a :term:`gate cutting` method (Harrow and Lowe,
      `arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) that uses one local ancilla per
      fragment running a Hadamard test, instead of a mid-circuit measurement. Gates crossing one
      boundary can be cut together with one shared ancilla per fragment. Exposed as
      ``gate_cutting_hadamard``.

    mixed cutting
      cutting some gates (measure-and-prepare or Hadamard) *and* some wires on the
      same circuit in one workflow (``circuit_cutting``).

    output distribution
      the probability of each measurement bitstring, ``{bitstring: probability}``.
      Qdislib's ``*_sampling`` functions reconstruct it from the subcircuits, as an
      alternative to a single :term:`expectation value`.

    packing
      running several small circuits side by side, on separate groups of qubits of one device,
      in one job (``qd.pack``).

    quantum processing unit
      (QPU) real quantum hardware. Qdislib runs the :term:`subcircuit`\ s on a QPU
      when you pass it as ``backend=``.

    quasiprobability overhead
      see :term:`sampling overhead`.

    qubit
      Quantum bit.

    reconstruction
      the post-processing step that recombines the results of the subcircuits into the result
      of the original, uncut circuit. Qdislib contracts one tensor per :term:`fragment` over
      the cuts.

    sampling overhead
      the extra sampling cost of :term:`circuit cutting`, ``gamma``. For the same worst-case
      precision the shots grow as ``gamma**2``, and the ``gamma`` of several cuts multiply.
      A CZ cut has ``gamma`` 3, and a wire cut 4.

    semantic cache
      a cache keyed by a ZX-reduced graph hash of a circuit, so equivalent circuits can be
      found. In Qdislib it is ``cache_key="semantic"``.

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
