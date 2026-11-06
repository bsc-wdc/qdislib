Concepts and calling convention
===============================

What a cut does
---------------

A cut replaces part of a circuit with a weighted sum of simpler operations that act on
one side only. Run over all the terms of that sum, the circuit falls apart into
**fragments**, independent pieces that are smaller than the circuit. Each fragment runs
as several **subcircuits**, one per term of the cuts that touch it, and the result of the
original circuit is **reconstructed** from theirs.

The price is the **sampling overhead** ``gamma``. For the same worst-case precision the
shots grow as ``gamma**2``, and the ``gamma`` of several cuts multiply. A CZ gate cut has
``gamma`` 3, so ten of them have ``gamma`` 59,049. Cutting is useful when a few cuts split
the circuit: keep the number of cuts low, and use ``find_cut`` to find where.

``cutting_overhead`` reports all of this without running anything:

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(5)
    circuit.ry(0.5, range(5))
    for q in range(4):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(5))

    print(qd.cutting_overhead(circuit, ["CZ_7"]))
    # {'subcircuits': 12, 'configurations': 6, 'fragments': 2, 'gamma': 3.0}

``configurations`` is the number of terms of the combined sum. Here the two fragments run
6 settings each, so 12 subcircuits.

Calling convention
------------------

Every cut type has the same four functions with the same arguments:

.. code-block:: text

    X(circuit, cuts, observable=None, *, backend=None, shots=1024, seed=None, sync=True, ...)
    X_sampling(circuit, cuts, *, backend=None, shots=1024, seed=None, sync=True, ...)
    X_subcircuits(circuit, cuts, observable=None, *, ...)
    X_sampling_subcircuits(circuit, cuts, *, ...)

where ``X`` is ``gate_cutting``, ``wire_cutting``, ``gate_cutting_hadamard`` or
``circuit_cutting``.

- ``X`` returns the expectation value of ``observable``, a Pauli string or a weighted
  Pauli sum; ``None`` is ``Z`` on every qubit.
- ``X_sampling`` returns the output distribution ``{bitstring: probability}``, qubit 0
  first in each key.
- ``X_subcircuits`` and ``X_sampling_subcircuits`` return the list of subcircuits, for the
  :doc:`three-step workflow <three_step>`.

The options have the same name and default everywhere:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Meaning
   * - ``backend=None``
     - Where the subcircuits run: an SDK name, a registered backend name or a backend
       object. ``None`` is the circuit's own SDK's simulator. See :doc:`backends`.
   * - ``shots=1024``
     - Shots per fragment circuit. ``None`` is exact (a statevector), on every cut type.
   * - ``seed=None``
     - Makes a simulator's shots reproducible. Each subcircuit gets its own stream.
   * - ``sync=True``
     - ``False`` returns a future under PyCOMPSs, to overlap several cuts. See :doc:`hpc`.
   * - ``optimized=False``
     - Each fragment runs once per setting of the cuts that touch it (gate, wire and
       ``circuit_cutting``). See :doc:`cut_types`.
   * - ``joint=True``
     - Joint Hadamard cut and joint wire measurement (wire, Hadamard and
       ``circuit_cutting``).
   * - ``hadamard=False``
     - In ``circuit_cutting``, cut the gates of a list with the Hadamard test.
   * - ``num_samples=None``
     - Draw this many configurations instead of running them all. See :doc:`shots`.
   * - ``error=False``
     - Return a standard error too. See :doc:`shots`.
   * - ``total_shots=None``
     - One budget for the whole call, instead of ``shots`` per circuit. See :doc:`shots`.
   * - ``qubits=None``
     - On ``_sampling``: the marginal over these qubits. See :doc:`distributions`.
   * - ``cache=None``
     - Reuse fragment results. See :doc:`caching`.
   * - ``parameter_values=None``
     - One result per set of values of a circuit's free parameters. See :doc:`parameters`.
   * - ``circuits_per_task=None``
     - Several fragment circuits per PyCOMPSs task. See :doc:`hpc`.

Only ``circuit``, ``cuts`` and ``observable`` go by position; the rest by keyword.

Naming what to cut
------------------

A gate is named by its gate name and its **position in the circuit, counting from 1**.
``"CZ_8"`` is the 8th instruction, a CZ. Any SDK's spelling works: ``"CX_7"``, ``"CNOT_7"``
and ``"cx_7"`` are the same gate, and naming a gate that is not there says what is there.
A wire cut is the edge between two consecutive gates on one qubit, ``("RY_3", "CZ_8")``.
:func:`~Qdislib.find_cut` returns these names, and :func:`~Qdislib.draw_cut` prints each
gate's index on its diagram.

Inside, every circuit becomes the same DAG, whose gates are named and defined as in
OpenQASM 3 (``stdgates.inc``: ``h``, ``x``, ``rz``, ``p``, ``u``, ``cx``, ``cz``, ``cp``,
``cu``, ``swap``, ``ccx``, ..., plus ``rxx``, ``ryy``, ``rzz``, ``iswap``). Each SDK's
conventions are converted on the way in, so the same circuit gives the same DAG in every
SDK.

Coming from Qdislib 1.x
-----------------------

The 1.x names still work in this release with a ``DeprecationWarning``, and go away in the
next one. :doc:`../whats_new` has the table of replacements.
