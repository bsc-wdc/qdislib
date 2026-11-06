.. _three-step-label:

The three-step workflow
=======================

The all-in-one calls run the cut, the circuits and the reconstruction. The three-step form splits
them, so the middle step is yours: choose a device per subcircuit, submit them with your own
tools, batch, or save them for later. The notebook :doc:`../notebooks/02_cutting_workflow` runs it.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(6)
    circuit.ry(0.5, range(6))
    for q in range(5):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(6))
    cut = ["CZ_9"]

    # 1. Cut.
    subcircuits = qd.gate_cutting_subcircuits(circuit, cut, "Z" * 6)
    # The same for wire_cutting_, gate_cutting_hadamard_ and circuit_cutting_subcircuits.

    # 2. Run each subcircuit. backend= picks where.
    values = [qd.expectation_value(s, shots=8192, backend="qiskit") for s in subcircuits]

    # 3. Recombine.
    print(subcircuits.reconstruct(values))

What the list holds
-------------------

Every item is a :class:`~Qdislib.Subcircuit`, whatever the cut type:

- ``s.circuit``: the native circuit, built when asked for, in the input circuit's SDK;
- ``s.num_qubits``, ``s.software`` and ``s.observable``;
- ``s.draws``: with ``num_samples``, how many drawn configurations use it, to size its shots.

Running an item, with ``qd.expectation_value(s, shots=..., backend=..., seed=...)`` or
``qd.sample_subcircuit(...)``, returns a ``FragmentResult`` ``(corr, qubits)``: the Z correlators the
reconstruction reads, and the circuit qubits they refer to. ``reconstruct`` takes one per item, in
order. Choosing the backend per subcircuit is how you route by size:

.. code-block:: python

    values = [
        qd.expectation_value(s, shots=4000, backend="qibo" if s.num_qubits <= 3 else "qiskit")
        for s in subcircuits
    ]
    print(subcircuits.reconstruct(values))

Seeing what is passed between the steps
---------------------------------------

The list, each item and each result print what they hold. ``subcircuits.describe(values)`` is a
table of every subcircuit with its fragment, its qubits (``q2'`` is a wire cut's prepared line
carrying qubit 2, ``a0`` a Hadamard ancilla), the cut terms it runs and its values named by the qubits
and Pauli they read, with the reconstructed value at the end. ``.rows`` has the data, and in Jupyter it
shows as a table.

.. code-block:: python

    print(subcircuits)
    print(subcircuits[0])
    print(values[0])
    print(subcircuits.describe(values))

.. code-block:: text

    <CutCircuits: 12 subcircuits, 2 fragments, gate ['CZ_9'], gamma 3, observable ZZZZZZ>
    <Subcircuit: fragment 0, 3 qubits [q0 q1 q2], CZ_9 side 1: RZ(-π/2)>
    FragmentResult(Z[q0,q1,q2]: -0.0167; 1 correlator)

A value is the subcircuit's own expectation of its part of the observable. With a Pauli sum there is one
per term of its group. For a distribution list each result holds all ``2**n`` correlators, so the views
show only ``<Z>`` of each qubit, and the most likely outcomes of the reconstructed distribution.

Just the circuits, in any SDK
-----------------------------

To get the circuits and nothing else, ask the list for them, in ``"qiskit"``, ``"qibo"``, ``"cudaq"``,
``"pennylane"`` or ``"qilisdk"``. There is one native circuit per subcircuit, in order, on the
subcircuit's qubits. A gate cut's measuring side has a mid-circuit measurement;
``defer_measurements=True`` moves it onto an extra qubit read at the end, for hardware without
mid-circuit measurement. An SDK that cannot express it (QiliSDK) gets that circuit anyway.

.. code-block:: python

    circuits = subcircuits.circuits("qiskit", defer_measurements=True)
    print(len(circuits), circuits[0].num_qubits, circuits[4].num_qubits)

Reading back your own counts
----------------------------

Counts of circuits you ran yourself come back with ``reconstruct_counts``: one ``{bitstring: count}``
per circuit, every qubit measured at the end. ``bit_order`` is ``"q0_first"``, the default, or
``"q0_last"``, the order of Qiskit's ``get_counts()``.

.. code-block:: python

    from qiskit_aer import AerSimulator

    simulator = AerSimulator(seed_simulator=1)
    counts = []
    for c in circuits:
        c = c.copy()
        c.measure_all()
        counts.append(simulator.run(c, shots=8192).result().get_counts())
    print(subcircuits.reconstruct_counts(counts, bit_order="q0_last"))

``qd.result_from_counts(s, counts)`` converts one subcircuit's counts into the result ``reconstruct``
reads, to mix your own runs with Qdislib's. ``keep_counts=False`` keeps only the correlators the list's
observable reads. The default keeps the counts, so other observables can be read back. Counts of a
circuit whose gate cut still measures mid-circuit are refused: build it with
``defer_measurements=True``.

.. code-block:: python

    mixed = [
        qd.result_from_counts(s, c, bit_order="q0_last") if i % 2 else qd.expectation_value(s, shots=8192)
        for i, (s, c) in enumerate(zip(subcircuits, counts))
    ]
    print(subcircuits.reconstruct(mixed))
