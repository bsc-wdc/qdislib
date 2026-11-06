Quickstart
==========

Installation
------------

Qdislib needs Python 3.10 or later and at least one quantum SDK, installed as an extra:

.. code-block:: bash

    pip install "qdislib[qiskit]"     # or [qibo], [cudaq], [pennylane], [qilisdk], [all]

``[all]`` is all of the above plus the semantic cache (``[cache]``). ``[iqm]`` for IQM hardware and ``[qpu]`` for the BSC QPU listener are separate. The library imports with no SDK installed, but you need one to run circuits.

Qdislib runs serially as it is. To run in parallel, install `PyCOMPSs <https://compss-doc.readthedocs.io/>`__ 3.3 or later, which ``pip`` does not do for you: see :doc:`guide/hpc`.

Or skip the install: every :doc:`notebook <examples>` has a link to open it in Colab or Binder.

A first cut
-----------

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(5)
    for qubit, angle in enumerate([0.8, 0.5, 0.6, 0.4, 0.3]):
        circuit.ry(angle, qubit)
    for qubit in range(4):
        circuit.cz(qubit, qubit + 1)

    cut = qd.find_cut(circuit, max_qubits=3, wire_cut=False)
    print(cut)                              # ['CZ_7']
    print(qd.cutting_overhead(circuit, cut))
    # {'subcircuits': 12, 'configurations': 6, 'fragments': 2, 'gamma': 3.0}

    value = qd.gate_cutting(circuit, cut, shots=8192, seed=1, error=True)
    print(value)                            # 0.441427 ± 0.0082
    print(qd.analytical_solution(circuit, "ZZZZZ"))   # 0.4440, exact

``CZ_7`` is the 7th instruction of the circuit, a CZ. With no observable, ``gate_cutting`` reads ``Z`` on every qubit. ``shots=None`` would give the exact value, and ``qd.draw_cut(circuit, cut)`` a picture of the cut.

Three steps
-----------

Split the call in three when you want to run the subcircuits yourself, choose a device for each, or inspect them:

.. code-block:: python

    # 1. Cut the circuit into independently executable subcircuits.
    subcircuits = qd.gate_cutting_subcircuits(circuit, cut, "ZZZZZ")

    # 2. Run each one. This loop is yours.
    values = [qd.expectation_value(s, shots=8192, seed=i) for i, s in enumerate(subcircuits)]

    # 3. Recombine.
    print(subcircuits.reconstruct(values))

``subcircuits.circuits("qiskit")`` gives the native circuits, in that or any other SDK, and ``backend=`` on ``expectation_value`` runs one anywhere: see :doc:`guide/three_step`.

Other cuts
----------

.. code-block:: python

    # a wire cut: the line of qubit 2, between two gates
    print(qd.wire_cutting(circuit, [("CZ_7", "CZ_8")], shots=8192, seed=1))

    # a Hadamard gate cut: no mid-circuit measurement, one ancilla per fragment
    print(qd.gate_cutting_hadamard(circuit, ["CZ_7"], shots=8192, seed=1))

    # any mix, and a search over all of them priced by cost
    cuts = qd.find_cut(circuit, max_qubits=3, mixed=True)
    print(qd.circuit_cutting(circuit, cuts, shots=8192, seed=1))

    # the output distribution, here over two qubits
    print(qd.circuit_cutting_sampling(circuit, cuts, shots=8192, seed=1, qubits=[0, 4]))

Next
----

- :doc:`notebooks/01_quickstart` is the same with results and figures, and starts the learning path.
- :doc:`user_guide` explains each part.
- :doc:`guide/hpc` runs it on a cluster.
