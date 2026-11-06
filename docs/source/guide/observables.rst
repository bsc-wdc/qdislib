Observables
===========

The notebook :doc:`../notebooks/05_observables` runs everything on this page.

``observable`` takes a per-qubit Pauli string over ``X``, ``Y``, ``Z`` and ``I``, with the first
letter for qubit 0. ``X`` and ``Y`` are measured by rotating into the Z basis automatically, on any
backend. It also takes a weighted sum of Pauli strings, as a list of ``(coefficient, pauli)`` pairs or a
``{pauli: coefficient}`` mapping:

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(4)
    circuit.ry(0.6, range(4))
    for q in range(3):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(4))
    cut = ["CZ_6"]

    print(qd.gate_cutting(circuit, cut, "XZIY", shots=None), qd.analytical_solution(circuit, "XZIY"))

    hamiltonian = [(0.5, "ZZII"), (1.0, "XIII"), (-0.3, "IIZI")]
    print(qd.gate_cutting(circuit, cut, hamiltonian, shots=None))
    print(qd.gate_cutting(circuit, cut, {"ZZII": 0.5, "XIII": 1.0, "IIZI": -0.3}, shots=None))
    print(qd.analytical_solution(circuit, hamiltonian))

The same forms work with every cut type, its ``_subcircuits`` form, and ``qd.analytical_solution``.

Terms run by group
------------------

The terms are not run one by one. They are put in **qubit-wise commuting groups**, where on every qubit
the terms of a group have the same Pauli or ``I``. Each group runs one set of fragment circuits, in its
basis, and every term of the group is read from it: above, ``ZZII`` and ``IIZI`` share one run and
``XIII`` has its own. An identity term is a constant and runs nothing. With many commuting terms, as in a
Z-diagonal Hamiltonian such as QAOA's, the whole sum is one run. The three-step lists have one part per
group.

The changelog gives an example of 7 terms in 3 groups plus an identity, which runs 36 fragment circuits
instead of 84.

.. code-block:: python

    print(len(qd.gate_cutting_subcircuits(circuit, cut, "ZZII")))
    print(len(qd.gate_cutting_subcircuits(circuit, cut, hamiltonian)))

Reading other observables from the same runs
--------------------------------------------

Results of ``sample_subcircuit`` run with shots keep their counts, so ``subs.reconstruct(results,
observable)`` reads any Pauli string in the list's basis from them. An exact run keeps only what the
list's observable needs, since all of a 30-qubit fragment's correlators would be ``2**30`` numbers.

.. code-block:: python

    subs = qd.gate_cutting_subcircuits(circuit, cut, "ZZZZ")
    results = [qd.sample_subcircuit(s, shots=20_000, seed=i) for i, s in enumerate(subs)]
    print(subs.reconstruct(results), subs.reconstruct(results, "ZIIZ"))
