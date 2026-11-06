Circuits with free parameters
=============================

The notebook :doc:`../notebooks/12_parameters` runs everything on this page.

When the same circuit runs for many sets of angles, such as the tiles of an image or the steps of a variational
loop, cut it **once** with symbolic Qiskit ``Parameter``\ s and give the values afterwards. Planning the cut is done
once, and each set of values runs its own circuits.

.. code-block:: python

    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter
    import Qdislib as qd

    theta, phi = Parameter("theta"), Parameter("phi")
    circuit = QuantumCircuit(4)
    circuit.ry(theta, range(4))
    circuit.cz(0, 1)
    circuit.cz(1, 2)
    circuit.cz(2, 3)
    circuit.ry(phi, range(4))
    cuts = {"gate": ["CZ_6"]}

    # all-in-one: one result per set of values, in a list
    values = qd.circuit_cutting(circuit, cuts, "ZZZZ", shots=4096, seed=1,
                                parameter_values=[[0.1, 0.2], [0.3, 0.4]])
    print(values)

    # three-step: bind the list, then run and reconstruct as usual
    subs = qd.circuit_cutting_subcircuits(circuit, cuts, "ZZZZ")
    print(subs.parameters)
    bound = subs.bind({theta: 0.1, phi: 0.2})     # or [0.2, 0.1], in subs.parameters order
    print(bound.reconstruct([qd.expectation_value(s, shots=4096) for s in bound]))

``parameter_values=`` works on every cutting function (gate, wire, Hadamard, ``circuit_cutting`` and their
``_sampling`` forms). A set of values is a list in the order of ``circuit.parameters``, which is sorted by name, or a dict. It combines with ``sync=False`` (one future per set),
``circuits_per_task=`` and ``error=``. Each set gets its own seed stream from ``seed``. The values are bound where each circuit
is built, so a fragment's template still travels to a worker once.

A **cut gate's** angle must be a number, since it sets the cut's terms and weights. An unbound list still gives its circuits,
``subs.circuits("qiskit")``, with their parameters, to bind yourself.
