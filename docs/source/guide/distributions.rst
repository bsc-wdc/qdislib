Output distributions
====================

The notebook :doc:`../notebooks/06_sampling` runs everything on this page.

The ``_sampling`` form of every cutting function returns the output distribution
``{bitstring: probability}``, with qubit 0 first in each key. It has the same all-in-one and three-step
forms as the expectation value.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(5)
    circuit.ry(0.6, range(5))
    for q in range(4):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(5))
    cut = ["CZ_7"]

    dist = qd.gate_cutting_sampling(circuit, cut, shots=100_000, seed=0)
    dist = qd.wire_cutting_sampling(circuit, [("CZ_7", "CZ_8")], shots=100_000, seed=0)
    dist = qd.gate_cutting_hadamard_sampling(circuit, cut, shots=100_000, seed=0)
    dist = qd.circuit_cutting_sampling(circuit, {"gate": cut}, shots=100_000, seed=0, optimized=True)

    # three-step: sample each subcircuit yourself, then reconstruct
    subs = qd.gate_cutting_hadamard_sampling_subcircuits(circuit, cut)
    dist = subs.reconstruct([qd.sample_subcircuit(s, shots=8192) for s in subs])
    print(sorted(dist.items(), key=lambda kv: -kv[1])[:3])

With the Hadamard cut every fragment is read in the Z basis anyway, so the distribution costs the same
circuits as one expectation value. ``shots=None`` gives the exact distribution on every cut type. A
distribution carries a quasiprobability overhead, and needs more shots than an expectation value.

Marginals: ``qubits=``
----------------------

When only some qubits matter, name them. Bit ``j`` of each key is qubit ``qubits[j]``. Each fragment
circuit then returns only what the marginal needs, so the cost and the result grow with the kept qubits, not
with the circuit. The changelog's example is a 28-qubit circuit read on 13 qubits, which returns
``2**13`` outcomes instead of ``2**28``.

.. code-block:: python

    print(qd.circuit_cutting_sampling(circuit, cut, shots=20_000, seed=1, qubits=[4, 0]))

More than 24 qubits: ``FactoredDistribution``
---------------------------------------------

A distribution has ``2**n`` outcomes, but a cut circuit's is a contraction of one small tensor per
fragment. Over more than 24 output qubits ``circuit_cutting_sampling`` returns a
:class:`~Qdislib.FactoredDistribution`, which keeps those tensors and contracts them when read. No circuit
runs again. Up to 24 qubits the output is the plain dict.

.. code-block:: python

    wide = QuantumCircuit(40)
    wide.ry(0.6, range(40))
    for q in range(39):
        wide.cz(q, q + 1)
    wide.rx(0.3, range(40))
    wide_cuts = {"gate": ["CZ_53", "CZ_66"]}

    d = qd.circuit_cutting_sampling(wide, wide_cuts, shots=None)
    print(d)
    print(d["0" * 40])         # one outcome (with shots, a quasi-probability)
    print(d.marginal([3, 17])) # a marginal, as a dict
    print(d.most_likely(2))    # the likeliest outcomes (a beam search)
    print(d.sample(3, seed=1)) # {bitstring: count}

When you know which qubits you read, ``qubits=`` is cheaper still, as each fragment circuit keeps only
those. ``error=`` and ``num_samples=`` build every outcome, so above 24 qubits they need ``qubits=`` and
refuse before anything runs otherwise.

Error bars
----------

``error=True`` returns a :class:`~Qdislib.Distribution`, the same dict with ``.errors`` per outcome and
``.shots``. See :doc:`shots`.

.. code-block:: python

    d = qd.circuit_cutting_sampling(circuit, cut, shots=20_000, seed=1, error=True, qubits=[0, 4])
    print(d, d.errors["00"])
