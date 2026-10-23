Quickstart
==========

Installation
------------

Qdislib needs at least one quantum SDK backend. Install it together with the
backend you want to use (the SDKs are optional extras):

.. code-block:: bash

    pip install "qdislib[qibo]"      # or [qiskit], [cudaq], or [all]

The library imports fine with no SDK installed (useful for documentation
builds), but you need one of ``qibo``, ``qiskit`` or ``cudaq`` to run circuits.

Circuit cutting in three steps
-------------------------------

The recommended way to use Qdislib is the **three-step workflow**: generate the
subcircuits, evaluate each one in your own loop, then recombine the results.
Owning the evaluation loop is the point — that is where you dispatch to a
simulator, a GPU, or a real QPU — while Qdislib handles the cutting and the
reconstruction.

.. code-block:: python

    import Qdislib as qd
    from qibo import models, gates

    # A small example circuit (works the same with Qiskit, CUDA-Q, PennyLane or QiliSDK circuits).
    circuit = models.Circuit(2)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.RY(0, 0.8))

    # 1. Cut the circuit into independently executable subcircuits.
    #    Use find_cut to discover good cut points automatically.
    subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"])

    # 2. Evaluate each subcircuit — this loop is yours to control.
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

    # 3. Recombine into the original circuit's expectation value.
    result = subcircuits.reconstruct(values)
    print(result)  # ≈ -0.717

``qd.wire_cutting_subcircuits(circuit, [("H_1", "CZ_2")])`` works the same way for
**wire cutting**, and so do the Hadamard cut and ``circuit_cutting`` below. Each
subcircuit ``s`` is written in the circuit's own SDK unless you pass ``software=``
(``"qibo"``, ``"qiskit"``, ``"cudaq"``, ``"pennylane"``, ``"qilisdk"``); ``s.circuit`` is the native
circuit. ``qd.expectation_value(s, shots=..., backend=...)`` runs it anywhere: the
SDK's simulator by default, or a QPU backend.

All-in-one
----------

If you'd rather not manage the loop, ``qd.gate_cutting(circuit, cuts, observable)``
(and ``wire_cutting``, ``gate_cutting_hadamard``, ``circuit_cutting``) run the cut,
execution and reconstruction in a single call. They all take the same keywords:
``backend`` (where to run), ``shots`` (``None`` is exact), ``seed`` and ``sync``.
Each also has a ``_sampling`` form that returns the output distribution. The
three-step workflow above is the flexible, production-oriented path.

Cutting any gate (Hadamard test)
--------------------------------

``gate_cutting`` cuts any two-qubit gate, but a non-``CZ``/``CX`` gate becomes
several ``CZ`` cuts. The ancilla-based **Hadamard-test** cut (Harrow & Lowe,
`arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) cuts *any* gate as a single
cut at the optimal overhead, with no mid-circuit measurement:

.. code-block:: python

    circuit = models.Circuit(2)
    circuit.add(gates.RY(0, 0.9))
    circuit.add(gates.RY(1, 0.4))
    circuit.add(gates.iSWAP(0, 1))

    val = qd.gate_cutting_hadamard(circuit, ["ISWAP_3"], "ZZ", shots=None)  # exact

    # three-step form, and Monte-Carlo sampling of the configurations for many cuts:
    subs = qd.gate_cutting_hadamard_subcircuits(circuit, ["ISWAP_3"], "ZZ")
    val = subs.reconstruct([qd.expectation_value(s, shots=4096) for s in subs])
    val = qd.gate_cutting_hadamard(circuit, ["ISWAP_3"], "ZZ", num_samples=5000)

Mixing any cuts
---------------

``circuit_cutting`` cuts gates (measure-and-prepare or Hadamard) **and** wires on the
same circuit, in any number, in a single fully parallel workflow:

.. code-block:: python

    cuts = {"gate": ["CZ_2"], "wire": [("RY_3", "CZ_4")]}
    val = qd.circuit_cutting(circuit, cuts, "ZZZ", optimized=True)
    qd.cutting_overhead(circuit, cuts, optimized=True)  # circuits and gamma, before running

See the :doc:`User guide <user_guide>` for the full API, GPU acceleration, and
distributed execution on PyCOMPSs.
