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

    # A small example circuit (works the same with Qiskit, CUDA-Q or PennyLane circuits).
    circuit = models.Circuit(2)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.RY(0, 0.8))

    # 1. Cut the circuit into independently executable subcircuits.
    #    Use find_cut to discover good cut points automatically.
    subcircuits = qd.gate_cutting_subcircuits(circuit, ["CZ_2"], software="qibo")

    # 2. Evaluate each subcircuit — this loop is yours to control.
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

    # 3. Recombine into the original circuit's expectation value.
    result = subcircuits.reconstruct(values)
    print(result)  # ≈ -0.717

``qd.wire_cutting_subcircuits(circuit, [("H_1", "CZ_2")], software="qibo")`` works
the same way for **wire cutting**. The ``software`` argument (``"qibo"``,
``"qiskit"`` or ``"cudaq"``) selects the backend the subcircuits are emitted for,
and must match the backend you evaluate them on.

All-in-one
----------

If you'd rather not manage the loop, ``qd.gate_cutting(...)`` and
``qd.wire_cutting(...)`` run the cut, execution and reconstruction in a single
call — convenient for quick experiments. The three-step workflow above is the
flexible, production-oriented path.

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

    val = qd.gate_cutting_hadamard(circuit, ["ISWAP_3"], "ZZ", software="qibo")

    # three-step form, and a Monte-Carlo variant for many cuts:
    subs = qd.gate_cutting_hadamard_subcircuits(circuit, ["ISWAP_3"], "ZZ")
    val = subs.reconstruct([qd.expectation_value(s) for s in subs])
    val = qd.gate_cutting_hadamard_sampling(circuit, ["ISWAP_3"], "ZZ", num_samples=5000)

Mixing gate and wire cutting
----------------------------

``gate_wire_cutting`` cuts some gates **and** some wires on the same circuit in a
single, fully parallel workflow:

.. code-block:: python

    val = qd.gate_wire_cutting(
        circuit,
        gates_cut=["CZ_2"],
        wire_cut=[("RY_3", "CZ_4")],
        observable="ZZZ",
        software="qibo",
    )

See the :doc:`User guide <user_guide>` for the full API, GPU acceleration, and
distributed execution on PyCOMPSs.
