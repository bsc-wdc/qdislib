Quickstart
==========

Installation
------------

Qdislib needs Python 3.10 or later and at least one quantum SDK, installed as an extra. With ``pip`` or `uv <https://docs.astral.sh/uv/>`__, which installs the same package and is typically much faster:

.. code-block:: bash

    pip install "qdislib[qiskit]"            # or [qibo], [cudaq], [pennylane], [qilisdk]
    uv pip install "qdislib[qiskit]"         # in an environment
    uv add "qdislib[qiskit]"                 # in a uv project
    uv run --with "qdislib[qiskit]" app.py   # one-off

A bare ``pip install qdislib`` has no SDK, so it can import Qdislib but cannot run a circuit. ``[all]`` is all the SDKs plus the semantic cache (``[cache]``). It is about 3.3 GB, mostly NVIDIA libraries for CUDA-Q, so most people do not want it. ``[iqm]`` for IQM hardware and ``[qpu]`` for the BSC QPU listener are separate.

.. note::

   ``pymetis``, one of ``find_cut``'s strategies, has no wheel for Linux arm64 or Windows, so it is not installed there. ``find_cut`` then warns and runs without METIS, and a wire cut can come out larger. ``pip install "qdislib[metis]"`` adds it, and needs a C compiler. Windows has been checked for wheel availability only, not run.

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
