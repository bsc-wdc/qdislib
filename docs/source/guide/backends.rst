Backends and SDKs
=================

The notebooks :doc:`../notebooks/04_backends` and :doc:`../notebooks/08_backend_internals` run most of this page.

SDKs
----

A circuit can come from Qibo, Qiskit, CUDA-Q, PennyLane or QiliSDK. Subcircuits can be written in any of them and run on any of them, whatever
the input was. Install what you use:

.. code-block:: bash

    pip install "qdislib[qiskit]"      # or [qibo], [cudaq], [pennylane], [qilisdk], [all]

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - SDK
     - ``backend=``
     - Simulator
   * - Qiskit
     - ``"qiskit"``
     - Aer
   * - Qibo
     - ``"qibo"``
     - Qibo's own
   * - CUDA-Q
     - ``"cudaq"``
     - CUDA-Q's
   * - PennyLane
     - ``"pennylane"``
     - ``default.qubit``, through PennyLane
   * - QiliSDK
     - ``"qilisdk"``
     - QiliSim (``extra={"gpu": True}`` for its GPU mode)

``shots=None`` does not use these. It is exact, on a small numpy statevector simulator inside Qdislib that applies the DAG's gates directly.

``qd.from_qasm(text, software=...)`` builds a circuit in any of them from OpenQASM 2, so any framework that exports it can come in.

Another SDK can be added in one place, ``Qdislib/utils/sdks.py``: a converter to and from the DAG, and a backend, registered with
``Qdislib.utils.sdks.register_sdk``.

Choosing a backend
------------------

``backend=`` is accepted by every cutting call and by ``expectation_value`` and ``sample_subcircuit``. It takes an SDK name (its simulator), a
registered backend name, or a backend object, used as it is. ``None`` is the circuit's own SDK.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(4)
    circuit.ry(0.5, range(4))
    for q in range(3):
        circuit.cz(q, q + 1)
    cut = ["CZ_6"]

    print(qd.available_backends())
    print(qd.gate_cutting(circuit, cut, shots=8192, seed=1, backend="qibo"))

    aer = qd.build_backend("AER_BACKEND", qd.BackendConfig(num_qubits=20))
    print(qd.gate_cutting(circuit, cut, shots=8192, seed=1, backend=aer))

A run that does something other than what was asked says so, once per process, with a ``UserWarning``: CUDA-Q running on ``qpp-cpu`` because its GPU target cannot be
used, a Hadamard cut of a general gate using the Pauli basis because KAK is not available, and ``shots=None`` with a registered backend or backend object, which is an exact Qiskit simulation, not a run on it.

Dense unitary gates
-------------------

A Qiskit ``UnitaryGate`` can stay in a cut circuit as it is. It runs on Qiskit/Aer, Qibo and Qdislib's exact simulator; CUDA-Q and PennyLane refuse it with a clear error. A circuit loaded with multi-qubit isometries, such as an MPS-loaded volume, is simulated with its dense gates instead of their CNOT decomposition: the changelog gives hundreds of gates instead of hundreds of thousands.

.. code-block:: python

    from qiskit.circuit.library import UnitaryGate
    from scipy.stats import unitary_group

    dense = QuantumCircuit(4)
    dense.append(UnitaryGate(unitary_group.rvs(4, random_state=1)), [0, 1])
    dense.cz(1, 2)
    dense.append(UnitaryGate(unitary_group.rvs(4, random_state=2)), [2, 3])

    print(qd.gate_cutting(dense, ["CZ_2"], "ZZZZ", shots=None), qd.analytical_solution(dense, "ZZZZ"))

.. _gpu-support-label:

GPUs
----

Aer's GPU simulator (``device="GPU"``, cuStateVec, which needs ``qiskit-aer-gpu``) and CUDA-Q's ``nvidia`` target run subcircuits on a GPU. Send only the large subcircuits there from the three-step loop.

.. not run: needs a GPU

.. code-block:: python

    gpu = qd.build_backend("AER_BACKEND", qd.BackendConfig(extra={"use_gpu": True}))
    value = qd.gate_cutting(circuit, cut, backend=gpu)

    subs = qd.gate_cutting_subcircuits(circuit, cut)
    values = [
        qd.expectation_value(s, backend=gpu if s.num_qubits >= 18 else None)
        for s in subs
    ]

Under PyCOMPSs these run in tasks that ask for a GPU: see :doc:`hpc`.

A backend of your own
---------------------

A new backend needs one method, ``evaluate_samples``: the final counts of a subcircuit, keyed qubit 0 first. Cutting derives the rest; a mid-circuit cut point is deferred onto an
extra qubit read at the end. For a device that runs Qiskit circuits, ``qd.sample_with_qiskit`` builds the circuit and reads the counts:

.. code-block:: python

    from qiskit_aer import AerSimulator

    class MyDevice(qd.BaseBackend):
        def evaluate_samples(self, dag, num_qubits, shots=1024, seed=None):
            simulator = AerSimulator(seed_simulator=seed)
            return qd.sample_with_qiskit(
                dag, num_qubits, lambda c: simulator.run(c, shots=shots).result().get_counts()
            )

    device = MyDevice("MY_DEVICE", qd.BackendConfig(num_qubits=20))
    print(qd.gate_cutting(circuit, cut, backend=device, seed=1))

These are optional, when they apply:

- ``evaluate_joint_counts`` if the device reports mid-circuit outcomes;
- ``evaluate_joint_counts_many`` and ``packs`` if it runs a list of circuits in one job;
- ``uses_gpu`` if it runs on a GPU;
- ``evaluate_expectation`` for ``qd.expectation_value`` and the 1.x paths;
- ``transpile`` with ``_execute_samples`` and ``_execute_estimator`` for ``QdBackEnd``'s queue.

``qd.register_backend(name, module, class_name)`` makes it available by name and imports it lazily, with a constructor taking a config. Under PyCOMPSs the backend object has to reach a task, so its module must be importable on the workers and the object picklable.

A pool of devices
-----------------

``qd.QdBackEnd(qpu_dict={...})`` holds a pool of backends and routes each whole circuit: it keeps those with enough qubits and tries the shortest queue first, falling back to the next if one fails. It is independent of cutting, where you choose the device per subcircuit with ``s.num_qubits`` in the three-step loop.

.. code-block:: python

    pool = qd.QdBackEnd(qpu_dict={"AER_BACKEND": {"num_qubits": 30}})
    bell = QuantumCircuit(2)
    bell.h(0)
    bell.cx(0, 1)
    bell.measure_all()
    print(pool.sample_circuit(bell, shots=1024).counts)
