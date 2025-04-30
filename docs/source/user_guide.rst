User guide
==========

This guide covers the basics of using Qdislib, as well as details on the main algorithms included in the library for circuit partitioning, simulation, and distributed execution.

Qdislib supports *gate cutting* and *wire cutting* techniques to decompose large quantum circuits and execute them on distributed and heterogeneous environments. All computations are parallelized transparently using `PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`_.

Example
-------

Below is a basic example demonstrating how to construct a quantum circuit and apply gate cutting using Qdislib:

.. code-block:: python

    from qibo import models, gates
    import qd as qd  # qd is Qdislib

    # Build a 10-qubit circuit
    def entire_circuit():
        c = models.Circuit(10)
        c.add(gates.H(0))
        c.add(gates.RX(1, 0.3))
        c.add(gates.RY(2, 0.4))
        c.add(gates.RZ(3, 0.5))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.CZ(2, 3))
        c.add(gates.CZ(3, 4))
        c.add(gates.CZ(4, 5))
        c.add(gates.CZ(5, 6))
        return c

    circuit = entire_circuit()

    # Find a gate to cut
    cut = qd.find_cut(circuit)

    # Apply gate cutting
    result = qd.gate_cutting(circuit, cut)
    print("Reconstructed expectation value:", result)

It is worth noting that, although the code above looks completely sequential,
all Qdislib algorithms and operations are internally parallelized using PyCOMPSs.

How to run Qdislib
------------------

Qdislib can be installed and used as a regular Python library. However,
because it relies on PyCOMPSs to parallelize its tasks, you must execute Qdislib applications
with the `runcompss` or `enqueue_compss` command-line tools provided by PyCOMPSs:

.. code-block:: bash

    runcompss my_Qdislib_application.py

Alternatively, inside Jupyter notebooks, you can use the `ipycompss` interface:

.. code-block:: python

    import ipycompss
    ipycompss.start(project="project.xml", resources="resources.xml", monitor=True, graph=True)

Refer to the :doc:`quickstart guide <quickstart>` and :doc:`examples <examples>` for details on PyCOMPSs setup and configuration.

Algorithms
----------

Qdislib includes two core circuit cutting strategies:

- **Gate Cutting**: Cuts selected two-qubit gates to divide a circuit into smaller, independently executable subcircuits. These are then simulated and post-processed to reconstruct the expectation value of the original circuit.

- **Wire Cutting**: Cuts "wires" (qubit connections) between gates to disconnect logical qubits. This requires more advanced processing, often involving state preparation and postselection logic.

### Automatic Cut Detection

To simplify the cutting process, Qdislib provides the `find_cut` utility to automatically suggest suitable locations for gate or wire cutting:

.. code-block:: python

    cut = qd.find_cut(circuit)                # Gate cutting suggestion
    wire_cut = qd.find_cut(circuit, gate_cut=True)  # Wire cutting suggestion

The returned list contains candidate gate names or gate pairs that are suitable cut points.

### Subcircuit Generation Only

In addition to executing the full gate or wire cutting workflow, Qdislib also allows users to **only extract the subcircuits** without performing expectation value reconstruction. This enables:

- Custom execution of subcircuits (e.g., submitting to specific quantum devices).
- Saving/loading subcircuits for delayed simulation.
- Manual post-processing outside Qdislib.

You can generate and retrieve subcircuits using:

.. code-block:: python

    # For gate cutting
    val, subcircuits = qd.gate_cutting(circuit, cut, return_subcircuits=True)

    # For wire cutting
    val, subcircuits = qd.wire_cutting(circuit, wire_cut, return_subcircuits=True)

    # Alternatively, use the explicit subcircuit-extraction functions:
    _, subcircuits = qd.gate_cutting_subcircuits(circuit, cut, return_subcircuits=True)
    _, subcircuits = qd.wire_cutting_subcircuits(circuit, wire_cut, return_subcircuits=True)

The returned `subcircuits` list contains all subcomponents resulting from the cutting process. These can be simulated independently using Qiskit, Qibo, or submitted to QPUs.

This modular structure gives users fine-grained control over hybrid quantum-classical workflows.


Resource Allocation
-------------------

Qdislib is designed to run efficiently across CPUs, GPUs, and quantum devices depending on the user’s configuration and available resources. By default, simulation is performed on CPU backends.

.. _gpu-support-label:

Using GPUs
----------

Qdislib includes support for GPU-accelerated execution using `cuQuantum` (with Qiskit) and `Qibojit` (with Qibo). GPU support is available for the following algorithms:

- `gate_cutting`
- `wire_cutting`

To enable GPU execution, set the following environment variable before running your application:

.. code-block:: bash

    export QDISLIB_GPU_AVAILABLE=True

When enabled, Qdislib offloads eligible subcircuit simulations to GPU, while still coordinating the overall computation through CPU. This GPU usage is transparent to the user: memory transfers are handled internally and automatically
