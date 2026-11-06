.. _error-bars-label:

Shots and error bars
====================

The notebook :doc:`../notebooks/11_shots_and_errors` runs everything on this page.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(6)
    circuit.ry(0.5, range(6))
    for q in range(5):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(6))
    cut = ["CZ_9"]
    hamiltonian = [(0.5, "Z" * 6), (1.0, "X" + "I" * 5)]

Error bars: ``error=``
----------------------

``error=True`` returns the value with its **standard error**, by bootstrap over the shots: every fragment
circuit's counts are resampled, and the ``num_samples`` draws within their strata, and the cheap
reconstruction is repeated, so it costs reconstruction time, not device time. The value is an
:class:`~Qdislib.Estimate`, a ``float`` that also carries ``.error`` and ``.shots``, so code that treats it as
a number keeps working.

.. code-block:: python

    value = qd.circuit_cutting(circuit, cut, hamiltonian, shots=4000, seed=1, error=True)
    print(value, value.error, value.shots)
    print(qd.gate_cutting(circuit, cut, "Z" * 6, shots=4000, seed=1, error=500))   # 500 resamples (default 200)

The ``_sampling`` functions take it too, see :doc:`distributions`.

According to the changelog it matches the spread of repeated runs; with ``num_samples`` and more configurations than draws it errs high
(up to about 40% in tests), as it resamples the draws and their shots separately. An exact run has ``error``
0.

One budget: ``total_shots=``
----------------------------

``total_shots=`` spends one budget across all the fragment circuits instead of ``shots`` per circuit, in two
rounds: a pilot, 20% of the budget, split by each circuit's weight in the reconstruction, then the rest by
weight times the spread the pilot measured (Neyman allocation). Identical circuits share theirs. In tests it reached the accuracy of an even split with 1.4-1.9 times fewer shots. Under
PyCOMPSs the master waits for the pilot round to decide the second. With ``num_samples`` it sets the shots per
draw instead.

.. code-block:: python

    value = qd.circuit_cutting(circuit, cut, hamiltonian, total_shots=10**6, seed=1, error=True)
    print(value)

Drawing configurations: ``num_samples=``
----------------------------------------

A cut is a weighted sum over configurations: ``6**n`` for ``n`` gate cuts, ``8**n`` for ``n`` wire cuts, products of
the Hadamard settings. By default every configuration runs. With ``num_samples=N``, Qdislib instead uses ``N``
configurations in proportion to ``|weight| / gamma`` and runs only the fragment circuits they use, each with
``shots`` per draw that uses it. This is the quasiprobability Monte Carlo estimate, and it is unbiased.

The draws are **allocated**, not drawn at random (Dai and Koczor: a share in proportion to the weight is never
worse than random draws). When all configurations fit in ``N`` each gets its share, which is the same as running
them all; otherwise the cuts that fit are allocated and the rest drawn.

The error has two parts. Which configurations are drawn contributes at most ``gamma / sqrt(N)``, **whatever the
shots**; the shot noise is the rest. When configurations do not all fit in ``N`` the first part usually dominates,
so more shots per draw do not help: for a given budget, use **more draws with fewer shots each**. ``shots=1`` is
fine. When they all fit, the error falls with the shots, as when running every configuration.

.. code-block:: python

    value = qd.gate_cutting(circuit, cut, "Z" * 6, shots=1, num_samples=50_000, seed=1)
    print(value)
    dist = qd.wire_cutting_sampling(circuit, [("CZ_9", "CZ_10")], shots=1, num_samples=50_000, seed=1)

    # three-step: each subcircuit's draws says how many draws use it
    subs = qd.circuit_cutting_subcircuits(circuit, {"gate": cut}, "Z" * 6, num_samples=5_000, seed=1)
    print(subs.reconstruct([qd.expectation_value(s, shots=s.draws) for s in subs]))

Use it when there are too many configurations to list, such as a dozen gate cuts (``6**12``) or a mix of gate and
wire cuts, since the cost is set by ``N``, not by the number of cuts; and when the weights are unequal, as in Hadamard cuts
of rotations, so the shots go where the weight is.

It does not change ``gamma``. With drawn configurations the shots for a given precision grow as ``gamma**2`` at
least, exactly ``(gamma**2 - <O>**2) / (1 - <O>**2)`` times the uncut circuit's, while running every configuration with its share of the
shots can need fewer.

The worst case: ``shots_for``
-----------------------------

``qd.shots_for(circuit, cuts, epsilon, delta=0.05)`` gives the shots that guarantee an error below ``epsilon`` with
probability ``1 - delta`` for *any* circuit with this cut (Hoeffding: ``2 gamma**2 / epsilon**2 * ln(2/delta)``). It is
a worst case, usually far more than ``error=True`` shows a run needs.

.. code-block:: python

    print(qd.shots_for(circuit, cut, 0.01))
