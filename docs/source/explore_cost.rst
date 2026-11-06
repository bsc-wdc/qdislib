Cutting and cost
================

What a cut costs, how Qdislib lowers it, and how it finds one.

.. raw:: html

   <link rel="stylesheet" href="_static/try.css">
   <link rel="stylesheet" href="_static/explore.css">
   <script src="_static/explore_engine.js"></script>
   <script src="_static/explore_data.js"></script>

What lowers the cost
--------------------

Each cut multiplies the work, so Qdislib has several ways to multiply it by less: optimized runs, joint cuts, golden cut points, and cuts that follow the gate's angle. *Stored from* :func:`~Qdislib.cutting_overhead` *and real runs.*

.. raw:: html

   <div id="qdx-low" class="qd-try"><p>This section needs JavaScript.</p></div>

Finding cuts
------------

Pick a circuit and a device size, and see the cut :func:`~Qdislib.find_cut` chooses and what it costs. *Stored from real runs.*

.. raw:: html

   <div id="qdx-find" class="qd-try"><p>This section needs JavaScript.</p></div>

Planning a run
--------------

What a number of cuts costs for each kind of cut, before anything runs: the subcircuits, and the shots for an almost exact result. *Stored from* :func:`~Qdislib.cutting_overhead` *and* :func:`~Qdislib.shots_for`.

.. raw:: html

   <div id="qdx-plan" class="qd-try"><p>This section needs JavaScript.</p></div>

.. raw:: html

   <script src="_static/explore.js"></script>

.. only:: latex

   The interactive page is at https://qdislib.readthedocs.io/en/latest/explore_cost.html.
