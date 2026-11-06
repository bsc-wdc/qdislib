.. Landing page for the poster QR code: https://qdislib.readthedocs.io/en/latest/try.html
   This path is printed on the poster: do not move this file.

Try Qdislib
===========

The companion to the poster *Qdislib: Implementation and Evaluation of Circuit Cutting and Sampling Techniques for Quantum Circuits*. Pick a cutting method and a number of cuts, and see what it costs and that the result comes back. Nothing to install.

.. To add the demo video: upload it (e.g. unlisted on YouTube), then remove the ".. " comment
   markers from the block below and replace VIDEO_ID.

   .. raw:: html

      <div style="position:relative;padding-bottom:56.25%;height:0;overflow:hidden;margin:1em 0;">
        <iframe src="https://www.youtube-nocookie.com/embed/VIDEO_ID" title="Qdislib demo"
                style="position:absolute;top:0;left:0;width:100%;height:100%;border:0;"
                allow="accelerometer; encrypted-media; picture-in-picture" allowfullscreen></iframe>
      </div>

.. raw:: html

   <link rel="stylesheet" href="_static/try.css">
   <div id="qd-try" class="qd-try"><p>This demo needs JavaScript. The notebook below has the same results.</p></div>
   <script src="_static/try_data.js"></script>
   <script src="_static/try.js"></script>

.. only:: latex

   The interactive demo is at https://qdislib.readthedocs.io/en/latest/try.html.

Run the real code
-----------------

The demo above shows stored results. The same four steps run live in a notebook, with a free, temporary environment and nothing to install on your device.

.. raw:: html

   <div class="qd-links">
     <a href="https://mybinder.org/v2/gh/bsc-wdc/qdislib/main?urlpath=lab%2Ftree%2Fnotebooks%2FQdislib_try_it.ipynb"><b>Binder, with real PyCOMPSs</b><span>No account. The subcircuits run as parallel tasks, and you see the task graph. Starting can take 1&ndash;3 minutes, or longer if the service is busy.</span></a>
     <a href="https://colab.research.google.com/github/bsc-wdc/qdislib/blob/main/notebooks/Qdislib_try_it.ipynb"><b>Colab</b><span>Needs a Google account. Its first cell installs Qdislib, and the subcircuits run one after another, without PyCOMPSs.</span></a>
     <a href="notebooks/Qdislib_try_it.html"><b>Just read it</b><span>The notebook with the results of a real run under COMPSs.</span></a>
   </div>

.. only:: latex

   :doc:`notebooks/Qdislib_try_it` is the notebook, with its results.

The poster's four steps
-----------------------

.. container:: qd-wrap

   .. code-block:: python

      from qiskit import QuantumCircuit
      import Qdislib as qd

      circuit = QuantumCircuit(3)
      circuit.ry(0.8, 0)
      circuit.cz(0, 1)
      circuit.ry(0.5, 1)
      circuit.cz(1, 2)

      # 1. Identify gates to cut
      cut = qd.find_cut(circuit, gate_cut=False, wire_cut=True)

      # 2. Generate subcircuits
      subcircuits = qd.wire_cutting_subcircuits(circuit, cut)

      # 3. Execute fragments
      values = [qd.expectation_value(s) for s in subcircuits]

      # 4. Reconstruct expectation value
      result = subcircuits.reconstruct(values)

Install it yourself
-------------------

.. code-block:: bash

   pip install "qdislib[qiskit]"

Qdislib runs serially as it is. For parallel runs add `PyCOMPSs <https://compss-doc.readthedocs.io/en/stable/Sections/01_Installation_Configuration.html>`__: see the :doc:`Quickstart <quickstart>`.

More
----

* Documentation: :doc:`the notebooks <examples>`, the :doc:`user guide <user_guide>`, and :doc:`PyCOMPSs and HPC <guide/hpc>`
* Source code: `github.com/bsc-wdc/qdislib <https://github.com/bsc-wdc/qdislib>`__
* Papers: `arXiv:2505.01184 <https://arxiv.org/abs/2505.01184>`__ and `SC '25 Workshops <https://doi.org/10.1145/3731599.3767547>`__
* Contact: `support-compss@bsc.es <mailto:support-compss@bsc.es>`__
