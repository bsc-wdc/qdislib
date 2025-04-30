.. Qdislib documentation master file, created by
   sphinx-quickstart on Tue Mar  5 09:41:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to Qdislib!
===================

.. image:: ../logo.png
   :alt: Qdislib Logo
   :align: left
   :scale: 100%

The Quantum Distributed Computing Library (Qdislib) provides distributed
Quantum algorithms ready to use as a library.
So far, Qdislib is highly focused on Quantum circuit execution in High Performance Computers (HPCs) and Quantum Machines.
However, other types of Quantum algorithms might be added in the future.
The main objective of Qdislib is to facilitate the execution of large Quantum circuits in distributed platforms, such as clusters, clouds, and supercomputers.

Qdislib has been implemented on top of
`PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`_ programming model,
and it is being developed by the
`Workflows and Distributed Computing <https://www.bsc.es/discover-bsc/organisation/scientific-structure/workflows-and-distributed-computing>`_
group of the `Barcelona Supercomputing Center <http://www.bsc.es>`_.

Documentation
-------------

* :doc:`Quickstart <quickstart>`
* :doc:`User guide <user_guide>`
* :doc:`API Reference <api_reference>`
* :doc:`Main <main>`
* :doc:`Modules <modules>`
* :doc:`Examples <examples>`
* :doc:`Development <development>`
* :doc:`Glossary <glossary>`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api_reference
   development
   examples
   glossary
   main
   modules
   quickstart
   user_guide

Source code
-----------

The source code of Qdislib is available online at `Github <https://github.com/bsc-wdc/Qdislib>`_.

Support
-------

If you have questions or issues about the Qdislib you can join us in `Slack <https://bit.ly/bsc-wdc-community>`_.

Alternatively, you can send us an e-mail to `support-compss@bsc.es <mailto:support-compss@bsc.es>`_.

Citing dislib
-------------

If you use Qdislib in a scientific publication, we would appreciate citations to the following paper:

\M. Tejedor, B. Casas, J. Conejero, A. Cervera-Lierta and R. M. Badia, "Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing" in *In Proceedings of The International Conference for High Performance Computing, Networking, Storage, and Analysis (SC'25)*, 2025, pp. 1-11

Bibtex
......

.. code:: latex

   @inproceedings{Qdislib,
               title       = {{Distributed Quantum Circuit Cutting for Hybrid Quantum-Classical High-Performance Computing}},
               author      = {Mar Tejedor and Berta Cervera and Javier Conejero and Alba Cervera-Lierta and Rosa M. Badia},
               booktitle   = {In Proceedings of The International Conference for High Performance Computing, Networking, Storage, and Analysis(SC'25).},
               pages       = {1-11},
               year        = {2025},
    }


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
