User guide
==========

This guide covers the basics of using Qdislib, and details on the different
algorithms included in the library.

Qdislib provides two major programming interfaces: an API to manage
quantum circuits in a distributed way, and an graph-based interface to
work with different circuit cutting and mapping algorithms.

The typical workflow in Qdislib consists of the following steps:

 1. Create a quantum circuit
 2. Convert to graph
 3. Apply circuit cutting
 4. Get result from the circuit executin in distributed environments

An example is as follows:

.. code:: python

    import Qdislib
    from Qdislib.api import circuit_cutting

    if __name__ == '__main__':
        # step 1
        qc = Qdislib.new_circuit()

        # step 2
        graph = Qdislib.convert(qc)

        # step 3
        subgraphs = circuit_cutting(graph, other_parameters)

        # step 4
        result = Qdislib.execute(subgraphs)

It is worth noting that, although the code above looks completely sequential,
all Qdislib algorithms and operations are paralellized using `PyCOMPSs
<https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`_.

How to run Qdislib
------------------

Qdislib can be installed and used as a regular Python library. However,
Qdislib makes use of PyCOMPSs directives internally to parallelize all the
computation. This means that applications using Qdislib need to be executed
with PyCOMPSs. This can be done with the ``runcompss`` or
``enqueue_compss`` commands:

.. code:: bash

    runcompss my_Qdislib_application.py

For more information on how to start running your Qdislib applications, refer
to the :doc:`quickstart guide <quickstart>`.


Resource allocation
-------------------

All Qdislib tasks are allocated a specific number of computational resources.
By default, each task receives one CPU. This number can be adjusted according
to the specific needs of the program by setting the environment
variable ComputingUnits before executing the script:

``export ComputingUnits=8``

The above example sets the number of CPUs available for each task to 8. This
is specifically useful for algorithms (for example those implemented in NumPy)
that automatically take advantage of fine-grained parallelism facilitated by a
higher number of computing units.

.. _gpu-support-label:

Using GPUs with CuPy
--------------------

Qdislib includes support for GPU using CuPy (CUDA) in the following algorithms:
 - circuit_cutting

In order to enable Qdislib to use GPUs for this algorithms the user must 
set the environment variable QDISLIB_GPU_AVAILABLE before executing the script:

``export QDISLIB_GPU_AVAILABLE=True``

The above example makes Qdislib use GPU for any tasks that can use CuPy acceleration.
Qdislib's implementation with cupy is completely transparent for the end user because
all the data is stored in the main memory and only transfered to GPU when necessary.
Also, for the rest of the algorithms that are not listed above they will just use usual
CPU computation. In this way all the resources will be used efficiently.

The required external library for this functionality is cupy >= 0.9.6
This functionality only works on computers with CUDA enabled devices.
