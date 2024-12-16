.. _formats-walks:

Walks
=====

A ``.walk`` file stores the walks (W lines) of a pangenome graph.

Unlike walks stored in the GFA format, the walks in the ``.walk`` format are designed to support faster querying from their node IDs.

Overview
~~~~~~~~
# TODO

Examples
~~~~~~~~
You can find an example of a ``.walk`` file without any extra fields in `tests/data/basic.walk <https://github.com/cast-genomics/panct/blob/main/tests/data/basic.walk>`_:

.. include:: ../../tests/data/basic.walk
   :literal:

And here's the corresponding GFA file:

.. include:: ../../tests/data/basic.gfa
   :literal:

Compressing and indexing
~~~~~~~~~~~~~~~~~~~~~~~~
If it isn't already, we encourage you to bgzip compress and index your ``.walk`` file whenever possible. This will reduce both disk usage and the time required to parse the file, but it is entirely optional. You can use the ``bgzip`` and ``tabix`` commands.

.. code-block:: bash

  bgzip file.walk
  tabix file.walk.gz

