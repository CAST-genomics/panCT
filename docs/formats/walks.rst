.. _formats-walks:


Walks
=====

A ``.walk`` file stores the walks (W lines) of a pangenome graph.

Compared to the way that walks are stored in the GFA format, this format is designed to support faster querying of the walks from their node IDs.

Overview
~~~~~~~~
# TODO

Examples
~~~~~~~~
You can find an example of a ``.walk`` file without any extra fields in `tests/data/basic.walk <https://github.com/cast-genomics/panct/blob/main/tests/data/basic.walk>`_:

.. include:: ../../tests/data/test.walk
   :literal:

Compressing and indexing
~~~~~~~~~~~~~~~~~~~~~~~~
If it isn't already, we encourage you to bgzip compress and index your ``.walk`` file whenever possible. This will reduce both disk usage and the time required to parse the file, but it is entirely optional. You can either use the :doc:`index command </commands/index>` or the ``bgzip`` and ``tabix`` commands.

.. code-block:: bash

  bgzip file.walk
  tabix -s 1 -b 2 -e 3 file.walk.gz
