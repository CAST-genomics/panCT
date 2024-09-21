.. _commands-complexity:


complexity
==========

Compute the complexity of a pangenome graph.

The ``complexity`` command outputs a file with complexity metrics for an entire graph or for a specified set of regions from a graph.

If a GFA file is provided, the whole graph is processed.

If a GBZ file is provided, you must specify a region or list of regions to process.

Usage
~~~~~
.. code-block:: bash

  panct complexity \
    --gfa GFAFILE \
    --out OUTFILE
    --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET]

Examples
~~~~~~~~
.. code-block:: bash

  panct complexity --gfa tests/data/test.gfa --out test.tab

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: complexity
