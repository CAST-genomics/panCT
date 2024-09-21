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

Additional examples
~~~~~~~~~~~~~~~~~~~

Below are additional examples based on the HPRC .gbz format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  # Run on a single region
  panct complexity --gbz hprc-v1.1-mc-grch38.gbz \
    --region chr11:119077050-119178859 --out test.tab \
    --metrics sequniq-normwalk,sequniq-normnode

  # Run on a file with a list of regions
  panct complexity --gbz hprc-v1.1-mc-grch38.gbz \
    --region regions.bed --out test.tab \
    --metrics sequniq-normwalk,sequniq-normnode

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: complexity
