.. _commands-complexity:


complexity
==========

Compute the complexity of a pangenome graph.

The ``complexity`` command outputs a file with complexity metrics for an entire graph or for a specified set of regions from a graph.

If a GFA file is provided, the whole graph is processed.

If a GBZ file is provided, you must specify a region or list of regions (as a BED file).

Formulas
~~~~~~~~
The complexity of a region is computed according to the following.

.. math::

  \sum_n \frac{|n|*p_n*(1-p_n)}L

For a node :math:`n` in the region, :math:`|n|` represents the length (in base pairs) of the sequence. So a node with sequence 'ATGAC' would have :math:`|n|=5`, for example.

The percent of sample haplotypes (aka "walks") that visit node :math:`n` is represented by :math:`p_n`.

:math:`L` can be computed in one of two ways:

1. If ``--metrics sequniq-normwalk`` is specified, :math:`L` is computed as the average length of all walks in the region
2. If ``--metrics sequniq-normnode`` is specified, :math:`L` instead represents the average length of all nodes in the region

Usage
~~~~~
.. code-block:: bash

  panct complexity \
    --region REGION or PATH \
    --metrics sequniq-normwalk,sequniq-normnode \
    --reference REFERENCE_ID \
    --out PATH \
    --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
    GFAFILE

Examples
~~~~~~~~
.. code-block:: bash

  panct complexity --out test.tab tests/data/basic.gfa

All files used in these examples are described :doc:`here </project_info/example_files>`.

Additional examples
~~~~~~~~~~~~~~~~~~~

Below are additional examples based on the HPRC .gbz format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  # Run on a single region
  panct complexity \
    --region chr11:119077050-119178859 --out test.tab \
    --metrics sequniq-normwalk,sequniq-normnode \
    hprc-v1.1-mc-grch38.gbz

  # Run on a file with a list of regions
  panct complexity \
    --region regions.bed --out test.tab \
    --metrics sequniq-normwalk,sequniq-normnode \
    hprc-v1.1-mc-grch38.gbz

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: complexity
