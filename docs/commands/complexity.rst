.. _commands-complexity:


complexity
==========

Compute the complexity of a pangenome graph.

The ``complexity`` command outputs a file with complexity metrics for an entire graph or for a specified set of regions from a graph.

If a GFA file is provided, the whole graph is processed.

If a GBZ file is provided, you must specify a region or list of regions (as a BED file).

..
  TODO: make a documentation page for the GBZ format and link to it from here

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

.. warning::
  You must have `gbz-base <https://github.com/jltsiren/gbz-base/tree/main>`_ installed if you are working with GBZ files.

  .. code-block:: bash

    conda install -c conda-forge aryarm::gbz-base
  ..
    TODO: Once gbz-base is published on bioconda, we should recommend installing it from there instead.
  

Output
~~~~~~
  ..
    TODO: Describe this output in the formats section of our docs once the output format has stabilized.

The output is a tab-separated file with the following columns:

1. **numnodes**: The number of nodes in the region
2. **total_length**: The total length of all nodes in the region
3. **numwalks**: The number of walks in the region
4. The complexity metrics requested by ``--metrics``. See the :ref:`formulas <formulas>` section for more information.

If the ``--region`` option is specified, there will be one line in the output for every region. Each line will also be prefixed by the following columns:

1. **chrom**: The chromosome of the region
2. **start**: The start position of the region
3. **end**: The end position of the region


Examples
~~~~~~~~
By default, tab-separated output is written to standard out.
.. code-block:: bash

  panct complexity tests/data/basic.gfa

If your input graph is in the GBZ format, you may also use the :code:`--region` option to select a specific region of the graph in the coordinates of the reference genome. Eternally, this uses the gbz-base library to first subset the GBZ to a smaller GFA file.
.. code-block:: bash

  panct complexity --region chrTest:0-1 tests/data/basic.gbz

You may also specify a list of regions as a BED file, instead. In this case, it might also be helpful to write output to a file.
.. code-block:: bash

  panct complexity --out basic.tsv --region tests/data/basic.bed tests/data/basic.gbz


All files used in these examples are described :doc:`here </project_info/example_files>`.

Additional examples
~~~~~~~~~~~~~~~~~~~

Below are additional examples based on the HPRC .gbz format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  # Run on a single region
  panct complexity \
    --region chr11:119077050-119178859 \
    --metrics sequniq-normwalk,sequniq-normnode \
    hprc-v1.1-mc-grch38.gbz

  # Run on a file with a list of regions
  panct complexity \
    --region regions.bed --out test.tsv \
    --metrics sequniq-normwalk,sequniq-normnode \
    hprc-v1.1-mc-grch38.gbz

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: complexity
