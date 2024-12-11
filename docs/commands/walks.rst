.. _commands-walks:


walks
=====

Extracts the walks of a pangenome graph.

The ``walks`` command outputs a file with walks for an entire graph.

Usage
~~~~~
.. code-block:: bash

  panct walks \
    --out OUTFILE \
    --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
    GFAFILE

Examples
~~~~~~~~
.. code-block:: bash

  panct walks --out test.tab tests/data/test.gfa

All files used in these examples are described :doc:`here </project_info/example_files>`.

Additional examples
~~~~~~~~~~~~~~~~~~~

Below are additional examples based on the HPRC .gbz format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  # Run on a file
  panct walks hprc-v1.1-mc-grch38.gfa

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: walks
