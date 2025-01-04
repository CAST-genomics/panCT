.. _commands-walks:


walks
=====

Extracts the walks of a pangenome graph.

The ``walks`` command outputs a file with walks for an entire graph.

Usage
~~~~~
.. code-block:: bash

  panct walks \
    --out PATH \
    --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
    GFAFILE

Examples
~~~~~~~~
To create a ``.walk.gz`` and ``.walk.gz.tbi`` file adjacent to a GFA, just specify the path to the GFA.

.. code-block:: bash

  panct walks tests/data/basic.gfa

You can also specify the output file path manually. If the output path does not have a ``.gz`` ending, an index will not be created automatically.

.. code-block:: bash

  panct walks --out basic.walk tests/data/basic.gfa

All files used in these examples are described :doc:`here </project_info/example_files>`.

Additional examples
~~~~~~~~~~~~~~~~~~~

Below are additional examples based on the HPRC .gfa format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  panct walks hprc-v1.1-mc-grch38.gfa

This may take some time. You can speed it up by allocating more memory and CPUs. The command should scale automatically with the resources that have been allocated to it.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: walks
