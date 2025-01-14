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

Below is an additional example based on the HPRC .gfa format graph (not included in this repo but available `here <https://github.com/human-pangenomics/hpp_pangenome_resources>`_). 

.. code-block:: bash

  panct walks hprc-v1.1-mc-grch38.gfa

This may take some time. You can speed it up by allocating more memory and CPUs. The command should scale automatically with the resources that have been allocated to it. On our system with 16 CPUs and 32 GB of RAM, we had to reserve 8 hours.

.. warning::
  The ``walks`` command may store thousands of small files (totaling nearly 100 GB) in your system's temporary directory. On some systems, this directory may not have the capacity for walks from the *entire* HPRC pangenome. So we recommend setting the ``$TMPDIR`` environment variable to a larger directory when calling the ``walks`` command.

  .. code-block:: bash

    TMPDIR=/scratch/$USER/job_$SLURM_JOBID panct walks hprc-v1.1-mc-grch38.gfa

  Ideally, you would choose a directory that is optimized for I/O-intensive work. In this example, we are using node-local scratch space on UCSD's TSCC or Expanse HPCs, for instance.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: walks
