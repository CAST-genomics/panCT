.. _commands-walks:


walks
=====

Extracts the walks of a pangenome graph.

The ``walks`` command outputs a file with walks for an entire graph.

.. warning::
  The ``sort`` command must be installed. We recommend having at least GNU sort 9.12 (``conda install 'conda-forge::coreutils>=9.12'``). Other versions of sort will be much slower and require more temporary disk space.

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

  panct walks hprc-v2.0-mc-chm13.gfa

This may take some time. You can speed it up by allocating more memory. (Allocating more CPUs won't help as much as more RAM.) On our system, we allocated 4 CPUs, 64 GB of RAM, and 20 hours of compute.

.. warning::
  The ``walks`` command may store thousands of small files (totaling up to ~100 GB) in your system's temporary directory. On some systems, this directory may not have the capacity for walks from the *entire* HPRC pangenome. So we recommend setting the ``$TMPDIR`` environment variable to a larger directory when calling the ``walks`` command.

  .. code-block:: bash

    TMPDIR=/scratch/$USER/job_$SLURM_JOBID panct walks hprc-v2.0-mc-chm13.gfa

  Ideally, you would choose a directory that is optimized for I/O-intensive work. In this example, we are using node-local scratch space on UCSD's TSCC (with `--nodelist tscc-11-[5,14,16,18,19,27]`) or Expanse HPCs, for instance.

Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: walks
