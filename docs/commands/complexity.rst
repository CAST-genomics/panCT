.. _commands-complexity:


complexity
==========

Compute the complexity of a pangenome graph.

The ``complexity`` command outputs a single float to ``stdout``.

Usage
~~~~~
.. code-block:: bash

  panct complexity \
  --verbosity [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET] \
  GFA

Examples
~~~~~~~~
.. code-block:: bash

  panct complexity tests/data/basic.gfa

All files used in these examples are described :doc:`here </project_info/example_files>`.


Detailed Usage
~~~~~~~~~~~~~~

.. click:: panct.__main__:typer_click_object
   :prog: panct
   :show-nested:
   :commands: complexity
