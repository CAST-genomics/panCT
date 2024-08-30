.. _manual-main:

panCT
=====

panCT is our pangenome complex trait toolkit, which implements command utilities for complex trait analyses using pangenomes

Commands
~~~~~~~~

* :doc:`panct complexity </commands/complexity>`: Compute complexity scores for a GFA file

Detailed information about each command can be found in the *Commands* section of our documentation. Examples there utilize a set of example files described :doc:`here </project_info/example_files>`.

Logging
~~~~~~~

All commands output log messages to standard error. The universal ``--verbosity`` flag controls the level of detail in our logging messages. By default, this is set to ``INFO``, which will yield errors, warnings, and info messages. To get more detailed messages, set it to ``DEBUG``. To get only error messages, set it to ``ERROR``. To get errors *and* warnings, set it to ``WARNING``. Refer to `the Python documentation on logging levels <https://docs.python.org/3/library/logging.html#levels>`_ for more information.

Contributing
~~~~~~~~~~~~

We gladly welcome any contributions to ``panct``!

Please read our :doc:`contribution guidelines </project_info/contributing>` and then submit a `Github issue <https://github.com/cast-genomics/panct/issues>`_.

Citing
~~~~~~

We do not yet have a manuscript for our toolkit.


.. toctree::
   :caption: Overview
   :name: overview
   :hidden:
   :maxdepth: 1

   project_info/installation
   project_info/example_files
   project_info/contributing

.. toctree::
   :caption: File Formats
   :name: formats
   :hidden:
   :maxdepth: 1

   formats/gfa.rst

.. toctree::
   :caption: Commands
   :name: commands
   :hidden:
   :maxdepth: 1

   commands/complexity.rst

.. toctree::
   :caption: API
   :name: api
   :hidden:
   :maxdepth: 1

   api/modules
