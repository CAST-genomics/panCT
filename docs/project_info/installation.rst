.. _project_info-installation:

============
Installation
============

.. warning::
   panCT is still under construction and has not yet been published to PyPI or bioconda, so it cannot yet be installed through ``pip`` or ``conda``. We will leave these docs here for later.

Using pip
---------

You can install ``panct`` from PyPI using ``pip``.

.. note::
   We recommend using ``pip >= 20.3``.

   .. code-block:: bash

      pip install 'pip>=20.3'

.. code-block:: bash

   pip install panct

Using conda
-----------

We also support installing ``panct`` from bioconda using ``conda``.

.. code-block:: bash

   conda install -c conda-forge -c bioconda panct

Installing the latest, unreleased version
-----------------------------------------
Can't wait for us to tag and release our most recent updates? You can install ``panct`` directly from the ``main`` branch of our Github repository using ``pip``.

.. code-block:: bash

   pip install --upgrade --force-reinstall git+https://github.com/cast-genomics/panct.git
