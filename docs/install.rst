.. _install:

Installing ``slitlessutils``
============================


Preparing your Local Environment
--------------------------------
We recommend using Anaconda to manage your ``slitlessutils`` environment.

You may want to consider installing ``slitlessutils`` in a new virtual or
conda environment to avoid version conflicts with other packages you may
have installed.

Start by creating an empty conda environment:

.. code-block:: bash

    conda create --name slitlessutils-env "python>=3.11"
    conda activate slitlessutils-env

You may also need to manually install precompiled ``llvmlite`` binaries before installing
``slitlessutils`` to avoid an error when building ``numba`` (used by ``slitlessutils`` to
improve performance) during the installation. This can be done by running ``conda install --channel=numba llvmlite``
in your conda environment before installing ``slitessutils``.


Installing ``slitlessutils``
----------------------------

There are several ways to install ``slitlessutils``:

* **Install from PyPI:** to install the latest stable version:

  .. code-block:: bash

    pip install slitlessutils


* **Install from GitHub:** to install the latest development version:

  .. code-block:: bash

    pip install git+https://github.com/spacetelescope/slitlessutils.git


* **Clone from GitHub:** the development version can also be cloned directly
  from GitHub:

  .. code-block:: bash

    git clone https://github.com/spacetelescope/slitlessutils.git
    cd slitlessutils
    pip install .

  This approach is only somewhat more difficult than installing from PyPI,
  but is easier to make local modifications and submit pull requests.
