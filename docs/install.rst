.. _installing:

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

    conda create --name slitlessutils-env "python>=3.10, <3.11"
    conda activate slitlessutils-env


Installing ``slitlessutils``
----------------------------

There are several ways to install ``slitlessutils``:


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
  but is easier to make local modifcations and submit pull requests.




Fetching Reference Data
-----------------------

There are a host of reference and configuration files needed by
``slitlessutils``, which have been staged at `slitlessutils data <BOX>`_.
You should download the `.tar.gz` file for instrument of interest,
copy it to the home directory, and unpack with:

.. code-block:: bash

   tar -xcvzf <FILENAME>.tar.gz

This will create a directory in the home as :code:`.slitlessutils`



Adjust Primary Defaults
-----------------------

The default configuration is in the file ``{$slitlessutils_config}defaults.cfg```

.. code-block:: python

   # load the config module
   import slitlessutils as su
   conf = su.config.Confing()

   # adjust an existing parameter (two options)
   conf.fluxscale = 2e-17
   conf['fluxunit'] = 'erg/(s*cm**2*micron)'

   # save config to disk
   conf.write("myconf.cfg")

   
