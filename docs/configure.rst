.. _configure:

Configuring ``Slitlessutils``
=============================

Many of the functions and classes within ``slitlessutils`` expose their
individual default settings, however in some cases it is important to
define *global* variables that govern calibration files or large-scale aspects
of the package.  Consequently, ``slitlessutils`` establishes these settings
in a singleton class :class:``~slitlessutils.Config`` that generally acts like 
a standard Python dictionary.  


``Slitlessutils`` Calibration Files
-----------------------------------
The reference files used to perform the spectral extraction and modeling with
``slitlessutils`` must be installed in a dot-directory in the user's home:
:file:`{$HOME}/.slitlessutils`.  When the :class:`~slitlessutils.Config` object is instantiated, 
it will check if this directory exists and valid reference files are populated.  If
no such directory is found, then it is created; if no such reference files
are found, then the most recent files will be automatically retrieved from
a public box directory.  However, users can programmatically retrieve and
use older versions of files:

.. code-block:: python

   # load the config module
   from slitlessutils import config
   cfg = config.Config()

   # download reference files
   reffile = cfg.retrieve_reffiles(refversion='1.0.1', update=True)

The :code:`update=True` flag will use these versions for the remainder of this
session.  One can list available reference libraries or swap between them:


.. code-block:: python
   
   # print all config files available
   cfg.help_refmanifest()

   # swap the reference file versions
   cfg.set_reffiles(refversion='1.0.1')


.. warning::
   Each time that ``slitlessutils`` is imported, it will use the most advanced
   version of the reference files that are cached in the user's home
   directory: :file:`{$HOME}/.slitlessutils`


Global Variables
----------------
To ensure consistency between several major subcomponents of ``slitlessutils``,
global variables are stored in the singleton configuration class.  These
variables are:

+----------------------+----------------+------------------------------------+
| variable name        | type           | Purpose                            |
|                      |                |                                    |
+======================+================+====================================+
| fluxscale            | float          | Numeric scale of the spectral      |
|                      |                | elements.  This is important to    |
|                      |                | avoid numeric underflow in some    |
|                      |                | cases.                             |
+----------------------+----------------+------------------------------------+
| fluxunits            | str            | The physical units of the spectral |
|                      |                | elements.                          |
+----------------------+----------------+------------------------------------+
| compression          | str            | The type of data compression used  |
|                      |                | in the HDF5 files as implemented   |
|                      |                | by `h5py <https://pypi.org/project/h5py/>`_.  |
+----------------------+----------------+------------------------------------+
| compression_opts     | int            | The compression level of used in   |
|                      |                | in the HDF5 files as implemented   |
|                      |                | by `h5py <https://pypi.org/project/h5py/>`_.  |
+----------------------+----------------+------------------------------------+


The default values are set in the file :file:`$HOME/.slitlessutils/<VERSION_NUMBER>/defaults.json`.
These versions can be accessed and/or adjusted programmatically as either
dictionary-like or attribute-like access and saved to a file for usage later:

.. code-block:: python

   # change value using dict-like access
   cfg['fluxscale'] = 1.

   # change value using attribute-like access
   cfg.fluxunits = 'erg/s/cm**2/micron'

   # save file to a local config
   cfg.write("myconf.json")

.. note::
   One can manually edit the defaults file, however new reference files packages will come
   with their own `defaults.json` file. Therefore, we recommend programmatically alter the
   settings to ensure consistency in results if the reference files are updated.
