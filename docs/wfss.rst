.. _wfss:
.. sectionauthor:: Russell Ryan <rryan@stsci.edu>
.. codeauthor:: Russell Ryan <rryan@stsci.edu>

Wide-Field Slitless Spectroscopy
================================

:term:`wide-field slitless spectroscopy` (WFSS) refers to either the act of obtaining a complete, unbiased spectroscopic view of the sky by means of transmissive and dispersive optical element (usually a :term:`grism` or :term:`prism`) **and** the data produced by such an observation.  This is to distinguish other ways of using these optical elements, such as spatial/drift scanning, spectro-polarimetry, or detailed transient spectroscopy, which are not implemented in ``slitlessutils``.  It is unclear if/when such modes will be implemented.

Below the two main *container* datatypes that one generally interacts with are the :class:`~slitlessutils.wfss.WFSS()` and :class:`~slitlessutils.wfss.WFSSCollection()` to hold a single WFSS file/image and a set of such objects, respectively.  Since some of the extraction techniques implemented in ``slitlessutils`` can be computationally expensive and/or use a significant amount of memory, care has been taken to limit the memory footprint of objects (for example, the image-objects do not instantiate the data until requested).


WFSS data (`~slitlessutils.core.wfss.data.WFSS()`)
--------------------------------------------------
A light weight object that describes a single WFSS image.  This object will emulate a ``dict``-like object, where the keyword/value pairs are the detector name and a :class:`~slitlessutils.wfss.WFSSDetector()` object, which is unlikely to be directly instantiated by a user.  The :class:`~slitlessutils.wfss.WFSS()` objects are generally instantiated by one of two classmethods:

* :func:`~slitlessutils.wfss.WFSS.simulated()`: load the WFSS for simulation.

* :func:`~slitlessutils.wfss.WFSS.observed()`: load the WFSS as an observed file.



WFSS Collection (`~slitlessutils.core.wfss.data.WFSSCollection()`)
------------------------------------------------------------------

While one can instantiate a single WFSS file using the above, it is generally more common to load many of the files as a collection, as a :class:`~slitlessutils.wfss.WFSSCollection()` is the primary input to many of the additional functions/modules.  The :class:`~slitlessutils.wfss.WFSSCollection()` will act like a ``dict``, where the keyword/value pairs are the dataset name and *file-loading key*.  These keys are ``dataclass`` that are for loading observed and simulated data:




Observed Data (`~slitlessutils.wfss.data.wfsscollection.ObservedData()`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load observed fits files, there are three classmethods

* :func:`~slitlessutils.wfss.WFSSCollection.from_glob()`:  Load an observed ``WFSSCollection`` from a glob token.

* :func:`~slitlessutils.wfss.WFSSCollection.from_list()`:  Load an observed ``WFSSCollection`` from a list of files.

* :func:`~slitlessutils.wfss.WFSSCollection.from_file()`:  Load an observed ``WFSSCollection`` from an ascii file that contains the fullpath to a file (one per line).


Example
~~~~~~~
.. code:: python

	import slitlessutils as su

	# load as a glob string
	data1 = su.wfss.WFSSCollection.from_glob('*flt.fits')

	# load from a list
	filenames = ['file1_flt.fits', 'file2_flt.fits', 'file3_flt.fits']
	data2 = su.wfss.WFSSCollection.from_list(filenames)

	# load from an ascii file that contains filenames
	data3 = su.wfss.WFSSCollection.from_file(filelist)


Simulated data (`~slitlessutils.wfss.data.wfsscollection.SimulatedData()`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load simulated fits files, there are two classmethods

* :func:`~slitlessutils.wfss.WFSSCollection.from_dataframe()`:  Load a simulated ``WFSSCollection`` from a `pandas.DataFrame() <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html>`_.  The table shall have columns and datatypes given in :numref:`simulkeys`.

* :func:`~slitlessutils.wfss.WFSSCollection.from_wcsfile()`:  Load a simulated ``WFSSCollection`` from a comma-separated value (CSV) file whose columns and datatypes are given in :numref:`simulkeys` and an example CSV file is shown below.

.. _simulkeys:
.. list-table:: Simulated Data Keywords
   :widths: 15 10 75
   :header-rows: 1
   :stub-columns: 0
   :width: 600

   * - Keyword
     - Datatype
     - Notes
   * - dataset
     - ``str``
     - the file basename (ie. the IPPPSOOT for HST files)
   * - ra
     - ``float``
     - the right ascension of the instrument's reference point (in degrees)
   * - dec
     - ``float``
     - the declination of the instrument's reference point (in degrees)
   * - orientat
     - ``float``
     - the position angle of the instrument at its reference point (in degree)
   * - telescope
     - ``str``
     - the name of the telescope (e.g. ``HST``)
   * - instrument
     - ``str``
     - the tokenized name of the instrument (e.g. ``WFC3IR`` or ``ACSWFC``)
   * - disperser
     - ``str``
     - the name of the dispersive optic (e.g. ``G800L`` or ``G102``)
   * - blocking
     - ``str``
     - | the name of the blocking filter.  For HST instruments, this is to be left blank,
       | or set to ``None``, but is reserved for future development to support JWST
       | instruments or ACS spectro-polarimetry.


Example WCS File (in csv format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. include:: include/wcs.csv
	:literal:
