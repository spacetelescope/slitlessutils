.. _wfss:
.. sectionauthor:: Russell Ryan <rryan@stsci.edu>
.. codeauthor:: Russell Ryan <rryan@stsci.edu>

Wide-Field Slitless Spectroscopy
================================

:term:`wide-field slitless spectroscopy` (WFSS) refers to either the act of obtaining a complete, unbiased spectroscopic view of the sky by means of tramissive and dispersive optical element (usually a :term:`grism` or :term:`prism`) **and** the data produced by such an observation.  This is to distinguish other ways of using these optical elements, such as spatial/drift scanning, spectro-polarimetry, or detailed transient spectroscopy, which are not implemented in ``slitlessutils``.  It is unclear if/when such modes will be implemented.

Below the two main *container* datatypes that one generally interacts with are the :class:`~slitlessutils.wfss.WFSS()` and :class:`~slitlessutils.wfss.WFSSCollection()` to hold a single WFSS file/image and a set of such objects.  Since some of the extraction techniques implemented in ``slitlessutils`` can be computationally expensive and/or use a significant amount of memory.  Therefore, care has been taken to limit the memory footprint of objects, for example the image-objects do not instantiate the data until requested.


WFSS data (`~slitlessutils.core.wfss.data.WFSS()`)
--------------------------------------------------
A light weight object that describes a single WFSS image.  This object will emulate a ``dict``-like object, where the keyword/value pairs are the detector name and a :class:`~slitlessutils.wfss.WFSSDetector()` object, which is unlikely to be directly instantiated by a user.  The :class:`~slitlessutils.wfss.WFSS()` objects are generally instantiated by one of two classmethods:

* :func:`~slitlessutils.wfss.WFSS.simulated()`: load the WFSS for simulation.

* :func:`~slitlessutils.wfss.WFSS.observed()`: load the WFSS as an observed file. 



WFSS Collections (`~slitlessutils.core.wfss.data.WFSSCollection()`)
-------------------------------------------------------------------

While one can instantiate a single WFSS file using the above, it is generally more common to load many of the files as a collection, as a :class:`~slitlessutils.wfss.WFSSCollection()` is the primary input to many of the additional functions/modules.  The :class:`~slitlessutils.wfss.WFSSCollection()`





Observed data
^^^^^^^^^^^^^
* :func:`~slitlessutils.wfss.WFSSCollection.from_glob()`

* :func:`~slitlessutils.wfss.WFSSCollection.from_list()`

* :func:`~slitlessutils.wfss.WFSSCollection.from_file()` (list of filenames)



Simulated data
^^^^^^^^^^^^^^

* :func:`~slitlessutils.wfss.WFSSCollection.from_dataframe()`

* :func:`~slitlessutils.wfss.WFSSCollection.from_wcsfile()`


.. list-table:: Title
   :widths: 25 25 50
   :header-rows: 1

   * - Heading row 1, column 1
     - Heading row 1, column 2
     - Heading row 1, column 3
   * - Row 1, column 1
     -
     - Row 1, column 3
   * - Row 2, column 1
     - Row 2, column 2
     - Row 2, column 3



.. include:: include/wcs.csv
	:literal: 






