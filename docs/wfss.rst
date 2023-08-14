.. _wfss:
.. sectionauthor:: Russell Ryan <rryan@stsci.edu>
.. codeauthor:: Russell Ryan <rryan@stsci.edu>

Wide-Field Slitless Spectroscopy
================================

:term:`wide-field slitless spectroscopy` (WFSS) refers to either the act of obtaining a complete, unbiased spectroscopic view of the sky by means of tramissive and dispersive optical element (usually a :term:`grism` or :term:`prism`) **and** the data produced by such an observation.  This is to distinguish other ways of using these optical elements, such as spatial/drift scanning, spectro-polarimetry, or detailed transient spectroscopy, which are not implemented in ``slitlessutils``.  It is unclear if/when such modes will be implemented.


WFSS data (`~slitlessutils.core.wfss.data.WFSS()`)
--------------------------------------------------



WFSS Collections (`~slitlessutils.core.wfss.data.WFSSCollection()`)
-------------------------------------------------------------------


Observed data
^^^^^^^^^^^^^


Simulated data
^^^^^^^^^^^^^^

.. include:: <ascii/wcs.csv>

print('dataset,ra,dec,orientat,telescope,instrument,disperser,blocking', file=fp)








