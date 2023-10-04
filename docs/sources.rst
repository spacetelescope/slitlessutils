.. _sources:

Spectroscopic Sources
=====================



Source (`~slitlessutils.sources.Source()`)
-----------------------------------------------

A :term:`source` is a single region on the sky that will be considered for processing.  A source can only be instantiated from two-dimensional images that describe the direct image and segmentation map.  The segmentation map describes the direct-image pixels that belong to the source, and so in the most basic terms, a source is the set of direct-image pixels:

.. math::
	\mathbb{S} = \left\{(x_1, y_1), (x_2,y_2), (x_3, y_3), ..., (x_n, y_n)\right\}

The direct image establishes weights that will be used for the cross-dispersion profile and metadata for the source (e.g. barycenter, brightness [including local-sky subtraction], etc.).  While these images may be determined from any filter, the best results will be when it is from the same camera and its transmission is contained within the spectroscopic bandpass (e.g. as F140W is to G141 for WFC3/IR).  

In addition to the metadata mentioned above, a source also describes the spectrum that is either the assumed or measured, in the case of :doc:`simulation <simulation>` or :doc:`extraction <extraction>`, respectively.  However, in either case, **source may contain multiple spectra** by decomposing a source into multiple :term:`dispersed regions<dispersed region>` (more below).  A dispersed region is a unique subset of the source pixels (denoted as :math:`\mathbb{S}_i` for the :math:`i^\mathrm{th}` dispersed region) and has a unique spectrum associated with it.  Obviously, the union of all dispersed regions is the same as the source set (:math:`\mathbb{S}=\bigcup_i\mathbb{S}_i`).  If a source contains a single dispersed region it is said to be a :term:`simple source`, which is contrast to a :term:`compound source` that contains many distinct spectral regions. 

.. note::
	Given increased (self-)contamination associated with a compound source, it can only be extracted if multiple orients are available and the :func:`sliltessutils.core.modules.extract.multi.Multi()` is used.  Additionally, care should be taken to avoid overly subdividing a source such that there are more "unknowns" than "knowns".

Since a source may contain multiple dispersed regions, the :func:`slitlessutils.core.sources.Source()` inherits from ``list``, which iterates over the dispersed regions. 


separate wavelength settings
talk about grism/prism dispersers





Dispersed Region (`~slitlessutils.sources.DispersedRegion()`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As stated above, a :term:`dispersed region` is defined as the unique subset of direct-image pixels that have the same spectrum, which is stored as an ``slitlessutils.core.photometry.SED`` object.  Importantly, the ``SED`` object contains the data one expects for a spectrum (e.g. wavelength, flux, uncertainty, etc.), the ``DispersedRegion`` object must contain additional metadata, including the direct-image pixels, the :term:`source ID`, and the :term:`region ID`.  By construction, every source has a unique :term:`source ID`: ``segid``, but each of its constituent ``DispersedRegion``s will be assigned a :term:`region ID`: ``regid``.  However, the ``regid`` always counts from 0 for each source, therefore the tuple ``(segid, regid)`` uniquely specifies an :term:``sed ID``.

Source Collection (`~slitlessutils.sources.SourceCollection()`)
--------------------------------------------------------------------

This is the primary data structure that users will interact with, which is meant to mimic the structure of the


Inputs are segmentation maps



.. _segmapexample:
.. figure:: images/animate_segmap.gif
	:align: center
	:alt: Animation of direct image and segmentation map

	Illustration of the direct image and (classic) segmentation map.


EXAMPLE



These definitions establish a *hierarchy*, where a ``SourceCollection`` (likely) contains many ``Source``\s that (potentially) contain many ``DispersedRegion``\s that (typically) contain many spectral elements (ie. wavelengths, fluxes, and uncertainties).  This hierarchy is show schematically in :numref:`hierarchy`, with the any :term:`compound source` highlighted in gray.

.. _hierarchy:
.. figure:: images/sourcecollection.png
	:align: center
	:alt: Schematic of source/spectra hierarchy

	Schematic representation of the source/spectra hierarchy with the primary inputs (segmentation map and direct image) shown.  A ``SourceCollection`` (purple box) is the primary way to instantiate a ``Source`` (blue circles), which contain any number of ``DispersedRegion``s (orange hexagons) that each contain one ``SED`` (red cylinder).  A :term:`compound source` is highlighted in gray.  




Notes on Extraction Parameters
------------------------------

The default extraction parameters are specified in the :doc:`instrument YAML files <instrumentfiles>`.  However, they can be programmatically changed at any of the level of the above hierarchy, and will be propagated to all of its children levels.

EXAMPLE




