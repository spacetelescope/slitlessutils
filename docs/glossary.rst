.. _glossary:

Glossary
========

.. glossary::
	cosmic ray
		A high energy particle that imparts a significant amount of charge.  Cosmic rays often appear as very sharp, discontinuous features that do not reappear in successive images.

	data-quality array
		A bitmask containing the flags that indicate different failure modes of a pixel.  In general, the *good* pixels have a data-quality value of 0.  Also called DQA.

	direct imaging
		Data collected in a standard broadband imaging filter contemporaneously with the WFSS data.  These data are often used to improve the astrometric information in the WFSS data, but may have utility in other ways (such as determining the cross-dispersion profile weights, specifying the extraction apertures, or estimating contamination via broadband colors).  See also :term:`post-imaging` or :term:`pre-imaging`.

	grism
	    A transmissive and dispersive spectral element that often has a (nearly) constant rate of dispersion.  A grism differs from a :term:`prism` by having an additional diffractive surface on one side, which results in the constant dispersion, little spatial offset between the :term:`spectral trace` and :term:`undispersed position`, and multiple spectral orders.  See also :term:`prism`.

	master-sky image
		A model of the sky background for a WFSS image.  In principle, one should have a separate master-sky image for each spectral component present in the sky background.  However, only the HST/WFC3-IR detector has multiple-components measured (see `WFC3_back_sub <https://github.com/NorPirzkal/WFC3_Back_Sub>`_).

	pick-off mirror
	    An optical element that redirects the light to the instrument in question. 

	pixel-dispersion table
		A look-up table that describes the weights that a :term:`direct image` pixel projects onto the pixels of a WFSS image/detector as a function of wavelength for each spectral order.  Due to the hierarchical nature of this transformation, these data are stored as `hierarchical data-format 5 (HDF) <https://www.hdfgroup.org/solutions/hdf5/>`_.  This intermediate data product is also referred to as a PDT.

	post-imaging
	 	The contemporaneously :term:`direct imaging` taken *after* the WFSS observation.  See also :term:`pre-imaging`.

	pre-imaging
	 	The contemporaneously :term:`direct imaging` taken *before* the WFSS observation.  See also :term:`post-imaging`.

	prism
		A transmissive and dispersive spectral element with a highly non-uniform rate of dispersion.  See :term:`grism` for the similarities/differences between the two.

	spectral dispersion
		The parametric curve governing the wavelength along the :term:`spectral trace`.  Sometimes called the *wavelength solution*.  

	spectral trace
		The observed position of the two-dimensional spectra on the detector.  

	undispersed position
	    The position a source would have in the absence of the spectral grating: :math:`(x_0,y_0)`.  Importantly, this is **not** equivalent to the zeroth-order spectral trace.

	wide-field slitless spectroscopy
		The broad term for the use of a transmissive and dispersive optic to provide a complete, unbiased spectroscopic view of a scene.  This term may also refer to the data product of a single exposure/file taken through one of these optics.  May also be called WFSS for short.

	wedge offsets
		Positional offsets due to the variations in the thickness of the optical elements (see `Sabbi 2012 <https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2012/WFC3-2012-01.pdf>`_).