.. _glossary:

Glossary
========

Since many terms may be used colloquially and/or have different definitions in other contexts, this Glossary provides a concrete definition for ambiguous terms.


.. glossary::
	
	attitude
		This may have many definitions, but in this context it refers to the direction and orientation that the telescope is pointed.  It is largely set by the CRVAL-keywords and the position angle, which is encoded in either the CD- or PC-matrices.  Here, this is considered synonymous with *pointing*.

	compound source
		A source that is decomposed into more than one :term:`dispersed region`.

	contamination model
		A model that describes the amount of light from an unrelated source that adversely affects the flux of the source in question.  These models are built on existing observations, usually broadband photometry, but can be spectroscopic data as well.  The concept of a *contamination model* only pertains to the :doc:`Single-Exposure Extraction <single>`, as the :doc:`Multi-Exposure Extraction <multi>` uses data at multiple position angles to mitigate contamination (see `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_).

	cosmic ray
		A high energy particle that imparts a significant amount of charge.  Cosmic rays often appear as very sharp, discontinuous features that do not reappear in successive images.

	damping target
		The vector of spectra that the sparse-least squares solutions will tend to minimize high-frequency noise.  See scipy implementation of the `LSQR <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html>`_ algorithm.

	data-quality array
		A bitmask containing the flags that indicate different failure modes of a pixel.  In general, the *good* pixels have a data-quality value of 0.  Also called DQA.

	direct imaging
		Data collected in a standard broadband imaging filter contemporaneously with the WFSS data.  These data are often used to improve the astrometric information in the WFSS data, but may have utility in other ways (such as determining the cross-dispersion profile weights, specifying the extraction apertures, or estimating contamination via broadband colors).  See also :term:`post-imaging` or :term:`pre-imaging`.

	dispersed region
		A subset of a source that has a single spectrum.  If a source has a single dispersed region, then it is said to be a :term:`simple source`.  Alternatively, a source that is decomposed into a many dispersed then it is a :term:`compound source`.

	grism
	    A transmissive and dispersive spectral element that often has a (nearly) constant rate of dispersion.  A grism differs from a :term:`prism` by having an additional diffractive surface on one side, which results in the constant dispersion, little spatial offset between the :term:`spectral trace` and :term:`undispersed position`, and multiple spectral orders.  See also :term:`prism`.

	master-sky image
		A model of the sky background for a WFSS image.  In principle, one should have a separate master-sky image for each spectral component present in the sky background.  However, only the HST/WFC3-IR detector has multiple-components measured (see `WFC3_back_sub <https://github.com/NorPirzkal/WFC3_Back_Sub>`_).

	pick-off mirror
	    An optical element that redirects the light to the instrument in question. 

	pixel-area map
		The relative area of each pixel with respect to the area of the reference pixel, which is given by the absolute value of the determinant of the Jacobian matrix.  This arises due to geometric distortion, and in the case of a SIP distortion model is a polynomial in the pixel coordinates.  Also called PAM.

	pixel-dispersion table
		A look-up table that describes the weights that the :term:`direct imaging` pixel projects onto the pixels of a WFSS image/detector as a function of wavelength for each spectral order.  Due to the hierarchical nature of this transformation, these data are stored as `hierarchical data-format 5 (HDF) <https://www.hdfgroup.org/solutions/hdf5/>`_.  This intermediate data product is also referred to as a PDT.

	post-imaging
	 	The contemporaneously :term:`direct imaging` taken *after* the WFSS observation.  See also :term:`pre-imaging`.

	pre-imaging
	 	The contemporaneously :term:`direct imaging` taken *before* the WFSS observation.  See also :term:`post-imaging`.

	prism
		A transmissive and dispersive spectral element with a highly non-uniform rate of dispersion.  See :term:`grism` for the similarities/differences between the two.

	regularization parameter
		A tunable parameter that governs the relative importance of fitting the data and damping high-frequency noise.  In the literature this will often be denoted by :math:`\lambda`, but of obvious confusion with wavelength is given the symbol :math:`\ell` in the ``slitlessutils`` discussion.  This may also be referred to as the *damping parameter*.

	segmentation map
		An image that describes which :term:`direct imaging` pixels belong each object, which effectively sets the extraction/simulation apertures and is used to initialize the :term:`dispersed region` for the sources.

	sensitivity curve
		The conversion between instrumental units (usually :math:`e^-/s`) to physical units (usually :math:`erg/s/cm^2/Å`), which is necessarily a function of wavelength.  

	signal-to-noise
		An empirical estimate of the quality of the data by comparing the measurement (the signal) to its corresponding uncertainty (the noise).  This may also be referred to as S/N or quoted as a *number of sigma* (:math:`n_{sig}`).

	simple source
		A source that has a single :term:`dispersed region`.

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

	world-coordinate system
		The complete description of the instrument layout on the sky, including the :term:`attitude` of the telescope, the relative position of the detectors, and their individual distortion models.  Also called WCS, and see also: `world-coordinate system <https://docs.astropy.org/en/stable/wcs/>`_.

	zeropoint
		The magnitude corresponding to 1 unit of flux (typically given as :math:`e^-/s`).
