.. _single:


Single-Orient Extraction (`~slitlessutils.modules.extract.Single()`)
======================================================================

The *single-orient extraction* refers to extraction of one-dimensional spectra from a collection of WFSS images using a methodology very similar to that of `hstaxe <https://hstaxe.readthedocs.io/en/latest/>`_, with a few notable exceptions.

* Each pixel in the :term:`direct imaging` is dispersed (with its unique trace and wavelength model) and the fractional pixel area is collected (see `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_).
* This represents a complete forward model, which can be used in optimization strategy (such as fitting SED models with MCMC or greedy algorithm).
* The :term:`contamination model` (easily) permits arbitrary complexity for the contaminating source spectra.


.. note::
	The Single-Orient Extraction is only applicable to a :term:`simple source`, as the (self-)contamination is too severe to be addressed with this algorithm.


.. _extsec:

Extraction Details
------------------

The single-orient extraction is essentially the optimal spectroscopy algorithm presented by `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_, with a few modifications related to details of WFSS data and/or the preprocessing that is expected.  The method is 


* For each WFSS image in the :doc:`WFSS Collection <wfss>`:
	
	* For each source in the :doc:`Source Collection <sources>`:
	
		#. load the PDTs for each :term:`direct imaging` pixel in the source
		#. :term:`decimate<decimation>` the PDTs over the unique combinations of WFSS pixel and wavelength
		#. For each unique :math:`x`-coordinate:
			
			- compute average wavelength (weighted by the forward-model profile) for the :math:`y`-pixels in this vertical slice
			- divide each :math:`y`-pixel in the WFSS image by their :doc:`flat-field <calib>`, :doc:`sensitivity curve <calib>`, :doc:`pixel-area map <simulation>`, ``fluxscale`` (see the :doc:`configuration object <configure>`), and the instantaneous dispersion for this average wavelength
			- compute the `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_ optimal parameters:

			.. math::

				\begin{eqnarray}
					f_{\lambda,i} &=& \frac{\sum_y P_{x,y}S_{x,y}/U_{x,y}^2}{\sum_y P_{x,y}P_{x,y}/U_{x,y}^2}\\
					u_{\lambda,i} &=& \sqrt{\frac{\sum_y P_{x,y}}{\sum_y P_{x,y}P_{x,y}/U_{x,y}^2}}\\
					c_{\lambda,i} &=& \frac{\sum_y P_{x,y}C_{x,y}/U_{x,y}^2}{\sum_y P_{x,y}P_{x,y}/U_{x,y}^2}
				\end{eqnarray}
				
			where :math:`f_{\lambda,i}`, :math:`u_{\lambda,i}`, and :math:`c_{\lambda,i}` are the optimal flux, uncertainty, and contamination, respectively for the :math:`i^\mathrm{th}` WFSS image.  Additionally, :math:`S_{x,y}`, :math:`U_{x,y}`, :math:`P_{x,y}`, and :math:`C_{x,y}` are the science, uncertainty, cross-dispersion profile, and contamination images (more on this below in :ref:`Contamination Model <contmodel>`), respectively.  ``Slitlessutils`` offers three choices for the cross-dispersion profile :math:`P_{x,y}`:
				* **uniform** This does no profile weighting and instead just sums the pixels within the aperture.  This is effectively the box-extraction in `hstaxe <https://hstaxe.readthedocs.io/en/latest/>`_
				* **forward** This uses the forward model to establish the cross dispersion weights.  
				* **data** This uses the science image, masked for the :term:`DQA <data-quality array>` as the weights.  This is effectively the `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_ algorithm.
			
			these can be selected by the keyword argument: ``profile``.  The default behavior is ``profile='data'``.
			
			- record these values in a temporary data structure used to combine the spectra from different WFSS images


This produces a single spectrum for each source for each WFSS image, and these spectra are combined in the next section.  The two key differences between this and the `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_ algorithm are (1) there is no iterative reassessment of either the profile (:math:`P_{x,y}`), the cosmic ray mask, or the pixel variances; and (2) the science image is not smoothed in the dispersion direction.


.. important::
	In the above algorithm, vertical slices in :math:`x` are taken as the HST WFSS modes very nearly disperse along the :math:`x`-axis (though WFC3/UVIS has significant curvature).  However, these summations in the above should be carried out over a fixed wavelength interval, but this is very similar to the method in `hstaxe <https://hstaxe.readthedocs.io/en/latest/>`_.  This assumption will be addressed in future releases.



.. _contmodel:

Contamination Model
-------------------

The :term:`contamination model` is initialized by converting the spectral traces for every object into a polygon from the `Shapely <https://shapely.readthedocs.io/en/stable/>`_ from the information in the PDTs.  If the contamination is requested, then ``slitlessutils`` will select all polygons that intersect with the spectral-trace polygon in question.  For those that intersect, then a simulated image is created (see :doc:`Simulation<simulation>` for details), however this is only done in a small postage stamp.  This simulation is the contamination image :math:`C_{x,y}` described in the above section.  Therefore, the quality of the contamination model is directly related to the quality of the available broadband photometry.  Lastly, this methodology is the same as the ``fluxcube`` settings in `hstaxe <https://hstaxe.readthedocs.io/en/latest/>`_.


.. important::
	The contamination will be computed if-and-only-if the ``mskorders`` keyword is set.  This can be either a single string for the orders to mask, the string ``'all'`` to mask all orders, or ``None`` to mask no orders.  The default behavior is ``mskorders='all'``.  

.. note::
	The class setting :code:`savecont=True` to the ``Single()`` module will save the two-dimensional contamination models as a multi-extension fits file to disk, where each extension will refer to a different :term:`segmentation ID <source ID>`.


.. _expcombo:

Exposure Combination
--------------------

The results from the :ref:`Extraction <extsec>` module are combined into a single one-dimensional spectrum for each source.

* For each source in the :doc:`Source Collection <sources>`:

	* bin the wavelengths according to the extraction wavelengths
	* initialize the weights as the inverse of the square of the uncertainties: :math:`w=1/u^2`.
	* Compute the number of non-zero weights for each wavelength :math:`n_{\lambda}`, and the weighted moments of the photometric data:

	.. math::

		f_{\lambda} &=& \frac{\sum_i f_{\lambda,i}\,w_{\lambda,i}}{\sum_i w_{\lambda,i}}\\
		u_{\lambda} &=& \frac{1}{\sqrt{\sum_i w_{\lambda,i}}}\\
		c_{\lambda} &=& \frac{\sum_i c_{\lambda,i}\,w_{\lambda,i}}{\sum_i w_{\lambda,i}}

	where :math:`f_{\lambda}`, :math:`u_{\lambda}`, and :math:`c_{\lambda}` are the averaged spectrum, uncertainty, and contamination model that ``slitlessutils`` reports for this source

	* Output the table of :math:`\lambda`, :math:`f_{\lambda}`, :math:`u_{\lambda}`, and :math:`c_{\lambda}`, and :math:`n_{\lambda}` into an output fits file, whose suffix will be ``x1d.fits``.


.. important::
	The extraction wavelengths are specified for each source, but can be globally set for each of them in the :doc:`Source Collections <sources>`.


Example
-------

See :file:`slitlessutils.examples.starfield` for a working example.
