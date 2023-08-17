.. _single:


Simple-Exposure Extraction (`~slitlessutils.modules.extract.Single()`)
======================================================================

The *simple-exposure extraction* refers to extraction of one-dimensional spectra from a collection of WFSS images using a methodology very similar to that of `hstaxe <https://hstaxe.readthedocs.io/en/latest/>`_, with a few notable exceptions.

* Each pixel in the :term:`direct imaging` is dispersed (with its unique trace and wavelength model) and the fractional pixel area is collected (see `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_).
* This represents a complete forward model, which can be used in optimization strategy (such as fitting SED models with MCMC or greedy algorithm).  
* The :term:`contamination model` (easily) permits arbitrary complexity for the contaminating source spectra.


.. note::
	The Simple-Exposure Extraction is only applicable to a :term:`simple source`, as the (self-)contamination is too severe to be addressed with this algorithm.


.. _extsec:

Extraction
----------

* For each WFSS image in the :doc:`WFSS Collection <wfss>`:
	* For each source in the :doc:`Source Collection <sources>`:
		#. Load the PDTs for each :term:`direct imaging` pixel in the source.
		#. Decimate the PDTs over the :term:`direct imaging` pixel
		#. Decimate over wavelength to get the cross-dispersion profile
		#. Compute average and range of wavelength in each WFSS pixel
		#. Divide each pixel in the WFSS image by the :doc:`flat-field <calib>`, :doc:`sensitivity curve <calib>`, :doc:`pixel-area map <simulation>`, and ``fluxscale`` (see the :doc:`configuration object <configure>`)
		#. Record these values in a temporary data structure

.. _expcombo:

Exposure Combination
--------------------

The results from the :ref:`Extraction <extsec>` module are combined into a single one-dimensional spectrum for each source.  





Flux Contamination
------------------
The :term:`contamination model` processes through the exact same steps in :ref:`Extraction <extsec>` and :ref:`Exposure Combination <expcombo>`. Additionally, the one-dimensional spectra output from the previous step are **NOT** contamination corrected, but rather the estimated model is provided. Therefore, users are free to heuristically adjust the contamination model *post facto*, ignore regions with egregious contamination, or any other post-extraction choice.  Finally, the contamination model will have the same units as the source spectra and uncertainties, which will be set by the :doc:`configuration object <configure>`.



