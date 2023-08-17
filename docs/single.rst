.. _single:


Simple-Exposure Extraction (`~slitlessutils.modules.extract.single`)
====================================================================

.. _extsec:

Extraction
----------

* For each WFSS image in the :doc:`WFSS Collection <wfss>`:
	* For each source in the :doc:`Source Collection <source>`:
		#. Load the PDTs for each :term:`direct imaging` pixel in the source.
		#. Decimate the PDTs over the :term:`direct imaging` pixel
		#. Decimate over wavelength to get the cross-dispersion profile
		#. Compute average and range of wavelength in each WFSS pixel
		#. Divide each pixel in the WFSS image by the flat-field, sensitivity curve, pixel-area map, and ``fluxscale`` (see :doc:`<config>`)
		#. Record these values in a temporary data structure

.. _expcombo:

Exposure Combination
--------------------



Flux Contamination
------------------
The :term:`contamination model` processes through the exact same steps in :ref:`Extraction <extsec>` and :ref:`Exposure Combination <expcombo>`.


