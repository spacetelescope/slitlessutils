.. _tabulation:


Tabulation in ``slitlessutils``
===============================

The most computationally expensive aspect of extracting or simulating WFSS observations comes from the forward-modeling every (relevant) pixel in the :term:`direct imaging` 


Creating a :term:`pixel-dispersion table` (PDT):

# For each WFSS file:
	# For each detector in the WFSS file:
		# For each source in the source collection:
			# For each pixel in the source:
				#. apply WCS transformation between direct image and WFSS image
				#. For each spectral order:
					#. For each tabulation wavelength:
						#. convert wavelength into parameter :math:`t` using the inverse of the dispersion relation
						#. evaluate the trace at the paramter :math:`t`
						#. compute fractional pixel area (see :numref:`animatedpixel` below)
						#. record the fractional pixel area for each WFSS pixel multiplied by the bandwidth from the tabulation wavelengths.


.. _animatedpixel:
.. figure:: images/pixel_animate.gif
   :align: center
   :alt: fractional pixel animation

   Dispersed pixel and fractional area calculations.  ``Slitlessutils`` uses `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_ to compute fractional pixel area on a dispersed image pixel grid (shown by colored polygons).  The known area of the input polygon (shown in blue) is :math:`0.64 \mathrm{pix}^2`.  




