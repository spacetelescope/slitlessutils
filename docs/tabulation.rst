.. _tabulation:


Tabulation in ``slitlessutils``
===============================


#. For each WFSS file:
	#. For each detector in the WFSS file:
		#. For each source in the source collection:
			#. for each pixel in the source:
				#. apply WCS transformation
				#. for each tabulation wavelength:
					#. Invert
					#. fractional pixel area



.. _animatedpixel:
.. figure:: images/pixel_animate.gif
   :align: center
   :alt: fractional pixel animation

	Dispersed pixel and fractional area calculations.  ``Slitlessutils`` uses `pypolyclip <>`_ to compute fractional pixel area on a dispersed image pixel grid (shown by colored polygons).  The known area of the input polygon (shown in blue) is 0.64 pix:math:`^2`.  




