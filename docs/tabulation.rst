.. _tabulation:


Tabulation in ``slitlessutils``
===============================

The most computationally expensive aspect of extracting or simulating WFSS observations comes from the forward-modeling every (relevant) pixel in the :term:`direct imaging`.  Therefore, ``slitlessutils`` only performs these calculations when requested and stores these intermediate results as a :term:`pixel-dispersion table` (PDT).  

* For each WFSS file:
	- For each detector in the WFSS file:
		> For each source in the source collection:
			+ For each pixel in the source:
				#. apply WCS transformation between direct image and WFSS image
				#. For each spectral order:
					#. For each tabulation wavelength:
						#. convert wavelength into parameter :math:`t` using the inverse of the dispersion relation
						#. evaluate the trace at the parameter :math:`t`
						#. compute fractional pixel area (see :numref:`animatedpixel` below)
						#. record the fractional pixel area for each WFSS pixel multiplied by the bandwidth from the tabulation wavelengths.

.. _animatedpixel:
.. figure:: images/pixel_animate.gif
   :align: center
   :alt: fractional pixel animation

   Dispersed pixel and fractional area calculations.  ``Slitlessutils`` uses `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_ to compute fractional pixel area on a dispersed image pixel grid (shown by colored polygons).  The known area of the input polygon (shown in blue) is :math:`0.64 \mathrm{pix}^2`.  



Given the hierarchical nature outlined in the above algorithm, the PDTs are stored as `hierarchical data-format 5 (HDF5) <https://www.hdfgroup.org/solutions/hdf5/>`_ and then can be viewed or manually edited with standard tools (e.g. `HDFView <https://www.hdfgroup.org/downloads/hdfview/>`_).  Now the process of extraction or simulation will begin by aggregating the PDTs from the appropriate :term:`direct imaging` pixels and spectral order, then :term:`decimating`



