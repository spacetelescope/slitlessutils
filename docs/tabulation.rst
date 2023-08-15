.. _tabulation:


Tabulation in ``slitlessutils`` (`~slitlessutils.modules.Tabulate()`)
=====================================================================


The most computationally expensive aspect of extracting or simulating WFSS observations comes from the forward-modeling every (relevant) pixel in the :term:`direct imaging`.  Therefore, ``slitlessutils`` only performs these calculations when requested and stores these intermediate results as a :term:`pixel-dispersion table` (PDT).  These PDTs will be written to a subdirectory :file:`tables/`, which can be 


Tabulation Algorithm
--------------------
To tabulate all the pixel transformations for a WFSS image and create a PDT, the algorithm iterates over all combinations of detector, source, direct-image pixel, spectral order, and wavelength:

:class:`~slitlessutils.modules.Tabulate()`

* INPUT: a WFSS file:
	- For each detector in the WFSS file:
		> For each source in the source collection:
			+ For each :term:`direct imaging` pixel :math:`(x_d,y_d)` in the source:
				#. convert the direct imaging pixel to the :term:`undispersed position` :math:`(x_0,y_0)` using the WCS information for both images
				#. For each spectral order:
					#. For each tabulation wavelength (see Note below):
						#. convert wavelength into parameter :math:`t` using the inverse of the dispersion relation
						#. evaluate the trace at the parameter :math:`t`
						#. compute fractional pixel area (see :numref:`animatedpixel` below)
						#. record an entry in the PDT as :math:`(x, y, l, a)`, where :math:`a` is the fractional pixel area that the :term:`direct imaging` pixel :math:`(x_d,y_d)` projects onto the WFSS image pixel :math:`(x,y)` at the wavelength index (see Note below)
* OUTPUT: a PDT written to disk.

.. note:: **The tabulation wavelengths**
	The tabulation wavelengths are assumed to be linearly assigned *for all spectroscopic modes*, so that the wavelength is given:

	.. math::
		\lambda(l) = \lambda_0 + \left(\lambda_1-\lambda_0\right)\left(\frac{l}{N-1}\right)

	where :math:`N = \mathrm{ceil}((\lambda_1-\lambda_0)/\delta\lambda)+1`, :math:`l\in(0,1,2,3,\ldots, N-1)` is the wavelength index, and :math:`\delta\lambda` is the sampling bandwidth.  The parameters :math:`(\lambda_0, \lambda_1, \delta\lambda)` are set the ``yaml`` files in the calibration directory: :file:`{$HOME}/.slitlessutils`.


.. _animatedpixel:
.. figure:: images/pixel_animate.gif
   :align: center
   :alt: fractional pixel animation

   Dispersed pixel and fractional area calculations.  ``Slitlessutils`` uses `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_ to compute fractional pixel area on a dispersed image pixel grid (shown by colored polygons).  The known area of the input polygon (shown in blue) is :math:`0.64 \mathrm{pix}^2`.  


Given the hierarchical nature outlined in the above algorithm, the PDTs are stored as `hierarchical data-format 5 (HDF5) <https://www.hdfgroup.org/solutions/hdf5/>`_ and then can be viewed or manually edited with standard tools (e.g. `HDFView <https://www.hdfgroup.org/downloads/hdfview/>`_).  Now the process of extraction or simulation will begin by aggregating the PDTs from the appropriate :term:`direct imaging` pixels and spectral order, and summing over unique triplets :math:`(x,y,l)`.  These indices are combined into a single, unique index by `np.ravel_multi_index <https://numpy.org/doc/stable/reference/generated/numpy.ravel_multi_index.html>`_ following: 

.. math:: 
	i = x + n_x\,y + n_x\,n_y\,l

where :math:`(n_x,n_y)` represents the dimensionality of the WFSS image.  This computation and summation is carried out by :func:`~slitlessutils.utilities.indices.decimate()`.  


Example
-------

This example loads :doc:`sources <source>` and :doc:`WFSS data <wfss>`, instantiates the tabulation module, and returns the names of the PDT files.  

.. code:: python

	import slitlessutils as su

	# instantiate source from a segmentation image
	sources = su.source.SourceCollection(segfile, imgfile)

	# instantiate the spectral images from all the files matching some filename
	data = su.wfss.WFSSCollection.from_glob('*_flt.fits')

	# instantiate the tabulation object
	tab = su.modules.Tabulate(ncpu=2)

	# call the tabulation method
	pdtnames = tab(data, sources)



