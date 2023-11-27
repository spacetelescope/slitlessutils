.. _tabulation:


Tabulation in ``slitlessutils`` (`~slitlessutils.modules.Tabulate()`)
=====================================================================


The most computationally expensive aspect of extracting or simulating WFSS observations comes from the forward-modeling every (relevant) pixel in the :term:`direct imaging`.  Therefore, ``slitlessutils`` only performs these calculations when requested and stores these intermediate results as a :term:`pixel-dispersion table` (PDT).  These PDTs will be written to a subdirectory :file:`tables/` by default, but can be changed by setting the ``path`` keyword argument.  


Tabulation Algorithm
--------------------
To tabulate all the pixel transformations for a WFSS image and create a PDT, the algorithm iterates over all combinations of detector, source, direct-image pixel, spectral order, and wavelength according to:

* **INPUT:** a WFSS file:
	
	- For each detector in the WFSS file:
	
		> For each source in the source collection:
	
			+ For each :term:`direct imaging` pixel :math:`(x_d,y_d)` in the source:
	
				#. convert the direct imaging pixel to the :term:`undispersed position` :math:`(x_0,y_0)` using the WCS information for both images
				
				#. For each spectral order:
	
					* For each tabulation wavelength (see Note below):

						a. convert wavelength into parameter :math:`t` using the inverse of the dispersion relation
						b. evaluate the trace at the parameter :math:`t`
						c. compute fractional pixel area (see :numref:`animatedpixel` below)
						d. record an entry in the PDT as :math:`(x, y, l, a)`, where :math:`a` is the fractional pixel area that the :term:`direct imaging` pixel :math:`(x_d,y_d)` projects onto the WFSS image pixel :math:`(x,y)` at the wavelength index (see Note below)

* **OUTPUT:** a PDT written to disk.

.. note:: **The tabulation wavelengths:**
	The tabulation wavelengths are assumed to be linearly assigned *for all spectroscopic modes*, so that the wavelength is given:

	.. math::
		\lambda(l) = \lambda_0 + \left(\lambda_1-\lambda_0\right)\left(\frac{l}{N-1}\right)

	where :math:`N = \mathrm{ceil}\left(\frac{\lambda_1-\lambda_0}{\delta\lambda}\right)+1`, :math:`l\in(0,1,2,3,\ldots, N-1)` is the wavelength index, and :math:`\delta\lambda` is the sampling bandwidth.  The parameters :math:`(\lambda_0, \lambda_1, \delta\lambda)` are set the ``yaml`` files in the calibration directory: :file:`{$HOME}/.slitlessutils/<VERSION>/instruments/`.


.. _animatedpixel:
.. figure:: images/pixel_animate.gif
   :align: center
   :alt: fractional pixel animation

   Dispersed pixel and fractional area calculations.  ``Slitlessutils`` uses `pypolyclip <https://github.com/spacetelescope/pypolyclip>`_ to compute fractional pixel area on a dispersed image pixel grid (shown by colored polygons).  The known area of the input polygon (shown in blue outline) is :math:`0.64~\mathrm{pix}^2`.


Given the hierarchical nature outlined in the above algorithm, the PDTs are stored as `hierarchical data-format 5 (HDF5) <https://www.hdfgroup.org/solutions/hdf5/>`_ and the can be viewed or manually edited with standard tools (e.g. `HDFView <https://www.hdfgroup.org/downloads/hdfview/>`_).

Quick Primer on HDF5
^^^^^^^^^^^^^^^^^^^^

The HDF5 format is a "high-performance data management and storage suite" (`The HDF Group <https://www.hdfgroup.org/solutions/hdf5/>`_) that emulates the a file directory structure, where directories are referred to as "groups" and files are "datasets".  Each of these structures (groups or datasets) may contain some "attributes" that are effectively dictionary like keyword/value pairs that generally contain metadata.  Users who wish to inspect the precise nature/layout of any HDF5 file may find the `HDFView <https://www.hdfgroup.org/downloads/hdfview/>`_ graphical-user interface useful, which can be used to view or manually edit the HDF5 file.  The HDFGroup offers this editor for free (after registration), but it is **not required** for using ``slitlessutils``.  


Example
^^^^^^^

This example loads :doc:`sources <sources>` and :doc:`WFSS data <wfss>`, instantiates the tabulation module, and returns the names of the PDT files.

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


Use Cases
---------

The overall philosophy of ``slitlessutils`` is to compute these tables *once* at the outset, and use them for all downstream analyses, as they only contain the geometry of the astrophysical scene and the instrument/detector layout.  The primary use within ``slitlessutils`` begins with aggregating the PDTs from the appropriate :term:`direct imaging` pixels and spectral order, and summing over unique triplets :math:`(x,y,l)`.  These indices are combined into a single, unique index by `np.ravel_multi_index <https://numpy.org/doc/stable/reference/generated/numpy.ravel_multi_index.html>`_ following:

.. math::
	i = x + n_x\,y + n_x\,n_y\,l

where :math:`(n_x,n_y)` represents the dimensionality of the WFSS image.  This computation and summation is carried out by :func:`~slitlessutils.utilities.indices.decimate()`.

.. note::
	The PDT files only contain information on the scene geometry, and the detector effects and astrophysical signals are included in later stages.  Therefore these files **only depend** on the :term:`world-coordinate system` and its calibration.
