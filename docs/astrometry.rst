.. _astrometry::

Astrometry
==========

To make full use of the WFSS data, it is necessary to take contemporaneous data in a standard imaging mode at the same position on the sky (colloquially called :term:`pre-imaging`, :term:`post-imaging`, or :term:`direct imaging`).  These data are often used to deduce the exact position of the WFSS data by comparing to either existing imaging or a standard astrometric catalog (ie. Gaia).  Therefore, it is critical to ensure that any tweaks to the world-coordinate system (WCS) of the direct imaging is propagated to the corresponding WFSS imaging.  Therefore ``slitlessutils`` provides functionality to sync the WCS between the direct and WFSS data.

The standard tools for updating the WCS keywords (ie. `drizzlepac and tweakreg <https://drizzlepac.readthedocs.io/en/latest/>`_) record previous WCS solutions in the `WCSNAME*` keywords, where `WCSNAME` refers to the *active* solution and earlier solutions are encoded with alphabetic characters (e.g. `WCSNAMEA`).  This last token is referred to as the *key*, and so one must specify the starting (``key0``) and the ending (``key``) key.


For much of the HST data, the direct imaging has been matched to some existing astrometric catalog (often Gaia dr3), but the WFSS data is left at a less-refined state.  

.. important::
	Astrometric corrections for the WFSS data from HST are almost always needed.

The corrections can be either downgrading the astrometry in the direct image to match that in the WFSS, or upgrading the WFSS to match that of the direct image.


Downgrading WCS (`~slitlessutils.core.preprocess.astrometry.downgrade_wcs()`)
-----------------------------------------------------------------------------
Not Implemented.



Upgrading WCS (`~slitlessutils.core.preprocess.astrometry.AffineTweak()`)
-------------------------------------------------------------------------

Here we find the affine transformation between any two WCSs from the direct image, and apply it to the equivalent WCS in the grism image.  For this, ``slitlessutils`` loads the ``CD`` matrix from a `astropy WCS <https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS>`_ object as:

.. math::
	\mathrm{CD} = \left(\begin{array}{cc}
			            \mathrm{CD}1_1 & \mathrm{CD}1_2 \\
			   			\mathrm{CD}2_1 & \mathrm{CD}2_2 \end{array}\right)	

If the initial and final ``CD``-matrices for the direct image are :math:`D_0` and :math:`D_1` (respectively), then the affine transformation matrix is given as:

.. math::
	A = D_1 D^{-1}_0

This matrix is then applied to the ``CD`` matrix from the WFSS image with the same WCSKEY as used to instantiate :math:`D_0`:

.. math::
	C_1' = A C_0

Similarly, we must adjust the ``CRVAL`` keywords, which are loaded from the WCS object as:

.. math::
	\mathrm{CRVAL} = \left(\begin{array}{c}\mathrm{CRVAL}1 \\ 
					\mathrm{CRVAL}2\end{array}\right)

Again, if :math:`d_0` and :math:`d_1` refer to the ``CRVAL``-vectors for the initial and final WCS solution, then the perturbation is:

.. math::
	\Delta = d_1 - d_0

which can be applied to the ``CRVAL`` vector from the WFSS image with the same WCS key used to instantiate :math:`d_0`:

.. math::
	d_1' = d_0 + \Delta

This affine tweaking can be derived from a direct image and applied to a WFSS image using ``slitlessutils``:

.. code:: python
	
	import slitlessutils as su

	tweak = su.core.preprocess.astrometry.AffineTweak(direct_image_filename)

	updated_wfss_filename = tweak(wfss_filename, inplace=False, newfile='new_wfssfile.fits')


