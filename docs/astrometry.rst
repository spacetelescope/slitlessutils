.. _astrometry:

Astrometry
==========

To make full use of the WFSS data, it is necessary to take contemporaneous data in a standard imaging mode at the same position on the sky (colloquially called :term:`pre-imaging`, :term:`post-imaging`, or :term:`direct imaging`).  These data are used to deduce the exact position of the WFSS data by comparing to either existing imaging or a standard astrometric catalog (ie. Gaia; see `Mack et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022acs..rept....3M/abstract>`_).  Therefore, it is critical to ensure that any tweaks to the :term:`world-coordinate system` (WCS) of the direct imaging are propagated to the corresponding WFSS imaging.  Therefore ``slitlessutils`` provides functionality to sync the WCS between the direct and WFSS data.

The standard tools for updating the WCS keywords (ie. `drizzlepac and tweakreg <https://drizzlepac.readthedocs.io/en/latest/>`_) record previous WCS solutions in the `WCSNAME*` keywords, where `WCSNAME` refers to the *active* solution and earlier solutions are encoded with alphabetic characters (e.g. `WCSNAMEA`).  This last token is referred to as the *key*, and so one must specify the starting (``key0``) and the ending (``key``) key.


For much of the HST data, the direct imaging has been matched to some existing astrometric catalog (often Gaia dr3), but the WFSS data is left at a less-refined state.

.. important::
	**Astrometric corrections for the WFSS data from HST are almost always needed.**

The corrections can be either downgrading the astrometry in the direct image to match that in the WFSS, or upgrading the WFSS to match that of the direct image.


.. note::
	It is the responsibility of the user to verify that the WCS in the WFSS and direct image(s) are on the same astrometric reference, which is usually (but not always) answered by the ``WCSNAME`` keyword.  However, if the WCS solutions are inconsistent, then any subsequent spectral extraction or modeling will be assuredly incorrect.


Downgrading WCS (`~slitlessutils.core.preprocess.astrometry.downgrade_wcs()`)
-----------------------------------------------------------------------------

Here we roll back the "active" WCS in a **direct image** to an earlier version, which is meant to match that of the WFSS data.  This is usually the simpler option, but it will hamper using other (likely deeper) direct imaging for selecting the spectral traces for extraction.  In this case, we simply take the ``CRVAL`` and ``CD`` keywords from a previous WCS solution and replace them in the active solution, but there is a direct API for this:

.. code:: python

	import slitlessutils as su

	newfile = su.core.preprocess.astrometry.downgrade_wcs('a_direct_image_file_flt.fits', key='A')

will replace the WCS solution from ``WCSNAMEA`` into the active solution (ie. ``WCSNAME``), and return the name of the file that contains this WCS.  One can also pass a list of filenames to :func:`downgrade_wcs()`, which would then return a list of filenames.

.. warning::
	We expect this method will change as the astrometry description changes.

Upgrading WCS (`~slitlessutils.core.preprocess.astrometry.upgrade_wcs()`)
-------------------------------------------------------------------------

Here we find the affine transformation between any two WCSs from the direct image, and apply it to the equivalent WCS in the grism image.  For this, ``slitlessutils`` loads the ``CD`` matrix from a `astropy WCS <https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS>`_ object as:

.. math::
	\mathrm{CD} = \left(\begin{array}{cc}
			            \mathrm{CD}1\_1 & \mathrm{CD}1\_2 \\
			   			\mathrm{CD}2\_1 & \mathrm{CD}2\_2 \end{array}\right)

If the initial and final ``CD``-matrices for the direct image are :math:`\mathrm{CD}_{0,d}` and :math:`\mathrm{CD}_{1,d}` (respectively), then the affine transformation matrix is given as:

.. math::
	A = \mathrm{CD}_{1,d} \mathrm{CD}^{-1}_{0,d}

This transformation matrix is then applied to the ``CD`` matrix from the WFSS image (:math:`\mathrm{CD}_0`):

.. math::
	\mathrm{CD}_1 = A \mathrm{CD}_0

Similarly, we must adjust the ``CRVAL`` keywords, which are loaded from the WCS object as:

.. math::
	\mathrm{CRVAL} = \left(\begin{array}{c}\mathrm{CRVAL}1 \\
					\mathrm{CRVAL}2\end{array}\right)

Again, if :math:`\mathrm{CRVAL}_{0,d}` and :math:`\mathrm{CRVAL}_{1,d}` refer to the ``CRVAL``-vectors for the initial and final WCS solution from the direct image, then the perturbation vector is:

.. math::
	\Delta \mathrm{CRVAL} = \mathrm{CRVAL}_{1,d} - \mathrm{CRVAL}_{0,d}

which can be applied to the ``CRVAL`` vector from the WFSS image with the same WCS key used to instantiate :math:`\mathrm{CRVAL}_0`:

.. math::
	\mathrm{CRVAL}_1 = \mathrm{CRVAL}_0 + \Delta \mathrm{CRVAL}

This affine tweaking is implemented in the convenience function :func:`upgrade_wcs()`:

.. code:: python

	import slitlessutils as su

	newfile = su.core.preprocess.astrometry.upgrade_wcs('direct_image_reference_flt.fits', 'wfss_image_flt.fits')

This will compute the affine transformation between ``WCSNAME`` and ``WCSNAMEA`` from ``direct_image_reference_flt.fits`` and apply it to the ``WSCNAMEA`` astrometry in ``wfss_image_flt.fits``.  It writes the updated image to a new file, whose name is returned as ``newfile``.  The second argument can also be a list or tuple, and each of those images will be similarly tweaked, and in which case ``newfile`` will be a similar list.

.. note::
	In both above cases, a new file will be written if ``inplace==False`` and ``newfile==None``, and rules for generating this name are given in ``utils.py``.
