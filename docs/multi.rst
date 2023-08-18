.. _multi:


Multi-Orient Extraction (`~slitlessutils.modules.extract.multi`)
================================================================

The multi-orient extraction was first developed by `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ and referred to as *LINEAR*.  This method was founded to address the contamination (or confusion) from overlapping spectral traces, where the detector records a weighted sum of the the spectra (see :numref:`confusion` below).  Therefore, this represents a fundamental degeneracy, which can only be broken with additional data.  The :doc:`single-orient extraction <single>` formalism uses the broadband data to inform the contamination model, whereas the *LINEAR* uses data at multiple orients provide a self-consistent, spectroscopic model for the entire scene.  

.. _confusion:
.. figure:: images/confusion.png
   :align: center
   :alt: Illustration of confusion

   An illustration of spectral contamination or confusion (taken from `Ryan, Casertano, & Pirzkal 2018 <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_).  Along the left side of each panel, they show a red and blue point source whose celestial positions are fixed (as shown the coordinate vane).  This shows that under certain orients (lower panel), the spectra will overlap leading to the spectral degeneracy.  However, if the telescope is reoriented, then the spectral traces separate, which provides the leverage to break this degeneracy.
