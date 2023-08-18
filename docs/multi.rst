.. _multi:


Multi-Orient Extraction (`~slitlessutils.modules.extract.multi`)
================================================================

The multi-orient extraction was first developed by `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ and referred to as *LINEAR*.  This method was founded to address the contamination (or confusion) from overlapping spectral traces, where the detector records a weighted sum of the the spectra (see :numref:`confusion` below).  Therefore, this represents a fundamental degeneracy, which can only be broken with additional data.  The :doc:`single-orient extraction <single>` formalism uses the broadband data to inform the contamination model, whereas the *LINEAR* uses data at multiple orients provide a self-consistent, spectroscopic model for the entire scene.  

