
WFSS Image Simulation
=====================


Introduction
------------






Methodology
-----------





Included Effects
^^^^^^^^^^^^^^^^


The uncertainty image is given by the

.. math::
   U = \frac{\sqrt{(I+S+D) t+R^2}}{t}

where :math:`I`, :math:`S` [#f1]_, and :math:`D` are the Poissonian noise
terms that represent the flux (in :math:`e^-`/s) from the simulated
science image, the specified sky level, and the dark rate,
respectively.  The read noise (in :math:`e^-`) is specified as
:math:`R`, and represents the lone Gaussian noise term.  The specified
exposure time (in s) is given by :math:`t`.  Therefore, the simulated
images will have an `ERR` extension will be populated with these values.

The `SCI` extension

.. math::
   p \sim \mathcal{P}(I+S+D,t)

   f \sim \mathcal{N}(0,R^2)

   

(P((I+S+D)*t) + N(0,R))/t-S-D


   

Excluded Effects
^^^^^^^^^^^^^^^^


.. rubric:: Footnotes
.. [#f1] Currently the sky background is assumed as a single constant
	 value, and adding in the `master-sky backgrounds
	 <background.rst>`_ are not yet implemented.
   
