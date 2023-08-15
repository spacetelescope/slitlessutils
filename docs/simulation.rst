.. _simulation:

WFSS Image Simulation
=====================


Methodology
-----------


The current implementation will *only* simulate the signal from the sources, while including noise from many effects.  To simulate an image



#. Tabulate the WFSS image with the :doc:`tabulation module <tabulation>`
#. Initialize noiseless science and uncertainty images


.. list-table:: Simulation Parameters
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Background 
     - :math:`e^-/s`
     - The notional background level, which is assumed to be constant across the detector
   * - Exposure time
     - :math:`s`
     - The exposure time





.. math::
   :name: sneqns

   \begin{eqnarray}
      S' &\sim& \mathcal{P}left(t\,(S+B+D)\right)/t - B - D + \mathcal{N}(0,R^2)\\
      U &=& \frac{\sqrt{(I+B+D) t+R^2}}{t} 
   \end{eqnarray}

 The science image(s) is en


:doc:`tabulation module <tabulation>`

The uncertainty image is given by the

.. math::
   U = \frac{\sqrt{(I+B+D) t+R^2}}{t}

where :math:`I`, :math:`B` [#f1]_, and :math:`D` are the Poissonian noise terms that represent the flux (in :math:`e^-`/s) from the simulated science image, the specified background level, and the dark rate, respectively.  The read noise (in :math:`e^-`) is specified as :math:`R`, and represents the lone Gaussian noise term.  The specified exposure time (in s) is given by :math:`t`.  Therefore, the simulated images will have an `ERR` extension will be populated with these values.

The `SCI` extension

.. math::
   p \sim \mathcal{P}(I+S+D)

   f \sim \mathcal{N}(0,R^2)

   
.. note::
   The WFC3/IR images are in units of :math:`e-`/s, while all the data for all other instruments will be in :math:`e-`.  



   

Excluded Effects
^^^^^^^^^^^^^^^^

The simulations provided by ``slitlessutils`` make several simplifying assumptions that will be reevaluated in future releases.  In order of relative importance of their adverse effect on the expected :term:`signal-to-noise` (S/N), these are:

* The sky background is assumed to be a single value, however as discussed in :doc:`the master sky <background>` belies this assumption.  Employing a realistic :term:`master-sky image` with a scale factor (:math:`\alpha`) is simply modifying :ref:`<sneqns>` to have :math:`B\rightarrow \alpha\,B_{x,y}`.  This assumption will give the illusion of a constant S/N over the detector, but the deviations from constant will depend on the how adopted level compares to the (large-scale) variations in the :term:`master-sky image`. Therefore this may introduce small systematic biases based on the position of the sources.

* The DQA is assumed to have no bad pixels flagged, which effectively *overestimates* the number of valid science pixels and perhaps slightly the S/N.

* The dark current is assumed to be a single value that applies uniformly to *all* pixels, yet real detectors have pixel-to-pixel variations.  Like the sky-background issue, this may introduce weak systematic, spatial biases.

* The :term:`attitude` is set by the user and assumed to be noiseless, but in practice there are systematic uncertainties in the accuracy of the :term:`world-coordinate system` (WCS).  In general, errors in the WCS result in a systematic wavelength shift (sometimes called the *wavelength zeropoint*) and/or flux losses.  However `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ show that these effects are very small compared for most HST observations and negligible compared to the spectro-photometric noise.  


.. rubric:: Footnotes
.. [#f1] Currently the sky background is assumed as a single constant value, and adding in the :doc:`master-sky backgrounds <background>` are not yet implemented.
   
