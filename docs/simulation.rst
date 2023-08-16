.. _simulation:

WFSS Image Simulation (`~slitlessutils.modules.Simulate()`)
===========================================================


Methodology
-----------


The current implementation will *only* simulate the signal from the sources, but includes the Poisson noise from the source(s), sky background, and dark current and Gaussian noise from the read noise.  To create a noiseless WFSS image:


#. Load :doc:`WFSSCollection <wfss>` and :doc:`sources <sources>`, 
#. Tabulate the WFSS image with the :doc:`tabulation module <tabulation>`
#. Initialize noiseless science as all zero: :math:`S_{x,y}=0` for all WFSS image pixels :math:`(x,y)`.
#. For each detector in the WFSS file:
      > For each source in the source collection:
         + For each :term:`direct imaging` pixel :math:`(x_d,y_d)` in the source:
            * load the PDT from the :class:`~slitlessutils.tables.PDTFile()`
            * append to a list
         + multiply the fractional pixel :math:`a_{x,y}` from the PDTs by :doc:`wavelength-dependent flat-field <calib>` :math:`F_{x,y}(\lambda)`, :doc:`sensitivity curve <calib#sensitivity>` :math:`s(\lambda)`, :term:`pixel-area map` :math:`A_{x,y}` (:ref:`see more below <pam>`), and the source spectrum :math:`f(\lambda)`

         .. math::
            s_{x,y,l} = a_{x,y}\,F_{x,y}(\lambda)\,\mathcal{S}(\lambda)\,f(\lambda)\, A_{x,y}

         + decimate the list of PDTs over the wavelength index (:math:`l`)
         + sum this decimated list into to noiseless science image:

         .. math::
            S_{x,y} \rightarrow S_{x,y} + s_{x,y}





.. _pam:
Pixel-Area Maps
~~~~~~~~~~~~~~~

The :term:`pixel-area map` (PAM) describes the relative pixel size due to distortions in the detector, which is given by the absolute value of the determinant of the Jacobian:
.. math::
   
    J = \left(\begin{array}{cc} \partial a/\partial x & \partial a/\partial y\\
      \partial b/\partial y & \partial b/\partial y\end{array}\right)




.. list-table:: User-Specified Simulation Parameters
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Background 
     - :math:`e^-/s`
     - | The notional background level, which is assumed to be constant across the 
       | detector
   * - Exposure time
     - :math:`s`
     - The exposure time



However, there are many other parameters required to simulate a WFSS image, and these are stored in ``yaml`` files in the configuration directory in :file:`{$HOME}/.slitlessutils`.  Most of these parameters are the subject of considerable calibration efforts, and as such, should probably not be adjusted if the results are to be trusted.

.. toctree::
  :titlesonly:
  :maxdepth: 1

  instrumentfiles.rst 




.. math::

   \begin{eqnarray}
      S' &\sim& \mathcal{P}\left(t\,(S+B+D)\right)/t - B - D + \mathcal{N}\left(0,R^2\right)\\
      U &=& \frac{\sqrt{(S+B+D) t+R^2}}{t} 
   \end{eqnarray}


 The science image(s) is en


:doc:`tabulation module <tabulation>`

The uncertainty image is given by the

.. math::
   U = \frac{\sqrt{(I+B+D)\,t+R^2}}{t}

where :math:`I`, :math:`B`,  and :math:`D` are the Poissonian noise terms that represent the flux (in :math:`e^-`/s) from the simulated science image, the specified background level, and the dark rate, respectively.  The read noise (in :math:`e^-`) is specified as :math:`R`, and represents the lone Gaussian noise term.  The specified exposure time (in s) is given by :math:`t`.  Therefore, the simulated images will have an `ERR` extension will be populated with these values.

The `SCI` extension

.. math::
   p \sim \mathcal{P}(I+S+D)

   f \sim \mathcal{N}(0,R^2)

   
.. note::
   The WFC3/IR images are in units of :math:`e-`/s, while all the data for all other instruments will be in :math:`e-`.  



   

Excluded Effects
^^^^^^^^^^^^^^^^

The simulations provided by ``slitlessutils`` make several simplifying assumptions that will be reevaluated in future releases.  In order of relative importance of their adverse effect on the expected :term:`signal-to-noise` (S/N), these are:

* The sky background is assumed to be a single value, however as discussed in :doc:`the master sky <background>` belies this assumption.  Employing a realistic :term:`master-sky image` with a scale factor (:math:`\alpha`) by modifying the source/uncertainty equations to have :math:`B\rightarrow \alpha\,B_{x,y}`.  This assumption will give the illusion of a constant S/N over the detector, but the deviations from constant will depend on the how adopted level compares to the (large-scale) variations in the :term:`master-sky image`. Therefore this may introduce small systematic biases based on the position of the sources.

* The DQA is assumed to have no bad pixels flagged, which effectively *overestimates* the number of valid science pixels and perhaps slightly the S/N.

* The dark current is assumed to be a single value that applies uniformly to *all* pixels, yet real detectors have pixel-to-pixel variations.  Like the sky-background issue, this may introduce weak systematic, spatial biases.

* The :term:`attitude` is set by the user and assumed to be noiseless, but in practice there are systematic uncertainties in the accuracy of the :term:`world-coordinate system` (WCS).  In general, errors in the WCS result in a systematic wavelength shift (sometimes called the *wavelength zeropoint*) and/or flux losses.  However `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ show that these effects are very small compared for most HST observations and negligible compared to the spectro-photometric noise.  


