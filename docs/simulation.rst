.. _simulation:

WFSS Image Simulation
=====================


Methodology
-----------


The current implementation will *only* simulate the signal from the sources, but includes the Poisson noise from the source(s), sky background, and dark current and Gaussian noise from the read noise.  To create a noiseless WFSS image:


#. Load :doc:`WFSSCollection <wfss>` and :doc:`sources <sources>`, 
#. Tabulate the WFSS image with the :doc:`tabulation module <tabulation>`
#. Initialize noiseless science as all zero: :math:`S=0`
#. For each detector in the WFSS file:
      > For each source in the source collection:
         + For each :term:`direct imaging` pixel :math:`(x_d,y_d)` in the source:
            * load the PDT from the :class:`~slitlessutils.tables.PDTFile()`
            * append to a list
         + multiply the fractional pixel (:math:`a`) from the PDTs by wavelength-dependent flat-field (:math:`F(\lambda)`), :term:`sensitivity curve` (:math:`\mathcal{S}(\lambda)`, :term:`pixel-area map` (:math:`A`, see more below), and the source spectrum (:math:`f(\lambda)`)

         .. math::
            s = a\,F(x,y,\lambda)\,\mathcal{S}(\lambda)\,f(\lambda)


         + decimate the list of PDTs over the WFSS image pixels :math:`(x,y)`.  



.. math::
   
   J = \left(\begin{array}{cc} \partial a/\partial x & \partial a/\partial y\\
      \partial b/\partial y & \partial b/\partial y\end{array}\right)

.. math::

   \abs\left(\det(J)\right) = \left|\frac{\partial a}{\partial x}\frac{\partial b}{\partial y} - \frac{\partial a}{\partial y}\frac{\partial b}{\partial x}\right|




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


.. list-table:: Instrument-Wide Settings
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Image units
     - :math:`e^-/s` or :math:`e^-`
     - The units of the images to be written.
   * - File suffix
     - ``flt`` or ``flc``
     - The file suffix in the HST parlance.
   * - Path
     - ``str``
     - | The relative path from the ``yaml`` file where the files for this 
       | instrument are stored.
   * - Focal-plane position
     - 3-elements
     - | The :math:`(v_2,v_3)` position of the reference point 
       | and :math:`v_{3y}` angle with respect the :math:`v_3`-axis.

.. list-table:: Instrument-Wide Grating/Blocking [#gbnote]_
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Master-Sky Image
     - 
     - The name of the master-sky image.
   * - Tabulation Parameters
     - ``dict``
     - | This contains the starting wavelength (``wave0``), ending 
       | wavelength (``wave1``), sampling frequency (``dwave``), 
       | units (usually ``angstrom``), and disptype.  
   * - Extraction Parameters [#extnote]_
     - ``dict``
     - | This contains the starting wavelength (``wave0``), ending 
       | wavelength (``wave1``), sampling frequency (``dwave``), 
       | units (usually ``angstrom``), and disptype.  


.. list-table:: Detector Settings [#detnote]_
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Focal-plane position
     - 3-elements
     - | :math:`(v_2,v_3)` position of the reference point 
       | :math:`v_{3y}` angle with respect the :math:`v_3`-axis
   * - Extension properties
     - 
     - | ``name``: the name of the extension (must be ``str``)
       | ``ver``: the version of the extension (must be ``int``)
       | ``dtype``: a valid ``np.dtype``
   * - Noise properties
     - 
     - | dark current in :math:`e^-/s`
       | readnoise in :math:`e^-`
   * - Detector dimensionality
     - 
     - | ``naxis``: 2-element list of size of detector (must be ``int``)
       | ``crpix``: 2-element list for reference position (can be ``float``)
       | ``scale``: 2-element list for pixel scale (can be ``float``)
   * - Distortion model
     - 
     - `SIP coefficients <https://docs.astropy.org/en/stable/wcs/note_sip.html>`_ should be a ``dict``
   * - Configuration files
     - 
     - The file name for each grating/blocking combination



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


.. rubric:: Footnotes
.. [#gbnote] These settings are set for each grating/blocking combination, and if no blocking filter exists, then it is set as the ``null`` variable in ``yaml``.
.. [#extnote] The extraction and tabulation settings need-not be the same.  Indeed, to encapsulate the non-linearity in the prism modes they will **NOT** be the same.
.. [#detnote] There should be a separate stanza like this for each detector in the instrument (e.g. such as the two CCDs in ACS-WFC).
