.. _simulation:

WFSS Image Simulation (`~slitlessutils.modules.Simulate()`)
===========================================================


Methodology
-----------


The current implementation will *only* simulate the signal from the sources, but includes the Poisson noise from the source(s), sky background, and dark current and Gaussian noise from the read noise.  To create a noiseless WFSS image:


#. Load :doc:`WFSS data<wfss>` and :doc:`spectral sources<sources>` as ``WFSSCollection()`` and ``SourceCollection()``, respectively.
#. Tabulate each WFSS image with the :doc:`tabulation module <tabulation>`
#. For each WFSS image:

  a. Initialize noiseless science as all zero: :math:`\tilde{S}_{x,y}=0` for all WFSS image pixels :math:`(x,y)`.
  b. For each detector in the WFSS file:  
    
    * For each source in the source collection:
      
      * For each :term:`direct imaging` pixel :math:`(x_d,y_d)` in the source:
    
        * load the PDT from the :class:`~slitlessutils.tables.PDTFile()`
        * append to a list

      * multiply the fractional pixel :math:`a_{(x_d,y_d)\rightarrow(x,y)}` area between the direct and WFSS image (and tabulated in the PDTs), :doc:`wavelength-dependent flat-field <calib>` :math:`F_{x,y}(\lambda)`, :doc:`sensitivity curve <calib>` :math:`T(\lambda)`, :term:`pixel-area map` :math:`P_{x,y}` (:ref:`see more below <pam>`), and the source spectrum :math:`f_{x_d,y_d}(\lambda)` associated with this direct-image pixel:

      .. math::

        s_{x,y,l} = a_{(x_d,y_d)\rightarrow(x,y)}\,\frac{I_{x_d,y_d}}{\sum\limits_{\mathbb{S}_i} I_{x_d,y_d}}\,F_{x,y}(\lambda)\,T(\lambda)\,f_{x_d,y_d}(\lambda)\, P_{x,y}\,\delta\lambda

      where :math:`I_{x_d,y_d}` is the direct-image brightness, :math:`\mathbb{S}_i` is the collection of direct-imaging pixels associated with this source (see :doc:`source description <sources>`), and :math:`\delta\lambda` is the tabulation bandwidth.

      * :term:`decimate<decimation>` the list of PDTs over the wavelength index (:math:`l`)
      * sum this decimated list into the noiseless science image:

      .. math::

        \tilde{S}_{x,y} \rightarrow \tilde{S}_{x,y} + s_{x,y}

To add noise to this noiseless image, ``slitlessutils`` requires two additional user-specified parameters :numref:`usertab`, which must be set when loading the `~slitlessutils.wfss.WFSSCollection()`.  See the discussion for the :doc:`WFSS data <wfss>` for more details.

.. _usertab:
.. list-table:: User-Specified Simulation Parameters
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword
     - Unit
     - Description
   * - Background
     - :math:`\mathrm{e}^-/\mathrm{s}`
     - | The notional background level, which is assumed to be constant across the
       | detector
   * - Exposure time
     - :math:`s`
     - The exposure time


Now the expected total number of electrons will be :math:`E = \left(\tilde{S}_{x,y}+B+D\right)\times t`, which is used to draw random Poisson and normal deviates for each pixel:

.. math::

  \begin{eqnarray}
    p_{x,y} &\sim& \mathcal{P}\left(E\right)\\
    g_{x,y} &\sim& \mathcal{N}\left(0,R^2\right)
  \end{eqnarray}

where :math:`B` and :math:`t` are the background and exposure time described in :numref:`usertab`, and :math:`D` and :math:`R` are the dark current and read noise (see the :doc:`instrument tables <instrumentfiles>`). Now the final, noised-science :math:`S_{x,y}` and uncertainty :math:`U_{x,y}` rate images are then:

.. math::

  \begin{eqnarray}
    S_{x,y} &=& \frac{p_{x,y}+g_{x,y}}{t} - B - D\\
    U_{x,y} &=& \sqrt{\left(\frac{\tilde{S}_{x,y} + B + D}{t}\right)+ \left(\frac{R}{t}\right)^2}
  \end{eqnarray}

both in units of :math:`\mathrm{e}^-/\mathrm{s}`.  The background rate and dark current are subtracted here to produce an image equivalent to a :doc:`sky-subtracted WFSS image <background>`.

.. important::
  Some detectors (e.g. WFC3/IR) record the images in :math:`\mathrm{e}^-/\mathrm{s}`, while others (e.g WFC3/UVIS, ACS/WFC, ACS/SBC) use :math:`e^-`.  This information is encoded in the instrument-specific ``yaml`` files in :file:`$HOME/.slitlessutils/<VERSION>/instruments` (but see also :doc:`instrument tables <instrumentfiles>`), which will modify the definitions for the final, noised images :math:`S_{x,y}` and :math:`U_{x,y}`.



Example
^^^^^^^

A full simulation is fairly involved, so please see the function :file:`slitlessutils.examples.starfield.py`, specifically the functions :func:`slitlessutils.examples.starfield.simulate_grisms()` and :func:`slitlessutils.examples.starfield.make_scene()`.


Additional Instrument Settings
------------------------------
There are many other parameters required to simulate a WFSS image, and these are stored in several ``yaml`` files in the configuration directory in :file:`{$HOME}/.slitlessutils/<VERSION>/instruments/`. Modifying these files is highly discouraged.

.. toctree::
  :titlesonly:
  :maxdepth: 1

  instrumentfiles.rst



Excluded Effects
----------------

The simulations provided by ``slitlessutils`` make several simplifying assumptions that will be reevaluated in future releases.  In order of relative importance of their adverse effect on the expected :term:`signal-to-noise` (S/N), these are:

* The sky background is assumed to be a single value, however as discussed in :doc:`the master sky <background>` belies this assumption.  Employing a realistic :term:`master-sky image` with a scale factor (:math:`\alpha`) by modifying the source/uncertainty equations to have :math:`B\rightarrow \alpha\,B_{x,y}`, where :math:`\alpha` becomes the tunable parameter to be set by the user.  The current implementation will give the illusion of a constant S/N over the detector, but the deviations from constant will depend on the how adopted level compares to the (large-scale) variations in the :term:`master-sky image`. Therefore this may introduce small systematic biases based on the position of the sources.

* The :term:`DQA<data-quality array>` is assumed to have no bad pixels flagged, which effectively *overestimates* the number of valid science pixels and perhaps slightly the S/N.

* The dark current and read noise are assumed to be a single values that applies uniformly to *all* pixels, yet real detectors may have pixel-to-pixel variations.  Like the sky-background issue, this may introduce weak systematic, spatial biases.

* The :term:`telescope attitude<attitude>` is set by the user and assumed to be noiseless, but in practice there are systematic uncertainties in the accuracy of the :term:`world-coordinate system` (WCS).  In general, errors in the WCS result in a systematic wavelength shift (sometimes called the *wavelength zeropoint*) and/or flux losses (akin to slit-losses in aperture spectroscopy).  However `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ show that these effects are very small compared for most HST observations and negligible compared to the spectro-photometric noise.


.. _pam:

An Aside on Pixel-Area Maps
---------------------------

The :term:`pixel-area map` (PAM) describes the relative pixel size due to distortions in the detector, which is given by the absolute value of the determinant of the Jacobian matrix.  In principle, the distortion can be described in many ways (e.g. look-up table), but ``slitlessutils`` currently assumes this will be described by `Simple-Imaging Polynomials (SIP) <https://docs.astropy.org/en/stable/wcs/note_sip.html>`_.  In which case, the Jacobian is simply:

.. math::

    J = \left(\begin{array}{cc} \partial a/\partial x & \partial a/\partial y\\
      \partial b/\partial x & \partial b/\partial y\end{array}\right)

where all of these partial derivatives are (by definition) polynomials of :math:`(x,y)`.  Therefore, the PAM becomes:

.. math::

  \begin{eqnarray}
    P_{x,y} &=& \left|\det(J)\right|\\
            &=& \left|\frac{\partial a}{\partial x}\frac{\partial b}{\partial y} - \frac{\partial b}{\partial x}\frac{\partial a}{\partial y}\right|.
  \end{eqnarray}

``Slitlessutils`` computes these values on-the-fly as needed.
