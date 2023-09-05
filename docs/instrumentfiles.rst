.. _instrumentfiles:

Instrument Files :file:`$HOME/.slitlessutils/instruments/`
==========================================================

Most of these parameters are the subject of considerable calibration efforts, and as such, should probably not be adjusted if the results are to be trusted.  


Instrument-Wide Settings
------------------------

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
     - The file suffix in the HST parlance for without and with the CTE corrections.
   * - Path
     - ``str``
     - | The relative path from the ``yaml`` file where the files for this 
       | instrument are stored.
   * - Focal-plane position
     - 3-elements
     - | The :math:`(v_2,v_3)` position of the reference point 
       | and :math:`v_{3y}` angle with respect the :math:`v_3`-axis.

Instrument-Wide Grating/Blocking Parameters
-------------------------------------------


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

Detector Parameters
-------------------

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
     - | dark current :math:`D` in :math:`e^-/s`
       | readnoise :math:`R` in :math:`e^-`
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


.. rubric:: Footnotes
.. [#gbnote] These settings are set for each grating/blocking combination, and if no blocking filter exists, then it is set as the ``null`` variable in ``yaml``.
.. [#extnote] The extraction and tabulation settings need-not be the same.  Indeed, to encapsulate the non-linearity in the prism modes they will **NOT** be the same.
.. [#detnote] There should be a separate stanza like this for each detector in the instrument (e.g. such as the two CCDs in ACS-WFC).

