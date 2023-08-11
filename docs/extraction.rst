.. _extraction:

Spectral Extraction
===================


Introduction
------------
Spectral extraction refers to the conversion of two-dimensional spectroscopic image(s) to a one-dimensional, fluxed spectrum.  While there are many discussions of this process in the literature (e.g. `Horne (1986) <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`_), many of these focus on *aperture*-based spectrographs (e.g. long-slits or fibers).  Certainly, many of the aspects of the reduction/analysis of these data are applicable or have corollaries to slitless spectrographs, but there lack of the aperture leads to many complicating factors.  


Preprocessing WFSS Images
-------------------------

Before one can extract WFSS data, there are several *preprocessing* steps that should considered.  Below is a table that lists the three major preprocessing steps that one should implement for their data:






* `astrometry <astrometry>`
* `background <background>`
* `cosmic ray rejection <cosmicrays>`


How to extract the data?
------------------------


* `single <single>`
* `multi <multi>`
