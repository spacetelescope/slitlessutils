.. _cosmicrays:

Cosmic Rays
===========

Introduction
------------

Cosmic rays (CRs) are high energy particles that may often impart a
large amount of charge that is unpredictable and has a
characteristically sharp spatial profile.  For the multiaccum infrared
detectors, the up-the-ramp sampling can flag the overwhelming majority
of these events and populate the data-quality arrays (DQAs), the CCDs
do not have any analogous flagging algorithms.  In this case, there
are effectively two options for identifying pixels that

#. deviate by some threshold between independent exposures; or

#. have the sharp spatial profile.

The former method is the canonical approach for standard imaging, but
for grism data, this approach gets additionally complicated. At
present, `slitlessutils` has tools to flag pixels based on the profile
shapes.


Edge Detection with Laplacians
------------------------------

Detecting sharp edges in imaging has long been a subject of computer
vision research, and one common technique is to identify regions where
the flux distribution (:math:`I`) changes concavity.  This can be
achieved by finding where the second derivative goes to zero, and the
second derivative is given as the Laplacian:

.. math::
   \nabla^2 I = \frac{\partial^2 I}{\partial x^2}+\frac{\partial^2 I}{\partial y^2}.

For a pixelated light distribution, the Laplacian must be extended for
finite differences, which is an approximation to the continuous case.
One such approximation is given by:

.. math::
   \nabla^2 I_{x,y} \approx \frac{I_{x+h,y}+I_{x-h,y}-4\,I_{x,y}+I_{x,y-h}+I_{x,y+h}}{h^2}

and for :math:`h=1`, this expression is concisely given as a simple image
convolution :math:`\nabla^2 I \approx K \ast I`.  `slitlessutils` offers
several forms for the Laplacian covolution kernel:

.. math::
   
   K_{3a} = \left(\begin{array}{rrr}  0 & -1 &  0 \\
   -1 & +4 & -1 \\
    0 & -1 &  0 \end{array}\right)

   K_{3b} = \left(\begin{array}{rrr} -1 & -1 & -1 \\
   -1 & +8 & -1 \\
   -1 & -1 &  -1 \end{array}\right)


   K_{3c} = \left(\begin{array}{rrr} +1 & -2 & +1 \\
   -2 & +4 & -2 \\
   +1 & -2 & +1 \end{array}\right)

   K_{5a} = \left(\begin{array}{rrrrr}  0 &  0 & -1 &  0 &  0 \\
    0 & -1 & -2 & -1 &  0 \\
   -1 & -2 & +16 & -2 & -1 \\
    0 & -1 & -2 & -1 &  0 \\
    0 &  0 & -1 &  0 &  0 \end{array}\right)

where :math:`K_{3a}` is the kernel for the approximation [#f1]_.
After convolving the image with the Laplacian kernel, pixels that
deviate more than :math:`n` times above their respective uncertainties
(:math:`U`) are considered as candidate CR pixels:

.. math::
   \left|\nabla^2 I\right| \geq n\, U

These candidate pixels are grouped based on their connectivity (see
`skimage.measure.label()`) and only groups with a minimum number of
pixels are kept.  Finally, the remaining groups can be grown using
standard dilation operations (see `skimage.morphology.dilation()`),
and several different footprints (`square`, `rectangle`, `diamond`,
`disk`, `octagon`, and `star` --- see the respective functions in
`skimage.morphology`).



Example
~~~~~~~
    
This are the kernels and can be envoked by the subscript, for example

.. code:: python
   	  
   import slitlessutils as su

   # not totally necessary, but this will engage the slitlessutils logger
   su.start_logging()

   # perform the master sky subtraction on the filename "grismfile"
   su.core.preprocess.crrej.laplace(grismfile, kernel='3a', inplace=True)

This will update the file in place, as the flag is set: :code:`inplace=True`.  See :numref:`animatedcrs` for an animation of how cosmic rays appear and then can be bilinearly-interpolated over.


.. _animatedcrs:
.. figure:: images/cr_animation.gif
   :width: 600
   :alt: Example for cosmic ray flagging and interpolation from convolution from a Laplacian kernel.

   Example of cosmic-ray flagging from convolution from a Laplacian kernel and bilinear 
   interpolation to highlight the differences.

      
AstroDrizzle Cosmic Ray flagging
--------------------------------

Coming Soon.





Examples
~~~~~~~~


      

    
    
.. rubric:: Footnotes
.. [#f1] It is worth mentioning that Laplacian kernels must share the
	 property that :math:`\sum_{i,j}K_{i,j}=0`.
