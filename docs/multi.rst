.. _multi:



Multi-Exposure Extraction (`~slitlessutils.modules.extract.multi`)
==================================================================

`Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_

.. math::
      \chi^2 = ||I_{\vartheta} - W_{\vartheta,\varphi} f_{\varphi}||^2   

.. math::
      \psi^2 = \chi^2 + \ell\,||W||_F^2||f-f_0||^2

where :math:`\ell` is the regularization parameter


.. math::
      ||W||_F^2 = \sum_i\sum_j w_{i,j}^2


Regularization Optimization
---------------------------




.. _lcurveexample:
.. figure:: images/starfield_multi_lcv.pdf
   :align: center
   :alt: Example regularization plot.

   The top panel shows the standard L-curve with the scaling factor of the `Frobenius norm <https://en.wikipedia.org/wiki/Matrix_norm>`_ to ensure that the regularization parameter :math:`\ell` is dimensionless, which is encoded in the color of the plot symbols (see colorbar at the very bottom).  The lower panel shows the `Menger curvature <https://en.wikipedia.org/wiki/Menger_curvature>`_ as a function of the logarithm (base 10) of the (dimensionless) regularization parameter.  The clear peak at :math:`\log\ell\sim-1.9` represents the sharp vertex in the L-curve at :math:`\sim(2.1,3.6)`.  This point is adopted as it represents a roughly "equal" trade-off between modeling the data (ie. the parameter on the x-axis) and damping high-frequency structure (ie. the parameter on the y-axis).






