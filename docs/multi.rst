.. _multi:


Multi-Orient Extraction (`~slitlessutils.modules.extract.multi`)
================================================================

The multi-orient extraction was first developed by `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ and referred to as *LINEAR*.  This method was founded to address the contamination (or confusion) from overlapping spectral traces, where the detector records a weighted sum of the the spectra (see :numref:`confusion` below).  Therefore, this represents a fundamental degeneracy, which can only be broken with additional data.  The :doc:`single-orient extraction <single>` formalism uses the broadband data to inform the contamination model, whereas the *LINEAR* uses data at multiple orients provide a self-consistent, spectroscopic model for the entire scene.  

.. _confusion:
.. figure:: images/confusion.png
   :align: center
   :alt: Illustration of confusion

   An illustration of spectral contamination or confusion (taken from `Ryan, Casertano, & Pirzkal 2018 <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_).  Along the left side of each panel, they show a red and blue point source whose celestial positions are fixed (as shown the coordinate vane).  This shows that under certain orients (lower panel), the spectra will overlap leading to the spectral degeneracy.  However, if the telescope is reoriented, then the spectral traces separate, which provides the leverage to break this degeneracy.


Although the formative effort was to break the degeneracy from contamination/confusion, additional advantages were identified. Namely, the spectral resolution of a WFSS mode, is set by properties of the grating and the size of the source.  In analogy to the relationship between slit width and spectral resolution for long-slit spectroscopy, the size of the source projected along the dispersion axis set the resolution --- where the larger the source, the lower the resolution.  Therefore, simply averaging spectra extracted at separate orients (such as described in :doc:`single-orient extraction <single>`) will result in biases, where the resolutions are different.  Hence, one can either smooth all the data to a common resolution before averaging (which loses spectral resolution) or analyze the high-resolution data by itself (which loses signal-to-noise).  But `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ find that the *LINEAR* method recovers the highest spectral resolution present in the data, at the cost of an increased computational resources and collecting data at multiple orients. 


Mathematical Foundation
-----------------------

The *LINEAR* framework acknowledges that the flux in a WFSS image pixel is the sum of all sources and wavelengths, weighted by factors related to the sources (e.g. the cross-dispersion profiles) and the detector (e.g. :doc:`sensitivity curve, flat-field <calib>`, or :doc:'pixel-area map <simulation>`):

.. math::

   S_{x,y,i} = \sum_{l}\sum{k=1}^{N_\mathrm{obj}} W_{x,y,i,l,k} f_{l,k}



The known and unknown indices are grouped together with `np.ravel_multi_index() <https://numpy.org/doc/stable/reference/generated/numpy.ravel_multi_index.html>`_ as:

.. math::

   \begin{eqnarray}
      \vartheta &=& x + n_x\,y+ n_x\,n_y\,i\\
      \varphi &=& l + n_l\,k
   \end{eqnarray}

Now the above matrix-equation is recast as:

.. math::
   
   S_{\vartheta} = \sum_\varphi W_{\vartheta,\varphi} f_{\varphi}.

But since this is an overconstrained problem, then the vector of unknowns :math:`f_{\varphi}` must be solved with optimization techniques:

.. math::

   \chi^2 = \sum_{\vartheta} \left(\frac{S_{\vartheta} - \sum_{\varphi} W_{\vartheta,\varphi}\,f_{\varphi}}{U_{\vartheta}}\right)^2

which is simplified by subsuming the WFSS uncertainties into the data and linear operator as:

.. math::

   \begin{eqnarray}
      S_{\vartheta} &\rightarrow& \frac{I_{\vartheta}}{U_{\vartheta}}\\
      W_{\vartheta,\varphi} &\rightarrow& \frac{W_{\vartheta,\varphi}}{U_{\vartheta}}
   \end{eqnarray}

so that now :math:`\chi^2 = ||I - W\,f||^2`.  Although this can be directly solved, the poor condition number of :math:`W` can amplify the input noise into the output result, which can be ameliorated by including a `regularization term <https://en.wikipedia.org/wiki/Ridge_regression>`_.  Additionally, for most WFSS observations, the linear operator :math:`W` will be extremely sparse, which permits specialized techniques to iteratively compute the unknown vector :math:`f_{\varphi}` without computing the pseudo-inverse of :math:`W`.  However, `Ryan, Casertano, & Pirzkal (2018) <https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract>`_ re-frame the regularization term so that the regularization parameter becomes dimensionless.

.. math::

   \psi^2 = \chi^2 + \ell\,\xi^2

where :math:`\ell` is the regularization parameter and :math:`\xi^2 = ||W||_F^2\,\sum\left(f_{\varphi}-f_{\varphi}^{(0)}\right)^2` with :math:`||W||_F` is the `Frobenius norm <https://en.wikipedia.org/wiki/Matrix_norm>`_ and :math:`f_{\varphi}^{(0)}` is the :term:`damping target`, which is initialized from the broadband data. 




.. _matrix:

Sparse Linear-Operator Construction
-----------------------------------




.. _solutions:

Sparse Least-Squares Solution
-----------------------------

There have been several algorithms devised to find the vector :math:`f_{\varph}` that minimizes the cost function for :math:`\psi^2`, and many have been implemented into the `scipy sparse solvers <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#module-scipy.sparse.linalg>`_.  However, ``slitlessutils`` is only organized to work with the two most common methods:

* LSQR: first presented by `Paige & Saunders (1982) <https://dl.acm.org/doi/10.1145/355984.355989>`_, is the standard tool for these types of linear systems. 
* LSMR: later developed by `Fong & Saunders (2011) <https://arxiv.org/abs/1006.0758>`_, and improves upon LSQR by generally converging faster.  


.. warning::
   Based on experimentation with the *LINEAR* work, the LSQR solver yields better results, and so it is set as the default sparse least-squares solver.


.. _regularization:

Regularization Optimization
---------------------------

As discussed above, the regularized least-squares introduces a tunable parameter that trades between modeling the data (ie. the :math:`\chi^2`-term) and damping the high frequency noise present in inverse problems (ie. the :math:`\xi^2`-term).  However, there have been heuristic approaches at "optimizing" the damping parameter :math:`\ell`, and the most common method is to consider a plot of :math:`\xi^2` versus :math:`\chi^2`, which often called the "L-curve" as when plotted as log-log, this will show a characteristic sharp resembling a capital-L (see :numref:`lcurveexample`).  It is widely accepted that the vertex of the L is represents a good compromise, and so there are several techinques to honing in on this critical point.  


#. Single-value: Accept a single value of the regularization parameter, and return the vector :math:`f_{\varphi}`.
#. Brute-force search: Define a linear grid of :\math:`\ell`, compute the Menger curvature at all points, and return the value of :math:`f_{\varphi}` that is associated with the maximizing value of :math:`\ell`.
#. Golden-ratio search: 

.. _lcurveexample:
.. figure:: images/starfield_multi_lcv.pdf
   :align: center
   :alt: Example regularization plot.

   The top panel shows the standard L-curve with the scaling factor of the Frobenius norm to ensure that the regularization parameter :math:`\ell` is dimensionless, which is encoded in the color of the plot symbols (see colorbar at the very bottom).  The lower panel shows the `Menger curvature <https://en.wikipedia.org/wiki/Menger_curvature>`_ as a function of the logarithm (base 10) of the (dimensionless) regularization parameter.  The clear peak at :math:`\log\ell\sim-1.9` represents the sharp vertex in the L-curve at :math:`\sim(2.1,3.6)`.  This point is adopted as it represents a roughly "equal" trade-off between modeling the data (ie. the parameter on the x-axis) and damping high-frequency structure (ie. the parameter on the y-axis).




