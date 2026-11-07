2.1 (unreleased)
================

New Features
------------

Bug Fixes
---------

Other Changes and Additions
---------------------------


2.0.1 (unreleased)
==================

Bug Fixes
---------

- Updated to use scipy sparse arrays instead of sparse matrices, due to upcoming
  deprecation of the latter. [#158]

- Fixed ``WFSSCollection.from_wcsfile`` failing with pandas 3.1 when the
  ``blocking`` field is empty. [#167]

Other Changes and Additions
---------------------------

- Bumped the minimum required versions of ``drizzlepac`` to 3.11,
  ``numpy`` to 2.0.1, ``pandas`` to 2.2.2, ``psutil`` to 5.9.6,
  ``pypolyclip`` to 1.1.0, ``pysiaf`` to 0.24.0, ``scikit-image`` to
  0.23.1, ``scipy`` to 1.13, and ``shapely`` to 2.0.4. [#168]


2.0 (2026-03-27)
================

New Features
------------

- Improved/enhanced support for UVIS, fixed bugs in single-ended mode, and
  added utility functions. [#119]

- Improved forward modeling performance. [#123]

- Improved matrix handling and implemented grouping of spectral traces
  for linear reconstruction in multi-orient extraction. [#136]

- Implemented support for subarray embedding. [#110, #130]

Bug Fixes
---------

- Fixed contamination bugs and bug in ``Band`` class. [#141]

- Fixed global background retrieval for subarrays. [#86]

Other Changes and Additions
---------------------------

- Added support for NumPy 2.0. [#70]

- Updated and improved documentation. [#133, #139, 142]

- Bumped python minimum version to 3.12. [#143]

- Added more example scripts. [#136]

1.0 (2023-11-28)
================

New Features
------------

- Added Zenodo link to the README and CITATION file.


0.9 (2023-11-28)
================

New Features
------------

- Initial implementation of the package. See documentation for features.
