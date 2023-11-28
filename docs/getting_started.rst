.. _gettingstarted:

Getting Started with ``slitlessutils``
======================================

The general philosophy of ``slitlessutils`` is to instantiate two :doc:`primary data structures <datastructs>` that emulate a ``dict`` and contain the collection of :doc:`sources <sources>` and :doc:`spectroscopic images <wfss>`.  These structures are generally passed to the various :doc:`models and preprocessing steps <modules>` to carry out the primary goals (e.g. simulation, extraction, etc.).