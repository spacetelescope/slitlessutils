.. _modules:


Computational Modules (`~slitlessutils.modules.Module()`)
=========================================================

`Slitlessutils` has several primary computational modules that carry out the most important calculations.  Each module will have the same mandatory arguments, but a potential host of keyword-arguments related to the unique nature of each module.  Since the calculations of each module can be potentially very computationally intensive and many of them are `embarrassingly parallel <https://en.wikipedia.org/wiki/Embarrassingly_parallel>`_, they are outfitted with `multi-threading capabilities <https://docs.python.org/3/library/multiprocessing.html>`_.  This is implemented in the parent class :class:`~slitlessutils.modules.Module()` that is inherited by each module. 

It is unlikely that a user should ever need instantiate the parent `Module()` class directly, but there are three key keyword-arguments that control the multiprocessing and temporary tables (see :doc:`Tabulating Module <tabulation>`) that one might need to pass to the parent via any of the child Modules:

.. _modulekeys:

.. list-table:: Module Keyword Arguments
  :widths: 15 10 75
  :header-rows: 1
  :stub-columns: 0
  :width: 600

  * - Keyword
    - Datatype
    - Notes
  * - ncpu
    - ``int`` or ``None``
    - the number of cpu threads to use.  If set to ``None``, then the total number on the system **minus one** will be used. Otherwise, use the integer as specified.  Default is ``None``.
  * - tables
    - ``str``
    - local path with respect to the current-working directory where the temporary tables are stored.




Primary Computational Modules
-----------------------------
.. toctree::
	:maxdepth: 1

	tabulation.rst
	extraction.rst
	simulation.rst
	grouping.rst
	regions.rst
	
