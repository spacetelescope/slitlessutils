.. _faq:

Frequently-Asked Questions
==========================


.. contents::
    :local:



I got this error/warning message...
-----------------------------------


**WARNING**: LOCAL JWST PRD VERSION PRDOPSSOC-059 CANNOT BE CHECKED AGAINST ONLINE VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is coming from `pysiaf <https://github.com/spacetelescope/pysiaf>`_ and indicates that either (1) you are not connected to the internet or (2) there was a change in ``pysiaf`` codebase.  You should verify your connectivity and/or upgrade the ``pysiaf`` with ``pip install pysiaf --upgrade``.  (You may find different three digit numbers instead of ``059``)



RuntimeError related to ``freeze_support()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specifically, the error is:

    RuntimeError: 
    An attempt has been made to start a new process before the
    current process has finished its bootstrapping phase.

    This probably means that you are not using fork to start your
    child processes and you have forgotten to use the proper idiom
    in the main module:
    
    if __name__ == '__main__':
        freeze_support()


    The "freeze_support()" line can be omitted if the program
    is not going to be frozen to produce an executable.


This is because you are running some script in the ``__main__`` of the file.  Try putting the commands inside a function declaration:

.. code: python
    
    def my_commands_here():  # doctest: +SKIP
        pass                 # doctest: +SKIP


    if __name__ == '__main__':  # doctest: +SKIP
        my_commands_here()      # doctest: +SKIP




Why do I....
------------

get a the same logging message printed multiple times?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is likely because there are multiple instances of the logger running for the separate threads.

